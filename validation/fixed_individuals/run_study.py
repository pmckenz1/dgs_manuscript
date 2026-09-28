"""Simulate 60 populations, fit nested fixed samples, and collect the heatmap table."""
from pathlib import Path
from collections import Counter
from functools import lru_cache
import argparse
import hashlib
import json
import subprocess
import time

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from selfdgs import dgs_probabilities
from selfdgs._version import __version__ as model_version

MODEL = HERE / 'fixed_individuals.slim'
RATES = tuple(round(i * .05, 2) for i in range(20))
CENSUS_SIZES = (100, 250, 500)
SAMPLE_SIZES = (20, 50, 100)
CANDIDATES = tuple(float(s) for s in np.linspace(0, .95, 20))
TASKS = tuple((N, s) for N in CENSUS_SIZES for s in RATES)
BASE_SEED = 20260915


def seed_for(N, s, stream):
    """Return deterministic simulation (stream 0) and sampling (stream 1) seeds."""
    bits = np.float64(s).view(np.uint64).item()
    sequence = np.random.SeedSequence([
        BASE_SEED, N, bits & 0xffffffff, bits >> 32, 0, stream
    ])
    return int(sequence.generate_state(1)[0]) + 1


def parameters(N, s):
    return dict(census_size=N, selfing_rate=s, n_loci=10000, locus_length=1000,
                mu=1e-7, within_locus_r=5e-8, burn_mult=10, seed=seed_for(N, s, 0))


def write_json(path, value):
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def simulate(params, directory, slim):
    """Reuse compatible completed runs; otherwise run the local SLiM model."""
    directory.mkdir(parents=True, exist_ok=True)
    version = subprocess.run([slim, '-v'], capture_output=True, text=True, check=True).stdout.strip()
    if not version.startswith('SLiM version 5.2,'):
        raise ValueError('Use SLiM 5.2 to reproduce the saved study.')
    signature = dict(parameters=params, slim_version=version,
                     script_sha256=hashlib.sha256(MODEL.read_bytes()).hexdigest())
    manifest = directory / 'simulation.json'
    outputs = [directory / 'population.vcf', directory / 'histories.csv']
    if manifest.exists():
        if json.loads(manifest.read_text())['signature'] != signature:
            raise ValueError(f'Incompatible existing simulation: {directory}; use another --output.')
        if all(p.exists() and p.stat().st_size for p in outputs):
            return
    defines = dict(params, vcf_path=str(outputs[0].resolve()), history_path=str(outputs[1].resolve()))
    command = [slim]
    for name, value in defines.items():
        command.extend(['-d', f'{name}={json.dumps(value)}'])
    command.append(str(MODEL))
    start = time.monotonic()
    with (directory / 'slim.stdout').open('w') as out, (directory / 'slim.stderr').open('w') as err:
        subprocess.run(command, stdout=out, stderr=err, check=True)
    if not all(p.exists() and p.stat().st_size for p in outputs):
        raise RuntimeError(f'SLiM output missing in {directory}')
    write_json(manifest, dict(signature=signature, elapsed_seconds=time.monotonic()-start))


def read_population(directory, params):
    """Read complete biallelic diploid genotypes in VCF sample-column order."""
    N = params['census_size']
    records, positions, names = [], [], None
    with (directory / 'population.vcf').open() as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                names = line.rstrip().split('\t')[9:]
            elif not line.startswith('#') and line.strip():
                fields = line.rstrip().split('\t')
                if ',' in fields[4]:
                    raise ValueError('Unexpected multiallelic record')
                gt = fields[8].split(':').index('GT')
                alleles = [v.split(':')[gt].replace('|', '/').split('/') for v in fields[9:]]
                if len(alleles) != N or any(len(a) != 2 or any(x not in ('0', '1') for x in a) for a in alleles):
                    raise ValueError('Expected complete diploid 0/1 genotypes')
                records.append([sum(map(int, a)) for a in alleles])
                positions.append(int(fields[1])-1)
    if names is None or len(names) != N:
        raise ValueError('VCF sample count differs from census size')
    positions = np.asarray(positions)
    length = params['locus_length']
    if len(set(positions)) != len(positions) or np.any(positions % (length+1) >= length) or np.any(positions >= params['n_loci']*(length+1)-1):
        raise ValueError('Unexpected genomic positions')
    histories = pd.read_csv(directory / 'histories.csv')
    if histories.vcf_column.tolist() != list(range(N)) or (histories.selfing_generations < 0).any():
        raise ValueError('Invalid or unresolved selfing histories')
    histories['vcf_sample'] = names
    return np.asarray(records, dtype=np.int8).reshape(-1, N), histories


@lru_cache(maxsize=3)
def log_probabilities(n):
    probabilities = [dgs_probabilities(s=s, N=1, n_diploids=n) for s in CANDIDATES]
    configs = sorted(probabilities[0])
    with np.errstate(divide='ignore'):
        logs = np.log([[p[g] for g in configs] for p in probabilities])
    return configs, logs


def fit_population(directory, params):
    """Hold the selected individuals fixed across loci; fit each nested subset."""
    saved = json.loads((directory / 'simulation.json').read_text())['signature']
    if saved['parameters'] != params or saved['script_sha256'] != hashlib.sha256(MODEL.read_bytes()).hexdigest():
        raise ValueError('Simulation parameters/model do not match this analysis')
    N, s = params['census_size'], params['selfing_rate']
    dosage, histories = read_population(directory, params)
    selected = np.random.default_rng(seed_for(N, s, 1)).choice(N, 100, replace=False)
    members = histories.iloc[selected].copy()
    members['selection_order'] = range(100)
    curves, spectra = [], []
    for n in SAMPLE_SIZES:
        sample = dosage[:, selected[:n]]
        triples = np.stack([(sample == k).sum(axis=1) for k in range(3)], axis=1)
        counts = Counter(map(tuple, triples.tolist()))
        # Conditional on segregating sites: exclude both monomorphic cells.
        counts.pop((n, 0, 0), None)
        counts.pop((0, 0, n), None)
        configs, logs = log_probabilities(n)
        values = np.array([counts[g] for g in configs], dtype=float)
        if values.sum() == 0 or values.sum() != sum(counts.values()):
            raise ValueError('Empty or unsupported genotype spectrum')
        observed = values > 0  # Avoid 0 * log(0).
        ll = logs[:, observed] @ values[observed]
        if not np.isfinite(ll).all():
            raise ValueError('Nonfinite likelihood')
        meta = dict(census_size=N, true_s=s, n_diploids=n, replicate=0, n_loci=params['n_loci'])
        curves.extend(dict(meta, candidate_s=candidate, loglik=float(value),
                           delta_loglik=float(value-ll.max())) for candidate, value in zip(CANDIDATES, ll))
        spectra.extend(dict(n_diploids=n, n0=g[0], n1=g[1], n2=g[2], count=count)
                       for g, count in sorted(counts.items()))
        members[f'selected_n{n}'] = members.selection_order < n
    members.to_csv(directory / 'sampled_individuals.csv', index=False)
    pd.DataFrame(spectra).to_csv(directory / 'observed_spectra.csv', index=False)
    write_json(directory / 'analysis.json', dict(model_version=model_version, sample_seed=seed_for(N, s, 1),
        candidate_grid=CANDIDATES, sample_sizes=SAMPLE_SIZES,
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    temporary = directory / 'likelihoods.csv.tmp'
    pd.DataFrame(curves).to_csv(temporary, index=False)
    temporary.replace(directory / 'likelihoods.csv')


def population_dir(output, N, s):
    return output / 'populations' / f'N{N}' / f's_{s:.2f}'


def collect(output):
    frames = []
    for N, s in TASKS:
        frame = pd.read_csv(population_dir(output, N, s) / 'likelihoods.csv')
        if len(frame) != 60 or not frame.census_size.eq(N).all() or not np.allclose(frame.true_s, s):
            raise ValueError(f'Incomplete or mismatched results for N={N}, s={s}')
        frames.append(frame)
    temporary = output / 'likelihoods.csv.tmp'
    pd.concat(frames, ignore_index=True).to_csv(temporary, index=False)
    temporary.replace(output / 'likelihoods.csv')
    print(f'Collected all 180 fits into {output / "likelihoods.csv"}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task-index', type=int, choices=range(60), metavar='0..59')
    parser.add_argument('--output', type=Path, default=HERE / 'results')
    parser.add_argument('--slim', default='slim')
    parser.add_argument('--analyze-only', action='store_true')
    parser.add_argument('--collect-only', action='store_true')
    parser.add_argument('--describe', action='store_true')
    args = parser.parse_args()
    if args.describe:
        print(pd.DataFrame([dict(task=i, **parameters(N, s), sample_seed=seed_for(N, s, 1))
                            for i, (N, s) in enumerate(TASKS)]).to_string(index=False))
        return
    if args.collect_only:
        collect(args.output)
        return
    tasks = TASKS if args.task_index is None else [TASKS[args.task_index]]
    for N, s in tasks:
        directory = population_dir(args.output, N, s)
        params = parameters(N, s)
        if not args.analyze_only:
            simulate(params, directory, args.slim)
        fit_population(directory, params)
        print(f'Finished N={N}, s={s:.2f}; fitted n=20, 50, 100', flush=True)
    if args.task_index is None:
        collect(args.output)


if __name__ == '__main__':
    main()
