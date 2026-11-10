"""Run paired 2D midpoint ensembles or time a 512-square pilot."""

import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
from time import perf_counter

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'jam'))
from jam_gas_film_2d import GasFilm2D


def run_pair(exe, model, steps, noisy, seed, scratch, stride):
    n = model.n
    h = np.ones((n, n))
    saved_steps = sorted(set(range(0, steps+1, stride)) | {1, steps})
    save = set(saved_steps)
    frames = [h]
    jam_seconds = 0.
    preparation_seconds = 0.
    rng = np.random.default_rng(seed)
    with tempfile.TemporaryFile(dir=scratch) as file:
        h.astype('<f8', copy=False).tofile(file)
        for j in range(steps):
            start = perf_counter()
            white = rng.standard_normal(h.shape) if j < noisy else None
            if white is not None:
                white.astype('<f8', copy=False).tofile(file)
            preparation_seconds += perf_counter()-start
            start = perf_counter()
            h = model.advance(h, white)
            jam_seconds += perf_counter()-start
            if j+1 in save:
                frames.append(h)
        file.seek(0)
        start = perf_counter()
        args = (n, model.length, model.attraction, model.theta, model.dt, steps, stride, noisy)
        process = subprocess.run([str(exe), *map(str, args)], stdin=file,
                                 capture_output=True, check=True)
        basilisk_seconds = perf_counter()-start
    jam = np.asarray(frames)
    basilisk = np.frombuffer(process.stdout, dtype='<f8').reshape(jam.shape)
    error = float(np.max(abs(jam-basilisk)))
    assert error < 1e-10
    return jam, basilisk, dict(jam_seconds=jam_seconds, basilisk_seconds=basilisk_seconds,
                              preparation_seconds=preparation_seconds, max_paired_error=error)


def continuous_rate(model):
    q = model.q
    surface = np.zeros_like(q)
    np.divide(1., 4*q, out=surface, where=q > 0)
    mobility = q**2*(1/.24+surface)
    return mobility*(3*model.attraction-.5*q*q)


def shell_rows(model, power, times):
    """The main theory uses the continuum spatial and temporal symbol."""
    q = model.q
    rate = continuous_rate(model)
    decay = -np.expm1(2*times[:, None, None]*rate)
    equilibrium = np.divide(2*model.theta, model.length**2*q**2,
                            out=np.zeros_like(q), where=q > 0)
    theory_on = equilibrium*decay
    theory_off = power[1, 0, 0]*np.exp(2*(times[:, None, None]-model.dt)*rate)
    shell = np.rint(q*model.length/(2*np.pi)).astype(int)
    max_shell = int(np.floor((np.pi/6)*model.length/(2*np.pi)-.5))
    counts = np.bincount(shell.ravel())
    rows = []
    for p, theory in enumerate((theory_on, theory_off)):
        for i, t in enumerate(times):
            spectra = [np.bincount(shell.ravel(), weights=field.ravel())
                       for field in (theory[i], *power[p, :, i])]
            for m in range(1, max_shell+1):
                rows.append((p, t, m, counts[m], 2*np.pi*m/model.length,
                             *(np.sqrt(s[m]/counts[m]) for s in spectra)))
    return rows


def save_spectra(path, rows):
    np.savetxt(path, rows, delimiter=',', comments='',
               header='protocol,time,shell,mode_count,q,rms_theory,rms_jam,rms_basilisk')


def benchmark(exe, out, model, replicas):
    records = {}
    steps = 32
    for label, noisy in (('continuous', steps), ('first_step', 1)):
        measurements = [run_pair(exe, model, steps, noisy, 61120+i, out, steps)[2]
                        for i in range(3)]
        medians = {key: statistics.median(m[key] for m in measurements)
                   for key in measurements[0]}
        print(label, medians, flush=True)
        records[label] = dict(pilot=measurements, median=medians,
                              projected_seconds=replicas*160/steps*
                              sum(medians[k] for k in ('jam_seconds', 'basilisk_seconds',
                                                       'preparation_seconds')))
    result = dict(grid=model.n, length=model.length, dt=model.dt, pilot_steps=steps,
                  target_steps=160, pairs_per_protocol=replicas, protocol=records)
    (out/'benchmark.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


def experiment(exe, out, model, replicas):
    steps = 160
    # Index 0 is the common first-step field, followed by the three plotted times.
    power = np.zeros((2, 2, 4, model.n, model.n))
    for i in range(replicas):
        for protocol, noisy in enumerate((steps, 1)):
            j, b, _ = run_pair(exe, model, steps, noisy, 61100+i, out, 10)
            fields = np.stack((j[[1, 2, 5, 17]], b[[1, 2, 5, 17]]))
            power[protocol] += abs(np.fft.fft2(fields-1.)/model.n**2)**2/replicas
        print(f'{i+1}/{replicas} paired seeds in both protocols', flush=True)
    times = model.dt*np.array([1, 10, 40, 160])
    save_spectra(out/'tcw.csv', [row for row in shell_rows(model, power, times) if row[1] > model.dt])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--n', type=int, default=512)
    parser.add_argument('--replicas', type=int, default=30)
    parser.add_argument('--benchmark', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.n < 16 or args.n > 512 or args.n & (args.n-1) or args.replicas < 1:
        parser.error('n must be a power of two from 16 to 512, replicas must be positive')
    model = GasFilm2D(args.n, 2.*args.n, 8., theta=1e-5)
    action = benchmark if args.benchmark else experiment
    action(args.executable.resolve(), args.output, model, args.replicas)


if __name__ == '__main__':
    main()
