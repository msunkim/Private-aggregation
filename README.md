# Private-aggregation (PnF)

Reference benchmarks for **Product-and-Factorization (PnF)**, a framework for
*exact* and *verifiable* private aggregation on a **single untrusted cloud**.

## Summary

PnF eliminates the non-collusion assumption that multi-server private-aggregation
designs require. Instead of a secret-permutation shuffle, contributors encrypt
their values, the untrusted cloud homomorphically **multiplies** the ciphertexts,
and the product is **decrypted and factorized over a field** to recover the inputs
as an unordered multiset. Privacy rests on **threshold decryption** by a small
committee (tolerating up to a specified number of corrupt members), not on
non-collusion among aggregators. Built on RLWE-based fully homomorphic encryption
with threshold decryption, PnF provides **unlinkability** and **sender
verifiability** (via a lightweight tag) from the base construction; standard
lattice zero-knowledge proofs add security against fully malicious parties, and the
construction is **post-quantum**. Each contributor and committee member performs
work independent of the population size; the only super-linear cost—the homomorphic
product—falls on the resource-rich cloud.

## Repository layout

```
bench/
  CMakeLists.txt          build for all targets below
  profile_b_bench.cpp     per-party round cost vs n  (enc / product / partial-decrypt / recovery)
  ksweep.cpp              committee-dependent cost vs committee size k
  sigma_bench.cpp         zero-knowledge proof cost (POWF / PORK / POCD: generate, verify, size)
  thread_scale.cpp        aggregator product time vs OpenMP thread count
  ext_spike.cpp           F_{t^2} longer-message extension: single-round feasibility
  ext_bench.cpp           F_{t^2}: aggregator product + recovery vs n
  ext_enc.cpp             F_{t^2}: per-sender encryption (two ciphertexts) vs n
  ctsize.cpp, diag.cpp, probe_t.cpp,
  profile_b128.cpp, sigma_spike.cpp, threshold_spike.cpp   supporting probes / diagnostics
```

## Requirements

- C++17, CMake ≥ 3.16
- [OpenFHE](https://github.com/openfheorg/openfhe-development) (BFV) — `OpenFHEConfig.cmake`
  expected under `/usr/local/lib/OpenFHE`
- [NTL](https://libntl.org/) + [GMP](https://gmplib.org/) — interpolation and root-finding over `F_t`
- OpenMP (`libomp`) — `CMakeLists.txt` points at the Homebrew path
  `/usr/local/opt/libomp/lib`; adjust if your install differs

## Build

```sh
cd bench
cmake -S . -B build
cmake --build build -j
```

## Run / reproduce

Each executable prints a table to stdout (some also archive a `results_*.csv`).
Run from `bench/build`:

| Executable | Args (defaults) | Reproduces |
|---|---|---|
| `profile_b_bench` | `[trials=10]` | per-party round cost vs `n ∈ {50,100,200,300,500,1000}` |
| `ksweep`          | `[trials=5]`  | committee cost vs `k ∈ {2,5,10,20}` at `n=1000` |
| `sigma_bench`     | `[N=32768] [trials=5]` | ZK proof generate/verify/size (POWF, PORK, POCD) |
| `thread_scale`    | `[n=300] [trials=3] [only_threads=0]` | product time vs threads `{1,2,4,8,12}` (`only_threads>0` measures just that count) |
| `ext_bench`       | `[trials=3]`  | `F_{t^2}` aggregator product + recovery vs `n ∈ {50..500}` |
| `ext_enc`         | `[trials=10]` | `F_{t^2}` per-sender encryption vs `n` |
| `ext_spike`       | `[n=8] [dpad=0]` | single-round `F_{t^2}` feasibility (correctness + per-phase timing) |

Examples:

```sh
cd bench/build
./profile_b_bench 10        # per-party round costs
./ksweep 5                  # committee-size sweep
./sigma_bench 32768 5       # ZK proof costs
./ext_bench 3               # F_{t^2} aggregator cost
./ext_enc 10                # F_{t^2} sender cost
./thread_scale 300 5 12     # product time at 12 threads, n=300
```

### Notes for stable timing

- Pin the thread count via the **environment at process start** and prevent the
  machine from sleeping:

  ```sh
  caffeinate -i env OMP_NUM_THREADS=12 ./profile_b_bench 10
  ```

  On many-core machines this matters: setting `OMP_NUM_THREADS` at startup binds
  the threads cleanly, whereas changing the count at runtime can place them poorly
  and inflate timings.

- For the parallelism scaling, measure each thread count in a **fresh process**,
  e.g. the single-thread baseline as
  `env OMP_NUM_THREADS=1 ./thread_scale 300 5 1`.

- Reported figures are means over the trial count (sample SD where shown);
  expect a few percent run-to-run variation depending on machine load.
