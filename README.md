# Private-aggregation (PnF)

Reference benchmarks for **Product-and-Factorization (PnF)**, a framework for
*exact* and *verifiable* private aggregation on a **single untrusted cloud**.

## Brief summary

PnF eliminates the **non-collusion assumption** that multi-server private-aggregation
designs require. Instead of a secret-permutation shuffle, contributors encrypt
their values, the untrusted cloud homomorphically multiplies the ciphertexts,
and the product is decrypted and factorized over a field to recover the inputs
as an unordered multiset. Built on RLWE-based fully homomorphic encryption
with threshold decryption, PnF provides **unlinkability** and **sender
verifiability** (via a lightweight tag) from the base construction; standard
lattice zero-knowledge proofs add security against fully malicious parties, and the
construction is post-quantum.

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

## Basic requirements

- C++17, CMake (higher than 3.16)
- [OpenFHE](https://github.com/openfheorg/openfhe-development) (BFV) — `OpenFHEConfig.cmake`
  expected under `/usr/local/lib/OpenFHE`
- [NTL](https://libntl.org/) + [GMP](https://gmplib.org/) — interpolation and root-finding over `F_t`
- OpenMP (`libomp`) — `CMakeLists.txt` points at the Homebrew path
  `/usr/local/opt/libomp/lib`; adjust if your install differs

## Build options

```sh
cd bench
cmake -S . -B build
cmake --build build -j
```

### Notes for isolation

- Pin the thread count via the environment at process start and prevent the
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
