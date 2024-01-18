# Smith-Waterman Algorithm Implementation

## Introduction

The Smith-Waterman (SW) algorithm is a fundamental local sequence alignment algorithm in bioinformatics. This project presents two versions of the SW developed in C++: a basic implementation and an optimized version using AVX2 instructions. 

For testing purposes, we use influenza A virus (H1N1) strain data, specifically the initial sequenced strains from the US CDC dated April 27, 2009:
- A/Texas/05/2009(H1N1)
- A/California/04/2009(H1N1)

## Getting Started

### Prerequisites
Ensure Python is installed for data fetching and a C++ compiler is available for building and running the algorithm.

### Fetching Data
Run the Python script `fetch_data.py` to obtain sequence data for the influenza strains:

```bash
python fetch_data.py
```

Retrieve the specified viral sequences and save them in the `data` directory with this script.

### Building and Running the Algorithm
Utilize the Makefile for straightforward compilation and execution. Build and run the main application with:

```bash
make run
```

Compile the C++ code and execute the Smith-Waterman algorithm using the fetched data through this command.

### Testing
Execute tests to validate the implementation's correctness against known sequences:

```bash
python test.py
```

Compare the Smith-Waterman algorithm's output in both basic and AVX2-optimized versions against expected results using this method.