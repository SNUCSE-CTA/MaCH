## Hypergraph Pattern Matching

**Efficient Hypergraph Pattern Matching via Match-and-Filter and Intersection Constraint**

## Dependencies
- g++ compiler with C++20 support

## Datasets
Datasets used in the experiments can be downloaded from: `https://www.cs.cornell.edu/~arb/data/`

The `example` directory contains the example data and query corresponding to Figure 1 of our manuscript.

## Usage
To build the binaries, use the following commands:
```bash
mkdir build
cd build
cmake .. && make
```
This will generate two executables: `MaCH` and `HypergraphPreprocessing`.

Build an index file from the node-label and hyperedge files:
```bash
./HypergraphPreprocessing -l [path_to_node_labels] -e [path_to_hyperedge_file] -o [path_to_output]
```

Using the generated index file, you can perform hypergraph pattern matching:
```bash
./MaCH -d [path_to_index_file] -q [path_to_query]
```

To print the embeddings, add the `-p` flag. (Note: Using `-p` may significantly slow down execution if there are many embeddings.)

## Example Usage 
```bash
./HypergraphPreprocessing -l ../example/node-labels-example.txt -e ../example/hyperedges-example.txt -o ../example/index-example.out 
./MaCH -d ../example/index-example.out  -q ../example/query-example.txt -p
```
The example above demonstrates the execution using the sample data in the example directory.

