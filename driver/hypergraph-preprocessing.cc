#include <iostream>
#include <algorithm>
#include <iomanip>
#include <any>

#include "DataStructure/HyperGraph.h"
#include "SubhypergraphMatching/DataHyperGraph.h"
#include "SubhypergraphMatching/PatternHyperGraph.h"
using namespace std;
using namespace GraphLib;

using SubHyperGraphMatching::DataHyperGraph;
using SubHyperGraphMatching::PatternHyperGraph;

int32_t main(int argc, char *argv[]) {
    std::string dataset = "mathoverflow-answers";
    std::string dataset_path = "../../data/";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'p':
                    dataset_path = argv[i + 1];
                    break;
                case 'd':
                    dataset = argv[i + 1];
                    break;
            }
        }
    }
    std::string hyperedge_path = dataset_path + dataset + "/hyperedges-" + dataset + ".txt";
    std::string labels_path = dataset_path + dataset + "/node-labels-" + dataset + ".txt";
    std::string output_path = "./"+ dataset + ".out";
    GraphLib::HyperGraph HG;
    HG.Preprocess(hyperedge_path, labels_path, output_path);
}
