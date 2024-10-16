#include "DataStructure/HyperGraph.h"
#include "SubhypergraphMatching/DataHyperGraph.h"
#include "SubhypergraphMatching/PatternHyperGraph.h"
using namespace std;
using namespace GraphLib;

using SubHyperGraphMatching::DataHyperGraph;
using SubHyperGraphMatching::PatternHyperGraph;

int32_t main(int argc, char *argv[]) {
    std::string hyperedge_path, labels_path, output_path;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'o':
                    output_path = argv[i + 1];
                    break;
                case 'e':
                    hyperedge_path = argv[i + 1];
                    break;
                case 'l':
                    labels_path = argv[i + 1];
                    break;
            }
        }
    }
    HyperGraph HG;
    HG.Preprocess(hyperedge_path, labels_path, output_path);
}
