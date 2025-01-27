#include <iostream>
#include "DataStructure/HyperGraph.h"
#include "SubhypergraphMatching/DataHyperGraph.h"
#include "SubhypergraphMatching/PatternHyperGraph.h"
#include "Base/Timer.h"
#include "SubhypergraphMatching/MaCH.h"
using namespace std;
using namespace GraphLib;

Timer timer;

using SubHyperGraphMatching::DataHyperGraph;
using SubHyperGraphMatching::PatternHyperGraph;

int32_t main(int argc, char *argv[]) {
    std::string data_file, query_file;
    SubHyperGraphMatching::SubHyperGraphMatchingOption opt;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'd':
                    data_file = argv[i + 1];
                    break;
                case 'q':
                    query_file = argv[i + 1];
                    break;
                case 'p':
                    opt.print_answer = true;
                    break;
            }
        }
    }
    if (!file_ok(data_file)) {
        fprintf(stderr, "Data file %s does not exist\n", data_file.c_str());
        return 1;
    }
    if (!file_ok(query_file)) {
        fprintf(stderr, "Query file %s does not exist\n", query_file.c_str());
        return 1;
    }
    PatternHyperGraph PG;
    PG.ReadPatternHyperGraph(query_file);
    PG.PrintStatistics("PatternGraph");
    Timer pre_timer, pre_timer_2; pre_timer.Start();
    DataHyperGraph HG;
    HG.LoadDataGraph(data_file, PG);
    pre_timer.Stop();
    HG.PrintStatistics("Extracted DataGraph");
    SubHyperGraphMatching::MaCH MaCH(&HG, &PG, opt);
}
