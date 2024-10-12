#pragma once
#include <algorithm>

#include "DataStructure/HyperGraph.h"

namespace GraphLib {
    namespace SubHyperGraphMatching {
        class PatternHyperGraph : public HyperGraph {
            std::map<std::vector<int>, int> hyperedge_label_map;
            std::unordered_map<int, int> vertex_label_map;
            std::map<std::vector<int>, int> intersection_label_map;


           public:
            PatternHyperGraph() {};
            ~PatternHyperGraph() {};
            PatternHyperGraph &operator=(const PatternHyperGraph &) = delete;
            PatternHyperGraph(const PatternHyperGraph &) = delete;
            void ReadPatternHyperGraph(const std::string &filename);

            int GetMappedVertexLabel(const int l) {
                if (vertex_label_map.find(l) == vertex_label_map.end())
                    return -1;
                else
                    return vertex_label_map[l];
            }
            int GetMappedHyperedgeLabel(const std::vector<int> &l) {
                if (hyperedge_label_map.find(l) == hyperedge_label_map.end())
                    return -1;
                else
                    return hyperedge_label_map[l];
            }
        };

        void PatternHyperGraph::ReadPatternHyperGraph(const string &filename) {
            std::cerr << "Read " << fileSize(filename.c_str()) << " bytes from "
                      << filename << endl;
            std::ifstream fin(filename);
            fin >> num_vertex >> num_edge;
            vertex_label.resize(num_vertex);
            for (int i = 0; i < num_vertex; i++) {
                fin >> vertex_label[i];
            }
            for (int i = 0; i < num_vertex; i++) {
                int l = vertex_label[i];
                if (vertex_label_map.find(l) == vertex_label_map.end()) {
                    vertex_label_map[l] = num_vertex_labels++;
                }
                vertex_label[i] = vertex_label_map[l];
            }

            std::string line;
            for (int i = 0; i < num_edge; i++) {
                fin >> line;
                auto edge = parse(line, ",");
                hyperedges.push_back(std::vector<int>());
                auto &E = hyperedges.back();
                for (auto &elem : edge) {
                    E.push_back(stoi(elem));
                }
                std::sort(E.begin(), E.end());
                E.erase(std::unique(E.begin(), E.end()), E.end());
                if (E.size() == 1) {
                    hyperedges.pop_back();
                    continue;
                }
                total_arity += E.size();
            }
            std::sort(hyperedges.begin(), hyperedges.end());
            hyperedges.erase(std::unique(hyperedges.begin(), hyperedges.end()),
                             hyperedges.end());
            num_edge = hyperedges.size();

            hyperedge_label.resize(GetNumHyperedges());
            for (auto &E : hyperedges) {
                hyperedge_signatures.push_back(std::vector<int>());
                for (auto &elem : E) {
                    hyperedge_signatures.back().push_back(vertex_label[elem]);
                }
                std::sort(hyperedge_signatures.back().begin(),
                          hyperedge_signatures.back().end());
            }
            for (int i = 0; i < GetNumHyperedges(); i++) {
                if (hyperedge_label_map.find(hyperedge_signatures[i]) ==
                    hyperedge_label_map.end()) {
                    hyperedge_label_map[hyperedge_signatures[i]] =
                        num_hyperedge_labels++;
                }
                hyperedge_label[i] =
                    hyperedge_label_map[hyperedge_signatures[i]];
            }
            BuildIncidenceList();
        }

    }  // namespace SubHyperGraphMatching
}  // namespace GraphLib
