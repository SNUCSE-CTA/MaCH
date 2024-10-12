#pragma once
#include <algorithm>

#include "DataStructure/HyperGraph.h"

namespace GraphLib {
    namespace SubHyperGraphMatching {
        class PatternHyperGraph : public HyperGraph {
            std::map<std::vector<int>, int> hyperedge_label_map;
            std::unordered_map<int, int> vertex_label_map;
            std::map<std::vector<int>, int> intersection_label_map;

            std::vector<std::vector<int>> NEC_partition;
            std::vector<int> NEC_idx, NEC_representative;
            std::vector<std::vector<int>> ContainedNEC;
            std::vector<std::vector<int>> NEC_incidence_list;
            std::vector<std::vector<int>> NEC_inverse_hyperedge_index;
            std::vector<std::vector<int>> inverse_NEC_index;
            std::vector<int> NEC_label;
            int num_NEC;

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
            int GetNECIndex(const int v) { return NEC_idx[v]; }
            int GetNECRepresentative(const int v) {
                return NEC_representative[v];
            }

            // Returns j, which the NEC-head of e[i] in hyperedge e is e[j]
            inline int GetNECReprIndexInEdge(const int e, const int i) {
                int v = hyperedges[e][i];
                int e_for_ei = inverse_hyperedge_index[e][i];
                int nec_head = NEC_representative[v];
                return inverse_vertex_index[nec_head][e_for_ei];
            }

            std::vector<int> &GetNECPartition(const int idx) {
                return NEC_partition[idx];
            }
            void BuildNEC();

            inline int GetNumNECs() { return num_NEC; }
            inline int GetContainedNEC(const int e, const int i) {
                return ContainedNEC[e][i];
            }
            inline int GetNECArity(const int e) {
                return ContainedNEC[e].size();
            }
            inline int GetNECDegree(const int idx) {
                return (int)NEC_incidence_list[idx].size();
            }
            inline std::vector<int> &GetNECIncidentHyperedges(const int idx) {
                return NEC_incidence_list[idx];
            }
            inline int GetNECIncidentHyperedge(const int n, const int i) {
                return NEC_incidence_list[n][i];
            }
            inline int GetNECLabel(const int idx) {
                return (int)NEC_label[idx];
            }
            inline int GetNECSize(const int idx) {
                return (int)NEC_partition[idx].size();
            }
            inline int GetNECInverseHyperedgeIndex(int i, int j) {
                return NEC_inverse_hyperedge_index[i][j];
            }
            inline int GetInverseNECIndex(int i, int j) {
                return inverse_NEC_index[i][j];
            }
            std::map<std::vector<int>, int> &GetIntersectionLabelMap() {
                return intersection_label_map;
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
            BuildNeighborIndex();
            BuildNEC();
        }

        void PatternHyperGraph::BuildNEC() {
            NEC_idx.resize(GetNumVertices());
            NEC_representative.resize(GetNumVertices());
            std::vector<int> vertex_order(GetNumVertices(), 0);
            std::iota(std::begin(vertex_order), std::end(vertex_order), 0);
            std::sort(std::begin(vertex_order), std::end(vertex_order),
                      [&](int i, int j) {
                          if (vertex_label[i] != vertex_label[j]) {
                              return vertex_label[i] < vertex_label[j];
                          }
                          return incidence_list[i] < incidence_list[j];
                      });
            NEC_partition.push_back({vertex_order[0]});
            for (int i = 1; i < GetNumVertices(); i++) {
                if (vertex_label[vertex_order[i]] !=
                        vertex_label[vertex_order[i - 1]] ||
                    incidence_list[vertex_order[i]] !=
                        incidence_list[vertex_order[i - 1]]) {
                    NEC_partition.push_back({vertex_order[i]});
                } else {
                    NEC_partition.back().emplace_back(vertex_order[i]);
                }
            }
            num_NEC = NEC_partition.size();
            for (int i = 0; i < NEC_partition.size(); i++) {
                for (auto x : NEC_partition[i]) {
                    NEC_idx[x] = i;
                    NEC_representative[x] = NEC_partition[i].front();
                }
            }
            NEC_incidence_list.resize(num_NEC);
            NEC_label.resize(num_NEC);
            for (int i = 0; i < NEC_partition.size(); i++) {
                NEC_incidence_list[i] =
                    incidence_list[NEC_partition[i].front()];
                NEC_label[i] = vertex_label[NEC_partition[i].front()];
            }
            ContainedNEC.resize(GetNumHyperedges());
            for (int e = 0; e < GetNumHyperedges(); e++) {
                for (int u_idx = 0; u_idx < GetArity(e); u_idx++) {
                    ContainedNEC[e].push_back(
                        NEC_idx[GetContainedVertex(e, u_idx)]);
                }
                std::sort(ContainedNEC[e].begin(), ContainedNEC[e].end());
                ContainedNEC[e].erase(
                    std::unique(ContainedNEC[e].begin(), ContainedNEC[e].end()),
                    ContainedNEC[e].end());
            }
            NEC_inverse_hyperedge_index.resize(GetNumHyperedges());
            for (int e = 0; e < GetNumHyperedges(); e++) {
                NEC_inverse_hyperedge_index[e].resize(GetNECArity(e));
                for (int n_idx = 0; n_idx < GetNECArity(e); n_idx++) {
                    int n = GetContainedNEC(e, n_idx);
                    for (int e_idx_from_n = 0; e_idx_from_n < GetNECDegree(n);
                         e_idx_from_n++) {
                        if (NEC_incidence_list[n][e_idx_from_n] == e) {
                            NEC_inverse_hyperedge_index[e][n_idx] =
                                e_idx_from_n;
                            break;
                        }
                    }
                }
            }

            inverse_NEC_index.resize(GetNumNECs());
            for (int n = 0; n < GetNumNECs(); n++) {
                inverse_NEC_index[n].resize(GetNECDegree(n));
                for (int e_idx = 0; e_idx < GetNECDegree(n); e_idx++) {
                    int e = GetNECIncidentHyperedge(n, e_idx);
                    for (int n_idx_from_e = 0; n_idx_from_e < GetNECArity(e);
                         n_idx_from_e++) {
                        if (ContainedNEC[e][n_idx_from_e] == n) {
                            inverse_NEC_index[n][e_idx] = n_idx_from_e;
                        }
                    }
                }
            }
        }

    }  // namespace SubHyperGraphMatching
}  // namespace GraphLib
