#pragma once
#include "Base/Timer.h"
#include "DataStructure/HyperGraph.h"
#include "SubhypergraphMatching/PatternHyperGraph.h"

namespace GraphLib {
    namespace SubHyperGraphMatching {
        class DataHyperGraph : public GraphLib::HyperGraph {
           protected:
            std::vector<std::vector<int>> vertex_by_labels;
            std::vector<std::vector<int>> hyperedges_by_label;
            std::vector<int> vertex_compression_map;
            std::vector<int> hyperedge_ids;
            std::unordered_map<int, int> vertex_compression_map_inv;

           public:
            std::vector<int> &GetHyperedgesByLabel(const int l) {
                return hyperedges_by_label[l];
            }
            void LoadDataGraph(std::string &index_file, PatternHyperGraph &P);
            int GetOrigHyperedgeId(const int idx) { return hyperedge_ids[idx]; }
        };

        void DataHyperGraph::LoadDataGraph(std::string &index_file,
                                           PatternHyperGraph &P) {
            if (!(file_ok(index_file))) {
                fprintf(stderr, "Datagraph not found\n");
                exit(1);
            }
            std::ifstream fin(index_file, std::ios::binary);
            std::vector<int> tmp_vertex_label;
            ReadVector(tmp_vertex_label, fin, true, -1);
            std::transform(tmp_vertex_label.begin(), tmp_vertex_label.end(),
                           tmp_vertex_label.begin(), [this, &P](int x) {
                               return P.GetMappedVertexLabel(x);
                           });
            num_vertex = tmp_vertex_label.size();
            vertex_compression_map.resize(num_vertex, -1);
            while (!fin.eof()) {
                std::vector<int> signature;
                bool flag = ReadVector(signature, fin, true, -1);
                if (!flag) break;
                std::transform(
                    signature.begin(), signature.end(), signature.begin(),
                    [this, &P](int x) { return P.GetMappedVertexLabel(x); });
                std::sort(signature.begin(), signature.end());
                int num_hyperedges_by_signature;
                fin.read(reinterpret_cast<char *>(&num_hyperedges_by_signature),
                         sizeof(num_hyperedges_by_signature));
                num_hyperedges_by_signature /= signature.size();
                int l = P.GetMappedHyperedgeLabel(signature);
                if (l == -1) {
                    unsigned long long skip_bytes = num_hyperedges_by_signature;
                    skip_bytes *= (signature.size() + 1);
                    skip_bytes *= sizeof(int);
                    fin.seekg(skip_bytes, std::ios::cur);
                    continue;
                }
                std::vector<int> hyperedges_by_signature;
                std::vector<int> hyperedge_ids_by_signature;
                ReadVector(hyperedges_by_signature, fin, true,
                           num_hyperedges_by_signature * signature.size());
                ReadVector(hyperedge_ids_by_signature, fin, true,
                              num_hyperedges_by_signature);
                total_arity += hyperedges_by_signature.size();
                for (auto &elem : hyperedges_by_signature) {
                    vertex_compression_map[elem] = 0;
                }
                for (int i = 0; i < hyperedges_by_signature.size();
                     i += signature.size()) {
                    hyperedges.emplace_back(
                        hyperedges_by_signature.begin() + i,
                        hyperedges_by_signature.begin() + i + signature.size());
                    hyperedge_ids.push_back(hyperedge_ids_by_signature[i / signature.size()]);
                }
            }
            num_vertex = 0;
            for (int i = 0; i < vertex_compression_map.size(); i++) {
                if (vertex_compression_map[i] >= 0) {
                    vertex_label.push_back(tmp_vertex_label[i]);
                    vertex_compression_map[i] = num_vertex++;
                }
            }
            for (int i = 0; i < hyperedges.size(); i++) {
                for (int j = 0; j < hyperedges[i].size(); j++) {
                    hyperedges[i][j] = vertex_compression_map[hyperedges[i][j]];
                }
            }
            for (int i = 0; i < hyperedges.size(); i++)
                hyperedges[i].push_back(hyperedge_ids[i]);
            std::sort(hyperedges.begin(), hyperedges.end());
            for (int i = 0; i < hyperedges.size(); i++) {
                hyperedge_ids[i] = hyperedges[i].back();
                hyperedges[i].pop_back();
            }
            num_edge = hyperedges.size();
            for (auto &E : hyperedges) {
                hyperedge_signatures.push_back(std::vector<int>());
                for (int x : E) {
                    hyperedge_signatures[hyperedge_signatures.size() - 1]
                        .push_back(vertex_label[x]);
                }
                std::sort(hyperedge_signatures[hyperedge_signatures.size() - 1]
                              .begin(),
                          hyperedge_signatures[hyperedge_signatures.size() - 1]
                              .end());
                hyperedge_label.push_back(
                    P.GetMappedHyperedgeLabel(hyperedge_signatures.back()));
            }

            num_vertex_labels = P.GetNumVertexLabels();
            num_hyperedge_labels = P.GetNumHyperedgeLabels();

            hyperedges_by_label.resize(GetNumHyperedgeLabels());
            for (int i = 0; i < GetNumHyperedges(); i++) {
                hyperedges_by_label[GetHyperedgeLabel(i)].push_back(i);
            }
            vertex_by_labels.resize(GetNumVertexLabels());
            for (int i = 0; i < GetNumVertices(); i++) {
                vertex_by_labels[GetVertexLabel(i)].push_back(i);
            }
            BuildIncidenceList();
        }

    }  // namespace SubHyperGraphMatching
}  // namespace GraphLib
