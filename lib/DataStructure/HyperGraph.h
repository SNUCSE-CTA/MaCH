#pragma once
#include <set>
#include <map>
#include <vector>

#include "Base/Base.h"
#include "Base/FileIO.h"
#include "Base/Timer.h"

using std::cout, std::endl;
namespace GraphLib {
    class HyperGraph {
       protected:
        int num_vertex = 0, num_edge = 0, total_arity = 0;
        std::vector<std::vector<int>> hyperedges;
        std::vector<std::vector<int>> incidence_list;

        int num_vertex_labels = 0, num_hyperedge_labels = 0;
        std::vector<int> vertex_label, hyperedge_label;
        std::vector<std::vector<int>> hyperedge_signatures;
        std::vector<int> hyperedge_id_map;

       public:
        HyperGraph() {};
        ~HyperGraph() {};

        void ReadHyperGraph(const std::string &filename);
        void LoadHyperGraphDataset(std::string dataset, std::string path);
        void BuildIncidenceList();

        inline int GetVertexLabel(int v) { return vertex_label[v]; }
        inline int GetHyperedgeLabel(int e) { return hyperedge_label[e]; }
        inline int GetDegree(int v) const { return incidence_list[v].size(); }
        inline int GetNumVertices() { return num_vertex; }
        inline int GetNumHyperedges() { return num_edge; }
        inline int GetNumVertexLabels() { return num_vertex_labels; }
        inline int GetNumHyperedgeLabels() { return num_hyperedge_labels; }
        inline std::vector<int> &GetHyperedge(int idx) {
            return hyperedges[idx];
        }
        inline std::vector<int> &GetIncidentHyperedges(const int idx) {
            return incidence_list[idx];
        }
        inline int GetDegree(const int idx) {
            return (int)incidence_list[idx].size();
        }
        inline int GetArity(const int idx) {
            return (int)hyperedges[idx].size();
        }

        void PrintStatistics(std::string name);

        void Preprocess(std::string &hyperedge_path,
                        std::string &vertex_label_path,
                        std::string &output_path);
    };

    void HyperGraph::BuildIncidenceList() {
        incidence_list.resize(num_vertex);
        for (int i = 0; i < num_edge; i++) {
            for (int &elem : hyperedges[i]) {
                incidence_list[elem].push_back(i);
            }
        }
    }

    void HyperGraph::PrintStatistics(std::string name) {
        //        cout<<"\e\[0;32mProgram version: " << __head_version <<
        //        "\e[0m"<<endl;
        fprintf(log_to, "\e\[0;31m[Statistics] %s \e[0m\n", name.c_str());
        fprintf(log_to, "\e\[0;31mV, E, TotalArity = %d, %d, %d \e[0m\n",
                num_vertex, num_edge, total_arity);
        fprintf(log_to, "\e\[0;31m#VLabel, #ELabel = %d, %d \e[0m\n",
                num_vertex_labels, num_hyperedge_labels);
        fprintf(log_to, "\n");
    }

    void HyperGraph::Preprocess(std::string &hyperedge_path,
                                std::string &vertex_label_path,
                                std::string &output_path) {
        if (!(file_ok(hyperedge_path) and file_ok(vertex_label_path))) {
            fprintf(stderr, "Datagraph not found\n");
            exit(1);
        }
        Timer index_build_timer;
        std::set<std::vector<int>> hyperedge_set;
        std::map<std::vector<int>, std::vector<int>> hyperedges_by_signatures;
        std::map<std::vector<int>, std::vector<int>> hyperedge_ids_by_signatures;
        std::ofstream fout(output_path, std::ios::binary);
        std::streampos start_pos = fout.tellp();
        std::ifstream fin(vertex_label_path);
        std::string line;
        while (getline(fin, line)) {
            vertex_label.push_back(stoi(line));
        }
        num_vertex = vertex_label.size();
        PrintVector(this->vertex_label, fout, true);
        fin = std::ifstream(hyperedge_path);
        while (getline(fin, line)) {
            auto edge = parse(line, ",");
            std::vector<int> current_hyperedge;
            for (auto &elem : edge) {
                current_hyperedge.push_back(stoi(elem) - 1);
            }
            index_build_timer.Start();
            std::sort(current_hyperedge.begin(), current_hyperedge.end(),
                      [this](int &a, int &b) -> bool {
                          if (vertex_label[a] == vertex_label[b]) return a < b;
                          return vertex_label[a] < vertex_label[b];
                      });
            current_hyperedge.erase(
                std::unique(current_hyperedge.begin(), current_hyperedge.end()),
                current_hyperedge.end());
            if (current_hyperedge.size() == 1) {
                continue;
            }
            if (hyperedge_set.contains(current_hyperedge)) {
                continue;
            }
            hyperedge_set.insert(current_hyperedge);
            std::vector<int> current_signature(current_hyperedge.size());
            std::transform(current_hyperedge.begin(), current_hyperedge.end(),
                           current_signature.begin(), [this](int &elem) -> int {
                               return vertex_label[elem];
                           });
            total_arity += current_hyperedge.size();
            hyperedges.emplace_back(current_hyperedge);
            auto &v = hyperedges_by_signatures[current_signature];
            v.insert(v.end(), current_hyperedge.begin(),
                     current_hyperedge.end());
            auto &v_id = hyperedge_ids_by_signatures[current_signature];
            v_id.push_back(hyperedges.size() - 1);
            index_build_timer.Stop();
        }
        for (auto [signature, h] : hyperedges_by_signatures) {
            PrintVector(signature, fout, true);
            PrintVector(h, fout, true);
            PrintVector(hyperedge_ids_by_signatures[signature], fout, true, false);
        }
        fprintf(stderr, "Wrote %.02lf KB to %s\n",
                (fout.tellp() - start_pos) / 1024.0, output_path.c_str());
        fprintf(stderr, "Preprocess time: %.02lf milliseconds\n",
                index_build_timer.GetTime());
    }
}  // namespace GraphLib
