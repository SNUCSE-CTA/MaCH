#pragma once
#include "Base/Base.h"
#include "Base/FileIO.h"
#include "Base/Timer.h"
#include <map>
#include <vector>

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

  // A[i][j] = k: Vertex i is k-th vertex of it's j-th incident hyperedge.
  std::vector<std::vector<int>> inverse_vertex_index;
  // A[i][j] = k: Hyperedge i is k-th incidence hyperedge of it's j-th vertex.
  std::vector<std::vector<int>> inverse_hyperedge_index;

public:
  HyperGraph() {};
  ~HyperGraph() {};

  void ReadHyperGraph(const std::string &filename);
  void LoadHyperGraphDataset(std::string dataset, std::string path);
  void BuildIncidenceList();
  void BuildNeighborIndex();

  inline int GetVertexLabel(int v) { return vertex_label[v]; }
  inline int GetHyperedgeLabel(int e) { return hyperedge_label[e]; }
  inline int GetDegree(int v) const { return incidence_list[v].size(); }
  inline int GetNumVertices() { return num_vertex; }
  inline int GetNumHyperedges() { return num_edge; }
  inline int GetTotalArity() { return total_arity; }
  inline int GetNumVertexLabels() { return num_vertex_labels; }
  inline int GetNumHyperedgeLabels() { return num_hyperedge_labels; }
  inline std::vector<std::vector<int>> &GetHyperedges() { return hyperedges; }
  inline std::vector<int> &GetHyperedge(int idx) { return hyperedges[idx]; }
  inline std::vector<std::vector<int>> &GetIncidenceList() {
    return incidence_list;
  }
  inline std::vector<int> &GetIncidentHyperedges(const int idx) {
    return incidence_list[idx];
  }
  inline int GetIncidentHyperedge(const int u, const int i) {
    return incidence_list[u][i];
  }
  inline int GetContainedVertex(const int e, const int i) {
    return hyperedges[e][i];
  }
  inline std::vector<std::vector<int>> &GetHyperedgeSignatures() {
    return hyperedge_signatures;
  }
  inline std::vector<int> &GetHyperedgeSignature(const int idx) {
    return hyperedge_signatures[idx];
  }
  inline int GetDegree(const int idx) {
    return (int)incidence_list[idx].size();
  }
  inline int GetArity(const int idx) { return (int)hyperedges[idx].size(); }

  inline int GetInverseVertexIndex(int i, int j) {
    return inverse_vertex_index[i][j];
  }
  inline int GetInverseHyperedgeIndex(int i, int j) {
    return inverse_hyperedge_index[i][j];
  }

  void PrintStatistics(std::string name);

  void Preprocess(std::string &hyperedge_path, std::string &vertex_label_path,
                  std::string &output_path);
};

// Read Benson's Hypergraph Format
void HyperGraph::LoadHyperGraphDataset(std::string dataset, std::string path) {
  std::string hyperedge_file =
      path + "/" + dataset + "/hyperedges-" + dataset + ".txt";
  std::string vertex_label_file =
      path + "/" + dataset + "/node-labels-" + dataset + ".txt";
  std::cout << hyperedge_file << " " << fileSize(hyperedge_file.c_str())
            << std::endl;
  std::cout << vertex_label_file << std::endl;
  std::ifstream fin(vertex_label_file);
  std::string line;
  while (getline(fin, line)) {
    vertex_label.push_back(stoi(line) - 1);
  }
  num_vertex = vertex_label.size();
  fin = std::ifstream(hyperedge_file);
  while (getline(fin, line)) {
    auto edge = parse(line, ",");
    hyperedges.push_back(std::vector<int>());
    auto &E = hyperedges.back();
    for (auto &elem : edge) {
      E.push_back(stoi(elem) - 1);
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
  BuildIncidenceList();
}

void HyperGraph::BuildIncidenceList() {
  incidence_list.resize(num_vertex);
  for (int i = 0; i < num_edge; i++) {
    for (int &elem : hyperedges[i]) {
      incidence_list[elem].push_back(i);
    }
  }
}

void HyperGraph::BuildNeighborIndex() {
  inverse_vertex_index.clear();
  inverse_hyperedge_index.clear();

  inverse_vertex_index.resize(GetNumVertices());
  for (int i = 0; i < GetNumVertices(); i++)
    inverse_vertex_index[i].resize(GetDegree(i));
  inverse_hyperedge_index.resize(GetNumHyperedges());
  for (int i = 0; i < GetNumHyperedges(); i++)
    inverse_hyperedge_index[i].resize(GetArity(i));
  std::vector<int> num_seen_(std::max(GetNumVertices(), GetNumHyperedges()), 0);
  for (int i = 0; i < GetNumVertices(); i++) {
    for (int j = 0; j < GetDegree(i); j++) {
      inverse_vertex_index[i][j] = num_seen_[incidence_list[i][j]]++;
    }
  }
  std::fill(num_seen_.begin(), num_seen_.end(), 0);
  for (int i = 0; i < GetNumHyperedges(); i++) {
    for (int j = 0; j < GetArity(i); j++) {
      inverse_hyperedge_index[i][j] = num_seen_[hyperedges[i][j]]++;
    }
  }
}

void HyperGraph::PrintStatistics(std::string name) {
  //        cout<<"\e\[0;32mProgram version: " << __head_version <<
  //        "\e[0m"<<endl;
  fprintf(log_to, "\e\[0;31m[Statistics] %s \e[0m\n", name.c_str());
  fprintf(log_to, "\e\[0;31mV, E, TotalArity = %d, %d, %d \e[0m\n", num_vertex,
          num_edge, total_arity);
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
  std::map<std::vector<int>, std::vector<int>> hyperedges_by_signatures;
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
                if (vertex_label[a] == vertex_label[b])
                  return a < b;
                return vertex_label[a] < vertex_label[b];
              });
    current_hyperedge.erase(
        std::unique(current_hyperedge.begin(), current_hyperedge.end()),
        current_hyperedge.end());
    if (current_hyperedge.size() == 1) {
      continue;
    }
    std::vector<int> current_signature(current_hyperedge.size());
    std::transform(current_hyperedge.begin(), current_hyperedge.end(),
                   current_signature.begin(),
                   [this](int &elem) -> int { return vertex_label[elem]; });
    total_arity += current_hyperedge.size();
    hyperedges.emplace_back(current_hyperedge);
    auto &v = hyperedges_by_signatures[current_signature];
    v.insert(v.end(), current_hyperedge.begin(), current_hyperedge.end());
    index_build_timer.Stop();
  }
  for (auto [signature, h] : hyperedges_by_signatures) {
    PrintVector(signature, fout, true);
    PrintVector(h, fout, true);
  }
  fprintf(stderr, "Wrote %.02lf KB to %s\n",
          (fout.tellp() - start_pos) / 1024.0, output_path.c_str());
  fprintf(stderr, "Preprocess time: %.02lf milliseconds\n",
          index_build_timer.GetTime());
}
} // namespace GraphLib
