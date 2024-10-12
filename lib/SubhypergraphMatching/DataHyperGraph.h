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
  std::unordered_map<int, int> vertex_compression_map_inv;

public:
  std::vector<int> &GetHyperedgesByLabel(const int l) {
    return hyperedges_by_label[l];
  }
  std::vector<int> &GetVerticesByLabel(const int l) {
    return vertex_by_labels[l];
  }
  void LoadDataGraph_old(std::string dataset, std::string path,
                         PatternHyperGraph &P);
  void LoadDataGraph(std::string &index_file, PatternHyperGraph &P);
  int GetOrigVertex(int id);
  std::vector<int> GetOrigHyperEdge(int id);
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
                 tmp_vertex_label.begin(),
                 [this, &P](int x) { return P.GetMappedVertexLabel(x); });
  fprintf(stderr, "Loaded %lu vertices\n", tmp_vertex_label.size());
  num_vertex = tmp_vertex_label.size();
  vertex_compression_map.resize(num_vertex, -1);
  while (!fin.eof()) {
    std::vector<int> signature;
    bool flag = ReadVector(signature, fin, true, -1);
    if (!flag)
      break;
    std::transform(signature.begin(), signature.end(), signature.begin(),
                   [this, &P](int x) { return P.GetMappedVertexLabel(x); });
    std::sort(signature.begin(), signature.end());
    int num_hyperedges_by_signature;
    fin.read(reinterpret_cast<char *>(&num_hyperedges_by_signature),
             sizeof(num_hyperedges_by_signature));
    num_hyperedges_by_signature /= signature.size();
    //                fprintf(stderr, "We will read %d hyperedges\n",
    //                num_hyperedges_by_signature);
    int l = P.GetMappedHyperedgeLabel(signature);
    if (l == -1) {
      unsigned long long skip_bytes = num_hyperedges_by_signature;
      skip_bytes *= signature.size();
      skip_bytes *= sizeof(int);
      fin.seekg(skip_bytes, std::ios::cur);
      //                    fprintf(stderr, "Because we could not find the
      //                    signature, we skip %lld bytes\n", skip_bytes);
      continue;
    }
    //                fprintf(stderr, "Found %d matching hyperedges\n",
    //                num_hyperedges_by_signature);
    std::vector<int> hyperedges_by_signature;
    ReadVector(hyperedges_by_signature, fin, true,
               num_hyperedges_by_signature * signature.size());
    total_arity += hyperedges_by_signature.size();
    for (auto &elem : hyperedges_by_signature) {
      vertex_compression_map[elem] = 0;
    }
    for (int i = 0; i < hyperedges_by_signature.size(); i += signature.size()) {
      hyperedges.emplace_back(hyperedges_by_signature.begin() + i,
                              hyperedges_by_signature.begin() + i +
                                  signature.size());
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
  std::sort(hyperedges.begin(), hyperedges.end());
  hyperedges.erase(std::unique(hyperedges.begin(), hyperedges.end()),
                   hyperedges.end());
  num_edge = hyperedges.size();
  for (auto &E : hyperedges) {
    hyperedge_signatures.push_back(std::vector<int>());
    for (int x : E) {
      hyperedge_signatures[hyperedge_signatures.size() - 1].push_back(
          vertex_label[x]);
    }
    std::sort(hyperedge_signatures[hyperedge_signatures.size() - 1].begin(),
              hyperedge_signatures[hyperedge_signatures.size() - 1].end());
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
  BuildNeighborIndex();

  if (DEBUG) {
    for (int i = 0; i < vertex_compression_map.size(); i++) {
      if (vertex_compression_map[i] >= 0) {
        vertex_compression_map_inv[vertex_compression_map[i]] = i;
      }
    }
  }
}

void DataHyperGraph::LoadDataGraph_old(std::string dataset, std::string path,
                                       PatternHyperGraph &P) {
  std::string hyperedge_file =
      path + "/" + dataset + "/hyperedges-" + dataset + ".txt";
  std::string vertex_label_file =
      path + "/" + dataset + "/node-labels-" + dataset + ".txt";
  if (!(file_ok(hyperedge_file) and file_ok(vertex_label_file))) {
    fprintf(stderr, "Datagraph not found\n");
    exit(1);
  }
  std::cerr << "Read " << fileSize(hyperedge_file.c_str()) << " bytes from "
            << hyperedge_file << endl;
  Timer readtimer;
  readtimer.Start();
  std::ifstream fin(vertex_label_file);
  std::string line;
  std::vector<int> tmp_vertex_label;
  while (getline(fin, line)) {
    int l = stoi(line);
    tmp_vertex_label.push_back(P.GetMappedVertexLabel(l));
  }
  num_vertex = tmp_vertex_label.size();
  vertex_compression_map.resize(num_vertex, -1);
  fin = std::ifstream(hyperedge_file);
  while (getline(fin, line)) {
    auto edge = parse(line, ",");
    std::vector<int> current_hyperedge;
    std::vector<int> current_signature;
    for (auto &elem : edge) {
      current_hyperedge.push_back(stoi(elem) - 1);
    }
    std::sort(current_hyperedge.begin(), current_hyperedge.end());
    current_hyperedge.erase(
        std::unique(current_hyperedge.begin(), current_hyperedge.end()),
        current_hyperedge.end());
    if (current_hyperedge.size() == 1) {
      continue;
    }
    for (auto &elem : current_hyperedge) {
      current_signature.push_back(tmp_vertex_label[elem]);
    }
    std::sort(current_signature.begin(), current_signature.end());
    int l = P.GetMappedHyperedgeLabel(current_signature);
    if (l == -1)
      continue;
    for (auto &elem : current_hyperedge) {
      vertex_compression_map[elem] = 0;
    }
    total_arity += current_hyperedge.size();
    hyperedges.push_back(current_hyperedge);
  }
  readtimer.Stop();
  std::cerr << std::fixed << "Graph read Time: " << readtimer.GetTime() << '\n';
  Timer pretimer;
  pretimer.Start();
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
  std::sort(hyperedges.begin(), hyperedges.end());
  hyperedges.erase(std::unique(hyperedges.begin(), hyperedges.end()),
                   hyperedges.end());
  num_edge = hyperedges.size();
  for (auto &E : hyperedges) {
    hyperedge_signatures.push_back(std::vector<int>());
    for (int x : E) {
      hyperedge_signatures[hyperedge_signatures.size() - 1].push_back(
          vertex_label[x]);
    }
    std::sort(hyperedge_signatures[hyperedge_signatures.size() - 1].begin(),
              hyperedge_signatures[hyperedge_signatures.size() - 1].end());
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
  BuildNeighborIndex();
  pretimer.Stop();
  std::cerr << std::fixed << "Preprocessing Time: " << pretimer.GetTime()
            << '\n';
}

int DataHyperGraph::GetOrigVertex(int i) {
  if (vertex_compression_map_inv.empty()) {
    fprintf(stderr, "GetOrigVertex called without DEBUG info built\n");
    exit(-1);
  }
  return vertex_compression_map_inv[i];
}

std::vector<int> DataHyperGraph::GetOrigHyperEdge(int id) {
  if (vertex_compression_map_inv.empty()) {
    fprintf(stderr, "GetOrigVertex called without DEBUG info built\n");
    exit(-1);
  }
  std::vector<int> ret;
  for (int x : hyperedges[id]) {
    ret.push_back(vertex_compression_map_inv[x]);
  }
  return ret;
}
} // namespace SubHyperGraphMatching
} // namespace GraphLib
