#pragma once
#include <set>
#include <unordered_set>

#include "SubhypergraphMatching/DataHyperGraph.h"
#include "SubhypergraphMatching/PatternHyperGraph.h"
#define LARGEGRAPH

namespace GraphLib::SubHyperGraphMatching {

    class CellSignature {
        DataHyperGraph *data = nullptr;
        PatternHyperGraph *query = nullptr;
        unsigned long long *data_profiles, *query_profiles;
        unsigned long long *qbit, *dbit;

       public:
        void Initialize() {
            qbit = (unsigned long long *)calloc(query->GetNumVertices(),
                                                sizeof(unsigned long long));
            dbit = (unsigned long long *)calloc(data->GetNumVertices(),
                                                sizeof(unsigned long long));
            data_profiles = (unsigned long long *)calloc(
                data->GetNumVertices(), sizeof(unsigned long long));
            query_profiles = (unsigned long long *)calloc(
                query->GetNumVertices(), sizeof(unsigned long long));

            unsigned long long offset = (1ull << query->GetNumHyperedges());

            for (int i = 0; i < data->GetNumVertices(); i++) {
                data_profiles[i] = (offset * data->GetVertexLabel(i));
            }
            for (int i = 0; i < query->GetNumVertices(); i++) {
                query_profiles[i] = (offset * query->GetVertexLabel(i));
            }
        }

        CellSignature(DataHyperGraph *data_, PatternHyperGraph *query_) {
            data = data_;
            query = query_;
            Initialize();
        };

        void AddQueryMapping(int qe) {
            int i = 0;
            for (int v : query->GetHyperedge(qe)) {
                qbit[i++] = query_profiles[v];
            }
            std::sort(qbit, qbit + query->GetArity(qe));
            return;
        }

        bool CheckMapping(int qe, int he) {
            int i = 0;
            for (int v : data->GetHyperedge(he)) {
                dbit[i++] = data_profiles[v];
            }
            std::sort(dbit, dbit + query->GetArity(qe));
            for (int i = 0; i < query->GetArity(qe); i++) {
                if (qbit[i] != dbit[i]) return false;
            }
            return true;
        }

        void RemoveQueryMapping(int qe) { return; }

        void AddMapping(int qe, int he) {
            for (int v : data->GetHyperedge(he)) {
                data_profiles[v] = (data_profiles[v] | (1ull << qe));
            }
            for (int v : query->GetHyperedge(qe)) {
                query_profiles[v] = (query_profiles[v] | (1ull << qe));
            }
        }

        void RemoveMapping(int qe, int he) {
            for (int v : data->GetHyperedge(he)) {
                data_profiles[v] = (data_profiles[v] - (1ull << qe));
            }
            for (int v : query->GetHyperedge(qe)) {
                query_profiles[v] = (query_profiles[v] - (1ull << qe));
            }
        }
    };

}  // namespace GraphLib::SubHyperGraphMatching
