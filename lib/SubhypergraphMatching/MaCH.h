#pragma once
#include <array>
#include <queue>
#include <stack>

#include "SubhypergraphMatching/CellSignature.h"
#include "SubhypergraphMatching/DataHyperGraph.h"
#include "SubhypergraphMatching/PatternHyperGraph.h"

namespace GraphLib::SubHyperGraphMatching {
    using SubHyperGraphMatching::CellSignature;

    struct SubHyperGraphMatchingOption {
        bool print_answer = false;
    };

    class MaCH {
       private:
        SubHyperGraphMatchingOption opt;
        CellSignature *cell_signature;
        SubHyperGraphMatching::DataHyperGraph *hyper_data;
        SubHyperGraphMatching::PatternHyperGraph *hyper_query;

        std::vector<std::vector<int>> query_edge_adj;
        std::vector<std::vector<int>> query_edge_adj_idx, query_edge_adj_label;

        int MaxStage;

        long long *checksum;

        int **cand;

        int **cand_idx;

        int **initial_cand_idx;

        int *initial_cand_sz_64;

        int **cand_sz;

        unsigned long long ***cand_exist;

        long long VqVg = 0, VqVgDu = 0, VqVgDuVg = 0;
        long long Vg_64 = 0, VqVg_64 = 0, VqVgDuVg_64 = 0;
        long long VqCu = 0, VqCuDu = 0;
        long long VqCu_64 = 0, VqCuDuCun_64 = 0;

        unsigned long long ****cand_nbr_exist;

        int ***cand_nbr_sz;

        int ****cand_nbr;

        int ***local_match;

        std::queue<std::pair<int, int>> queue_removed;
        bool **in_queue_removed;

        int *global_match_of_v;

        double *hyperedge_weight;

        int *unmapped_edges;
        int *unmapped_edges_idx;
        int unmapped_edges_sz;

        long long num_embedding = 0;
        long long traversed_node = 0;
        long long conflict = 0;
        long long conflict2 = 0;
        long long num_check_intersection = 0;
        long long num_check_cleaning = 0;

        double estimated_embedding = 0;
        long double progress = 0;

        long long num_restore_local_match = 0;
        long long num_correct_local_match = 0;

        std::vector<int> matched;
        std::stack<int> matched_order;
        std::vector<int> unmapped_deg;

       public:
        MaCH(GraphLib::SubHyperGraphMatching::DataHyperGraph *hyper_data_,
               GraphLib::SubHyperGraphMatching::PatternHyperGraph *hyper_query_,
               SubHyperGraphMatchingOption option);

        ~MaCH();

        inline bool isCandidate(int u, int v, int stage) const {
            return cand_sz[stage][u] > cand_idx[u][v] && cand_idx[u][v] >= 0;
        }

        void BuildInitialHCS();

        bool FilterHCS(int stage = 0);

        void GetQueryConnectivity();

        void Initialize();

        void InitializeMemory();

        bool BuildHyperCandidateSpace();

        void BuildConnection();

        void PrintCSStatistics(int level = 0, int stage = 0);

        inline int GetLogSize(int n);

        bool RestoreLocalMatch(int u, int v, int e_idx, int stage);

        bool backtracking();

        inline void PushCand(int *cand, int *idx, int &sz, int v);

        inline void PushCand(int *cand, int &sz, int v);

        inline int SwapAndPop(int *cand, int *idx, int &sz, int pos);

        inline int SwapAndPop(int *cand, int &sz, int pos);

        inline void EmptyQueue();

        void mapping(int e, int f) {
            matched_order.push(e);
            cell_signature->AddMapping(e, f);
            global_match_of_v[f] = e;
            matched[e] = f;
            for (int e_ : query_edge_adj[e]) {
                unmapped_deg[e_]--;
            }
            SwapAndPop(unmapped_edges, unmapped_edges_idx, unmapped_edges_sz,
                       unmapped_edges_idx[e]);
        }

        void unmapping_last() {
            int e = matched_order.top();
            matched_order.pop();
            cell_signature->RemoveMapping(e, matched[e]);
            global_match_of_v[matched[e]] = -1;
            matched[e] = -1;
            for (int e_ : query_edge_adj[e]) {
                unmapped_deg[e_]++;
            }
            unmapped_edges_sz++;
        }

        int debugging = 0;
    };

    MaCH::MaCH(
        DataHyperGraph *hyper_data_,
        PatternHyperGraph *hyper_query_,
        SubHyperGraphMatchingOption option) {
        opt = option;
        hyper_data = hyper_data_;
        hyper_query = hyper_query_;
        cell_signature = new CellSignature(hyper_data, hyper_query);

        Timer InitTime;
        InitTime.Start();
        Initialize();
        InitTime.Stop();
        Timer AfterMemory;
        AfterMemory.Start();
        bool SuccessHCS = BuildHyperCandidateSpace();
        AfterMemory.Stop();
        PrintCSStatistics(1, 0);
        fprintf(log_to, "HCSBuildingAndFileringTime: %.02lf\n",
                AfterMemory.GetTime());
        traversed_node = 0;
        num_embedding = 0;
        Timer BtTimer;
        BtTimer.Start();

        if (SuccessHCS) backtracking();

        long long ans = num_embedding;

        BtTimer.Stop();
        fprintf(log_to, "SearchTime: %.02lf\n", BtTimer.GetTime());
        fprintf(log_to, "TotalTime: %.02lf\n",
                AfterMemory.GetTime() + BtTimer.GetTime());
        fprintf(log_to, "#Embeddings: %lld\n", ans);
    }

    MaCH::~MaCH() {}

    void MaCH::PrintCSStatistics(int level, int stage) {
        fprintf(log_to, "\e\[0;31m[CS Statistics after Filtering]\n");
        int num_vertices = 0;
        for (int u = 0; u < hyper_query->GetNumHyperedges(); u++) {
            num_vertices += cand_sz[stage][u];
            if (level >= 1) {
                fprintf(log_to, "Hyperedge %d: Matching %d Hyperedges\n",
                        u, cand_sz[stage][u]);
                if (level >= 2) {
                    for (int i = 0; i < cand_sz[stage][u]; i++) {
                        fprintf(log_to, "%d ", cand[u][i]);
                    }
                    fprintf(log_to, "\n");
                }
            }
        }
        fprintf(log_to, "Number of Hyperedges: %d \e[0m\n", num_vertices);
    }

    void MaCH::GetQueryConnectivity() {
        std::vector<int> visited(hyper_query->GetNumHyperedges());
        query_edge_adj =
            std::vector<std::vector<int>>(hyper_query->GetNumHyperedges());
        query_edge_adj_idx = std::vector<std::vector<int>>(
            hyper_query->GetNumHyperedges(),
            std::vector<int>(hyper_query->GetNumHyperedges()));
        query_edge_adj_label = std::vector<std::vector<int>>(
            hyper_query->GetNumHyperedges(),
            std::vector<int>(hyper_query->GetNumHyperedges()));
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            for (int v : hyper_query->GetHyperedge(e)) {
                for (int e_ : hyper_query->GetIncidentHyperedges(v)) {
                    if (visited[e_] || e == e_) continue;
                    query_edge_adj[e].push_back(e_);
                    visited[e_] = 1;
                    query_edge_adj_label[e][e_] =
                        hyper_query->GetVertexLabel(v);
                }
            }
            for (int i = 0; i < query_edge_adj[e].size(); i++) {
                int e_ = query_edge_adj[e][i];
                query_edge_adj_idx[e][e_] = i;
                visited[e_] = 0;
            }
        }
    }

    void MaCH::Initialize() {
        GetQueryConnectivity();

        MaxStage = 1;
        Vg_64 = (hyper_data->GetNumHyperedges() - 1) / 64 + 1;

        matched = std::vector<int>(hyper_query->GetNumHyperedges(), -1);
        unmapped_deg = std::vector<int>(hyper_query->GetNumHyperedges());
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            unmapped_deg[e] = query_edge_adj[e].size();
        }

        unmapped_edges =
            (int *)calloc(hyper_query->GetNumHyperedges(), sizeof(int));
        unmapped_edges_idx =
            (int *)calloc(hyper_query->GetNumHyperedges(), sizeof(int));

        unmapped_edges_sz = hyper_query->GetNumHyperedges();
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            unmapped_edges[e] = e;
            unmapped_edges_idx[e] = e;
        }

        initial_cand_sz_64 =
            (int *)calloc(hyper_query->GetNumHyperedges(), sizeof(int));

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            int e_label = hyper_query->GetHyperedgeLabel(e);
            auto &Hyperedges_with_same_label =
                hyper_data->GetHyperedgesByLabel(e_label);
            int num_same_label = Hyperedges_with_same_label.size();
            initial_cand_sz_64[e] =
                (hyper_data->GetHyperedgesByLabel(e_label).size() - 1) / 64 + 1;

            MaxStage += 1;
            VqCu += num_same_label;
            VqVg_64 += Vg_64;
            for (int f = 0; f < hyper_data->GetNumHyperedges(); f++) {
                VqVg++;
                VqVgDu += query_edge_adj[e].size();
                VqVgDuVg +=
                    query_edge_adj[e].size() * hyper_data->GetNumHyperedges();
                VqVgDuVg_64 += query_edge_adj[e].size() * Vg_64;
            }
        }

        hyperedge_weight =
            (double *)calloc(hyper_query->GetNumHyperedges(), sizeof(double));
        cand_sz = (int **)calloc(MaxStage, sizeof(int *));
        int *cand_sz_1 = (int *)calloc(
            MaxStage * hyper_query->GetNumHyperedges(), sizeof(int));

        for (int s = 0; s < MaxStage; s++) {
            cand_sz[s] = cand_sz_1;
            cand_sz_1 += hyper_query->GetNumHyperedges();
        }

        cand = (int **)calloc(hyper_query->GetNumHyperedges(), sizeof(int *));
        int *cand_1 = (int *)calloc(VqCu, sizeof(int));

        cand_idx =
            (int **)calloc(hyper_query->GetNumHyperedges(), sizeof(int *));
        int *cand_idx_1 = (int *)calloc(VqVg, sizeof(int));
        memset(cand_idx_1, -1, VqVg * sizeof(int));

        initial_cand_idx =
            (int **)calloc(hyper_query->GetNumHyperedges(), sizeof(int *));
        int *initial_cand_idx_1 = (int *)calloc(VqVg, sizeof(int));

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            int e_label = hyper_query->GetHyperedgeLabel(e);
            int max_num_of_cand =
                hyper_data->GetHyperedgesByLabel(e_label).size();
            cand[e] = cand_1;
            cand_1 += max_num_of_cand;
            cand_idx[e] = cand_idx_1;
            cand_idx_1 += hyper_data->GetNumHyperedges();
            initial_cand_idx[e] = initial_cand_idx_1;
            initial_cand_idx_1 += hyper_data->GetNumHyperedges();
        }

        global_match_of_v =
            (int *)calloc(hyper_data->GetNumHyperedges(), sizeof(int));
        memset(global_match_of_v, -1,
               sizeof(int) * hyper_data->GetNumHyperedges());

        in_queue_removed =
            (bool **)calloc(hyper_query->GetNumHyperedges(), sizeof(bool *));
        bool *in_queue_removed_1 = (bool *)calloc(
            hyper_query->GetNumHyperedges() * hyper_query->GetNumHyperedges(),
            sizeof(bool));
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            in_queue_removed[e] = in_queue_removed_1;
            in_queue_removed_1 += hyper_query->GetNumHyperedges();
        }
    }

    inline bool MaCH::BuildHyperCandidateSpace() {
        BuildInitialHCS();
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            if (cand_sz[0][e] == 0) return false;
        }
        InitializeMemory();
        BuildConnection();
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            if (cand_sz[0][e] == 0) return false;
        }
        return FilterHCS();
    }

    inline void MaCH::PushCand(int *cand, int *idx, int &sz, int v) {
        idx[v] = sz;
        cand[sz++] = v;
    }

    inline void MaCH::PushCand(int *cand, int &sz, int v) { cand[sz++] = v; }

    inline int MaCH::SwapAndPop(int *cand, int *idx, int &sz, int pos) {
        sz--;
        if (sz == pos)
            return sz;
        else {
            std::swap(cand[sz], cand[pos]);
            std::swap(idx[cand[sz]], idx[cand[pos]]);
            return sz;
        }
    }

    inline int MaCH::SwapAndPop(int *cand, int &sz, int pos) {
        sz--;
        if (sz == pos)
            return sz;
        else {
            std::swap(cand[sz], cand[pos]);
            return sz;
        }
    }

    inline void MaCH::EmptyQueue() {
        while (!queue_removed.empty()) {
            auto [u, e_idx] = queue_removed.front();
            queue_removed.pop();
            in_queue_removed[u][e_idx] = 0;
        }
    }

    void MaCH::BuildInitialHCS() {
        int min_e = -1;
        int min_sz = __INT_MAX__;
        fprintf(log_to, "# Hyperedges with the same signature\n");
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            int e_label = hyper_query->GetHyperedgeLabel(e);
            int num_same_label =
                hyper_data->GetHyperedgesByLabel(e_label).size();
            fprintf(log_to, "Hyperedge %d: %d\n", e, num_same_label);
            if (min_sz > num_same_label) {
                min_sz = num_same_label;
                min_e = e;
            }
        }

        int min_e_label = hyper_query->GetHyperedgeLabel(min_e);

        for (int &f : hyper_data->GetHyperedgesByLabel(min_e_label)) {
            PushCand(cand[min_e], cand_idx[min_e], cand_sz[0][min_e], f);
        }

        std::priority_queue<std::pair<int, int>> pq;

        pq.emplace(-cand_sz[0][min_e], min_e);
        if (cand_sz[0][min_e] == 1) {
            mapping(min_e, cand[min_e][0]);
        }

        std::vector<int> visited(hyper_query->GetNumHyperedges());
        visited[min_e] = 1;
        std::vector<int> data_visited(hyper_data->GetNumHyperedges());
        std::stack<int> adj_data;
        while (!pq.empty()) {
            auto [_sz, e] = pq.top();
            pq.pop();
            for (int i = 0; i < cand_sz[0][e]; i++) {
                int f = cand[e][i];
                if (matched[e] == -1) {
                    cell_signature->AddMapping(e, f);
                }
                for (int en_idx = 0; en_idx < query_edge_adj[e].size();
                     en_idx++) {
                    int en = query_edge_adj[e][en_idx];
                    if (visited[en] == 1) {
                        continue;
                    }
                    cell_signature->AddQueryMapping(en);
                    int en_label = hyper_query->GetHyperedgeLabel(en);
                    for (int v : hyper_data->GetHyperedge(f)) {
                        if (hyper_data->GetVertexLabel(v) !=
                            query_edge_adj_label[e][en])
                            continue;
                        for (int fn : hyper_data->GetIncidentHyperedges(v)) {
                            if (hyper_data->GetHyperedgeLabel(fn) != en_label)
                                continue;
                            if (data_visited[fn]) continue;
                            if (isCandidate(en, fn, 0)) continue;
                            data_visited[fn] = true;
                            adj_data.push(fn);
                            if (cell_signature->CheckMapping(en, fn)) {
                                PushCand(cand[en], cand_idx[en], cand_sz[0][en],
                                         fn);
                            }
                        }
                    }
                    cell_signature->RemoveQueryMapping(en);
                    while (!adj_data.empty()) {
                        int fn = adj_data.top();
                        adj_data.pop();
                        data_visited[fn] = false;
                    }
                }
                if (matched[e] == -1) {
                    cell_signature->RemoveMapping(e, f);
                }
            }
            for (int en_idx = 0; en_idx < query_edge_adj[e].size(); en_idx++) {
                int en = query_edge_adj[e][en_idx];
                if (visited[en] == 1) continue;
                visited[en] = 1;
                if (cand_sz[0][en] == 1) {
                    cell_signature->AddQueryMapping(en);
                    bool valid = cell_signature->CheckMapping(en, cand[en][0]);
                    cell_signature->RemoveQueryMapping(en);
                    if (valid) mapping(en, cand[en][0]);
                    pq.emplace(0, en);
                } else
                    pq.emplace(-cand_sz[0][en], en);
            }
        }

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            if (matched[e] != -1) continue;
            cell_signature->AddQueryMapping(e);
            for (int i = cand_sz[0][e] - 1; i >= 0; i--) {
                if (!cell_signature->CheckMapping(e, cand[e][i])) {
                    SwapAndPop(cand[e], cand_idx[e], cand_sz[0][e], i);
                }
            }
            cell_signature->RemoveQueryMapping(e);
        }
    }

    void MaCH::InitializeMemory() {
        memcpy(initial_cand_idx[0], cand_idx[0], VqVg * sizeof(int));
        MaxStage = 1;
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            initial_cand_sz_64[e] = (cand_sz[0][e] - 1) / 64 + 1;
            MaxStage += 1;
        }

        VqCu = 0;
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            VqCu_64 += initial_cand_sz_64[e];
            VqCu += cand_sz[0][e];
        }

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            for (int en_idx = 0; en_idx < query_edge_adj[e].size(); en_idx++) {
                int en = query_edge_adj[e][en_idx];
                VqCuDuCun_64 += (long long) cand_sz[0][e] * initial_cand_sz_64[en];
                VqCuDu += cand_sz[0][e];
            }
        }

        local_match =
            (int ***)calloc(hyper_query->GetNumHyperedges(), sizeof(int **));
        int **local_match_1 = (int **)calloc(VqCu, sizeof(int *));
        int *local_match_2 = (int *)calloc(VqCuDu, sizeof(int));
        memset(local_match_2, 0, VqCuDu * sizeof(int));

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            local_match[e] = local_match_1;
            local_match_1 += cand_sz[0][e];
            for (int f = 0; f < cand_sz[0][e]; f++) {
                local_match[e][f] = local_match_2;
                local_match_2 += query_edge_adj[e].size();
            }
        }

        cand_nbr_sz =
            (int ***)calloc(hyper_query->GetNumHyperedges(), sizeof(int **));
        int **cand_nbr_sz_1 = (int **)calloc(VqCu, sizeof(int *));
        int *cand_nbr_sz_2 = (int *)calloc(VqCuDu, sizeof(int));
        memset(cand_nbr_sz_2, 0, VqCuDu * sizeof(int));

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            cand_nbr_sz[e] = cand_nbr_sz_1;
            cand_nbr_sz_1 += cand_sz[0][e];
            for (int f = 0; f < cand_sz[0][e]; f++) {
                cand_nbr_sz[e][f] = cand_nbr_sz_2;
                cand_nbr_sz_2 += query_edge_adj[e].size();
            }
        }

        cand_exist = (unsigned long long ***)calloc(
            MaxStage, sizeof(unsigned long long **));
        unsigned long long **cand_exist_1 = (unsigned long long **)calloc(
            MaxStage * hyper_query->GetNumHyperedges(),
            sizeof(unsigned long long *));
        unsigned long long *cand_exist_2 = (unsigned long long *)calloc(
            MaxStage * VqCu_64, sizeof(unsigned long long));

        for (int s = 0; s < MaxStage; s++) {
            cand_exist[s] = cand_exist_1;
            cand_exist_1 += hyper_query->GetNumHyperedges();
            for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
                cand_exist[s][e] = cand_exist_2;
                cand_exist_2 += initial_cand_sz_64[e];
            }
        }

        cand_nbr_exist = (unsigned long long ****)calloc(
            hyper_query->GetNumHyperedges(), sizeof(unsigned long long ***));
        unsigned long long ***cand_nbr_exist_1 =
            (unsigned long long ***)calloc(VqCu, sizeof(unsigned long long **));
        unsigned long long **cand_nbr_exist_2 =
            (unsigned long long **)calloc(VqCuDu, sizeof(unsigned long long *));
        unsigned long long *cand_nbr_exist_3 = (unsigned long long *)calloc(
            VqCuDuCun_64, sizeof(unsigned long long));

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            cand_nbr_exist[e] = cand_nbr_exist_1;
            cand_nbr_exist_1 += cand_sz[0][e];
            for (int f = 0; f < cand_sz[0][e]; f++) {
                cand_nbr_exist[e][f] = cand_nbr_exist_2;
                cand_nbr_exist_2 += query_edge_adj[e].size();
                for (int en_idx = 0; en_idx < query_edge_adj[e].size();
                     en_idx++) {
                    int en = query_edge_adj[e][en_idx];
                    cand_nbr_exist[e][f][en_idx] = cand_nbr_exist_3;
                    cand_nbr_exist_3 += initial_cand_sz_64[en];
                }
            }
        }
    }

    void MaCH::BuildConnection() {
        int min_e = -1;
        int min_sz = __INT_MAX__;
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            if (cand_sz[0][e] < min_sz) {
                min_e = e;
                min_sz = cand_sz[0][e];
            }
            for (int f_idx = 0; f_idx < cand_sz[0][e]; f_idx++) {
                cand_exist[0][e][f_idx / 64] |= (1ULL << (f_idx % 64));
            }
        }
        std::vector<int> data_visited(hyper_data->GetNumHyperedges());
        std::stack<int> adj_data;

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            for (int en : query_edge_adj[e]) {
                if (e > en) continue;
                int e_idx = query_edge_adj_idx[en][e];
                int en_idx = query_edge_adj_idx[e][en];
                if (matched[e] != -1 && matched[en] != -1) {
                    cand_nbr_exist[e][0][en_idx][0] |= 1;
                    cand_nbr_exist[en][0][e_idx][0] |= 1;
                    cand_nbr_sz[e][0][en_idx] = 1;
                    cand_nbr_sz[en][0][e_idx] = 1;
                } else {
                    bool swapped = false;
                    if (matched[en] != -1) {
                        std::swap(e, en);
                        std::swap(e_idx, en_idx);
                        swapped = true;
                    }
                    for (int f_idx = 0; f_idx < cand_sz[0][e]; f_idx++) {
                        int f = cand[e][f_idx];
                        if (matched[e] == -1) cell_signature->AddMapping(e, f);
                        cell_signature->AddQueryMapping(en);
                        for (int v : hyper_data->GetHyperedge(f)) {
                            if (hyper_data->GetVertexLabel(v) !=
                                query_edge_adj_label[e][en])
                                continue;
                            for (int fn :
                                 hyper_data->GetIncidentHyperedges(v)) {
                                if (!isCandidate(en, fn, 0)) continue;
                                int fn_idx = initial_cand_idx[en][fn];
                                if (data_visited[fn]) continue;
                                data_visited[fn] = true;
                                adj_data.push(fn);
                                if (cell_signature->CheckMapping(en, fn)) {
                                    int e_idx = query_edge_adj_idx[en][e];
                                    cand_nbr_exist[e][f_idx][en_idx]
                                                  [fn_idx / 64] |=
                                        (1ULL << (fn_idx % 64));
                                    cand_nbr_exist[en][fn_idx][e_idx]
                                                  [f_idx / 64] |=
                                        (1ULL << (f_idx % 64));
                                    cand_nbr_sz[e][f_idx][en_idx]++;
                                    cand_nbr_sz[en][fn_idx][e_idx]++;
                                }
                            }
                        }

                        while (!adj_data.empty()) {
                            int fn = adj_data.top();
                            adj_data.pop();
                            data_visited[fn] = false;
                        }
                        cell_signature->RemoveQueryMapping(en);
                        if (matched[e] == -1)
                            cell_signature->RemoveMapping(e, f);
                    }
                    if (swapped) {
                        std::swap(e, en);
                        std::swap(e_idx, en_idx);
                    }
                }
            }
        }
        long long total_conns = 0;
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            for (int f_idx = 0; f_idx < cand_sz[0][e]; f_idx++) {
                for (int en_idx = 0; en_idx < query_edge_adj[e].size();
                     en_idx++) {
                    total_conns += cand_nbr_sz[e][f_idx][en_idx];
                }
            }
        }
        cand_nbr =
            (int ****)calloc(hyper_query->GetNumHyperedges(), sizeof(int ***));
        int ***cand_nbr_1 = (int ***)calloc(VqCu, sizeof(int **));
        int **cand_nbr_2 = (int **)calloc(VqCuDu, sizeof(int *));
        int *cand_nbr_3 = (int *)calloc(total_conns, sizeof(int));

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            cand_nbr[e] = cand_nbr_1;
            cand_nbr_1 += cand_sz[0][e];
            for (int f = 0; f < cand_sz[0][e]; f++) {
                cand_nbr[e][f] = cand_nbr_2;
                cand_nbr_2 += query_edge_adj[e].size();
                for (int en_idx = 0; en_idx < query_edge_adj[e].size();
                     en_idx++) {
                    cand_nbr[e][f][en_idx] = cand_nbr_3;
                    cand_nbr_3 += cand_nbr_sz[e][f][en_idx];
                }
            }
        }
        memset(cand_nbr_sz[0][0], 0, VqCuDu * sizeof(int));

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            for (int f = 0; f < cand_sz[0][e]; f++) {
                for (int en_idx = 0; en_idx < query_edge_adj[e].size();
                     en_idx++) {
                    int en = query_edge_adj[e][en_idx];
                    local_match[e][f][en_idx] = initial_cand_sz_64[en];
                    for (int fn_idx_64 = initial_cand_sz_64[en] - 1;
                         fn_idx_64 >= 0; fn_idx_64--) {
                        unsigned long long temp =
                            cand_exist[0][en][fn_idx_64] &
                            cand_nbr_exist[e][f][en_idx][fn_idx_64];
                        if (temp) {
                            local_match[e][f][en_idx] = fn_idx_64;
                            for (int i = 0; i < 64; i++) {
                                if (temp & (1ULL << i)) {
                                    int t = cand[en][fn_idx_64 * 64 + i];
                                    cand_nbr[e][f][en_idx]
                                            [cand_nbr_sz[e][f][en_idx]++] = t;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            for (int en : query_edge_adj[e]) {
                int e_idx = query_edge_adj_idx[en][e];
                if (!in_queue_removed[en][e_idx]) {
                    in_queue_removed[en][e_idx] = true;
                    queue_removed.emplace(en, e_idx);
                }
            }
        }
    }

    bool MaCH::RestoreLocalMatch(int e, int f, int en_idx, int stage) {
        int en = query_edge_adj[e][en_idx];
        int f_idx = initial_cand_idx[e][f];
        int &fn_64 = local_match[e][f_idx][en_idx];
        while (fn_64 < initial_cand_sz_64[en]) {
            num_restore_local_match++;
            if (cand_nbr_exist[e][f_idx][en_idx][fn_64] &
                cand_exist[stage][en][fn_64]) {
                return true;
            }
            fn_64++;
        }
        return false;
    }

    bool MaCH::FilterHCS(int stage) {
        while (!queue_removed.empty()) {
            auto [e, en_idx] = queue_removed.front();
            queue_removed.pop();
            int en = query_edge_adj[e][en_idx];

            if (matched[e] != -1) continue;
            bool filtered = false;
            for (int i = cand_sz[stage][e] - 1; i >= 0; i--) {
                int f = cand[e][i];
                int f_idx = initial_cand_idx[e][f];
                if (!RestoreLocalMatch(e, f, en_idx, stage)) {
                    conflict++;
                    cand_exist[stage][e][f_idx / 64] -= (1ULL << (f_idx % 64));
                    if (SwapAndPop(cand[e], cand_idx[e], cand_sz[stage][e],
                                   i) == 0) {
                        EmptyQueue();
                        return false;
                    }
                    filtered = true;
                }
            }
            if (filtered) {
                for (int en2 : query_edge_adj[e]) {
                    if (en2 == en) continue;
                    int e_idx = query_edge_adj_idx[en2][e];
                    if (!in_queue_removed[en2][e_idx]) {
                        in_queue_removed[en2][e_idx] = true;
                        queue_removed.emplace(en2, e_idx);
                    }
                }
            }
            in_queue_removed[e][en_idx] = false;
        }
        return true;
    }

    bool MaCH::backtracking() {
        auto print_embedding = [&](int stage, int last_e = -1) {
            if (last_e == -1) {
                fprintf(log_to, "  Embedding Found {");
                for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
                    fprintf(log_to, "(%d, %d)%c", e,
                        hyper_data->GetOrigHyperedgeId(matched[e]),
                            " }"[e == hyper_query->GetNumHyperedges() - 1]);
                }
                fprintf(log_to, "\n");
            } else {
                for (int i = 0; i < cand_sz[stage][last_e]; i++) {
                    fprintf(log_to, "  Embedding Found {");
                    for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
                        int f = e == last_e ? cand[e][i] : matched[e];
                        fprintf(log_to, "(%d, %d)%c", e,
                                hyper_data->GetOrigHyperedgeId(f),
                                " }"[e == hyper_query->GetNumHyperedges() - 1]);
                    }
                    fprintf(log_to, "\n");
                }
            }
        };

        std::vector<std::vector<int>> adj_list(
            hyper_query->GetNumHyperedges(),
            std::vector<int>(hyper_query->GetNumHyperedges(), -1));
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            for (int j = 0; j < query_edge_adj[e].size(); j++) {
                int en = query_edge_adj[e][j];
                adj_list[e][en] = j;
            }
        }
        int initial_stage = 0;

        std::stack<std::array<int, 4>> something;
        something.push({initial_stage, 0, 0,
                        hyper_query->GetNumHyperedges() - unmapped_edges_sz});

        while (!something.empty()) {
            if (++traversed_node % 1'000'000 == 0) {
                fprintf(log_to, "Traversed Nodes: %lld, Embeddings: %lld\n",
                        traversed_node, num_embedding);
            }
            auto [stage, divided_e, sz, num_matched] = something.top();
            something.pop();

            while (matched_order.size() > num_matched) {
                unmapping_last();
            }
            bool valid = true;
            if (!(stage == 0 && sz == 0)) {
                memcpy(cand_sz[stage], cand_sz[stage - 1],
                       sizeof(int) * hyper_query->GetNumHyperedges());

                cand_sz[stage][divided_e] = 1;
                int f = cand[divided_e][0];
                int f_idx = initial_cand_idx[divided_e][f];

                for (int en_idx = 0; en_idx < query_edge_adj[divided_e].size();
                     en_idx++) {
                    int en = query_edge_adj[divided_e][en_idx];
                    if (matched[en] != -1) continue;

                    if (cand_sz[stage][en] <
                        cand_nbr_sz[divided_e][f_idx][en_idx] / 2) {
                        for (int i = cand_sz[stage][en] - 1; i >= 0; i--) {
                            num_check_cleaning++;
                            int fn = cand[en][i];
                            int fn_idx = initial_cand_idx[en][fn];
                            int e_idx = query_edge_adj_idx[en][divided_e];
                            if (!(cand_nbr_exist[en][fn_idx][e_idx]
                                                [f_idx / 64] &
                                  (1ULL << (f_idx % 64)))) {
                                SwapAndPop(cand[en], cand_idx[en],
                                           cand_sz[stage][en], i);
                            }
                        }
                    } else {
                        int sz1 = 0;
                        for (int i = 0;
                             i < cand_nbr_sz[divided_e][f_idx][en_idx]; i++) {
                            num_check_cleaning++;
                            int fn = cand_nbr[divided_e][f_idx][en_idx][i];

                            if (isCandidate(en, fn, stage)) {
                                int idx = cand_idx[en][fn];
                                if (idx != sz1) {
                                    std::swap(cand[en][sz1], cand[en][idx]);
                                    std::swap(cand_idx[en][cand[en][sz1]],
                                              cand_idx[en][cand[en][idx]]);
                                }
                                sz1++;
                            }
                        }
                        cand_sz[stage][en] = sz1;
                    }
                    if (cand_sz[stage][en] == 0) {
                        valid = false;
                        break;
                    }
                }
            }

            if (valid) {
                if (stage != 0) {
                    mapping(divided_e, cand[divided_e][0]);
                    for (int i = 0; i < unmapped_edges_sz; i++) {
                        int en = unmapped_edges[i];
                        cell_signature->AddQueryMapping(en);
                        for (int j = cand_sz[stage][en] - 1; j >= 0; j--) {
                            int fn = cand[en][j];
                            bool poss = cell_signature->CheckMapping(en, fn);
                            num_check_intersection++;
                            if (!poss) {
                                conflict2++;
                                int n = SwapAndPop(cand[en], cand_idx[en],
                                                   cand_sz[stage][en],
                                                   cand_idx[en][fn]);
                                if (n == 0) {
                                    EmptyQueue();
                                    valid = false;
                                    break;
                                }
                            }
                        }
                        cell_signature->RemoveQueryMapping(en);

                        if (!valid) {
                            break;
                        }
                    }

                    if (!valid) {
                        if (stage != 0) {
                            SwapAndPop(cand[divided_e], cand_idx[divided_e],
                                       cand_sz[stage - 1][divided_e], 0);
                        }
                        continue;
                    }
                }

                long long now_sz = __LONG_LONG_MAX__;
                int now_e = -1;
                long long now_deg = -1;

                for (int i = 0; i < unmapped_edges_sz; i++) {
                    int e = unmapped_edges[i];
                    if (std::pair<int, long long>{-unmapped_deg[e],
                                                  cand_sz[stage][e]} <
                        std::pair<int, long long>{-now_deg, now_sz}) {
                        now_sz = cand_sz[stage][e];
                        now_deg = unmapped_deg[e];
                        now_e = e;
                    }
                }
                for (int i = 0; i < unmapped_edges_sz; i++) {
                    int e = unmapped_edges[i];
                    if (cand_sz[stage][e] == 1) {
                        now_e = e;
                    }
                }

                if (unmapped_edges_sz == 0) {
                    num_embedding++;
                    if (opt.print_answer) {
                        print_embedding(stage, -1);
                    }
                } else if (unmapped_edges_sz == 1) {
                    num_embedding += cand_sz[stage][now_e];
                    if (opt.print_answer and cand_sz[stage][now_e]) {
                        print_embedding(stage, now_e);
                    }
                } else {
                    for (int i = 0; i < cand_sz[stage][now_e]; i++) {
                        something.push({stage + 1, now_e, i,
                                        hyper_query->GetNumHyperedges() -
                                            unmapped_edges_sz});
                    }
                }
            }
            if (stage != 0) {
                SwapAndPop(cand[divided_e], cand_idx[divided_e],
                           cand_sz[stage - 1][divided_e], 0);
            }
        }
        return true;
    }
}  // namespace GraphLib::SubHyperGraphMatching
