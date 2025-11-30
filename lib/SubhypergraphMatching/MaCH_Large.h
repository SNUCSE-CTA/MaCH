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
        
        std::vector<int> data_visited;
        std::stack<int> adj_data;

        std::vector<std::vector<int>> data_visited2;
        std::vector<std::stack<int>> adj_data2;

        int MaxStage;

        long long *checksum;

        int ** cand;

        int ** cand_idx;

        int ** cand_sz;

        unsigned long long ***cand_exist;

        long long VqVg = 0, VqVgDu = 0, VqVgDuVg = 0;
        long long Vg_64 = 0, VqVg_64 = 0, VqVgDuVg_64 = 0;
        long long VqCu = 0, VqCuDu = 0;
        long long VqCu_64 = 0, VqCuDuCun_64 = 0;


        int***cand_nbr_sz;

        int ****cand_nbr;



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
        GraphLib::SubHyperGraphMatching::PatternHyperGraph *hyper_query_, SubHyperGraphMatchingOption option);
        ~MaCH();

        inline bool isCandidate(int u, int v, int stage) const { return cand_sz[stage][u] > cand_idx[u][v] && cand_idx[u][v] >= 0; }
        void BuildInitialHCS();
        bool FilterHCS(int stage = 0);
        void GetQueryConnectivity();
        void Initialize();
        void InitializeMemory();
        bool BuildHyperCandidateSpace();
        void BuildConnection();

        void PrintCSStatistics(int level = 0, int stage = 0);

        bool RestoreLocalMatch(int u, int v, int e_idx, int stage);

        bool backtracking(int e = -1, int f = -1);


        inline void PushCand(int *cand, int *idx, int &sz, int v);
        inline void PushCand(int *cand, int &sz, int v);
        inline int SwapAndPop(int *cand, int *idx, int ***cand_nbr,  int **cand_nbr_sz, int &sz, int pos);
        inline int SwapAndPop(int *cand, int *idx, int &sz, int pos);
        inline int SwapAndPop(int *cand, int &sz, int pos);

        inline void EmptyQueue();

        void mapping(int e, int f){
            matched_order.push(e);
            cell_signature->AddMapping(e, f);
            global_match_of_v[f] = e;
            matched[e] = f;
            for(int e_: query_edge_adj[e]){
                unmapped_deg[e_]--;
            }
            SwapAndPop(unmapped_edges, unmapped_edges_idx, unmapped_edges_sz, unmapped_edges_idx[e]);
        }

        void unmapping_last(){
            int e = matched_order.top(); matched_order.pop();
            cell_signature->RemoveMapping(e, matched[e]);
            global_match_of_v[matched[e]] = -1;
            matched[e] = -1;
            for(int e_: query_edge_adj[e]){
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
        PrintCSStatistics(1,0);

        fprintf(log_to, "HCSBuildingAndFileringTime: %.02lf\n", AfterMemory.GetTime());

        traversed_node = 0;
        num_embedding = 0;
        Timer BtTimer;
        BtTimer.Start();
        if(SuccessHCS) backtracking();
        long long ans = num_embedding;

        BtTimer.Stop();
        fprintf(log_to, "SearchTime: %.02lf\n", BtTimer.GetTime());
        fprintf(log_to, "TotalTime: %.02lf\n", AfterMemory.GetTime()+BtTimer.GetTime());
        fprintf(log_to, "#Embeddings: %lld\n", ans);

    }

    MaCH::~MaCH() {}


    void MaCH::PrintCSStatistics(int level, int stage) {
        fprintf(log_to, "\e\[0;31m[CS Statistics]\n");
        int num_vertices = 0;
        for (int u = 0; u < hyper_query->GetNumHyperedges(); u++) {
            num_vertices += cand_sz[stage][u];
            if (level >= 1) {
                int label = hyper_query->GetHyperedgeLabel(u);
                fprintf(log_to, "Hyperedge %d : Matching %d Hyperedges\n",u,cand_sz[stage][u]);
                if(level>=2){
                    for(int i=0; i<cand_sz[stage][u]; i++){
                        fprintf(log_to, "%d ", cand[u][i]);
                    }
                    fprintf(log_to, "\n");
                }
            }
        }
        fprintf(log_to,"Number of Hyperedges: %d \e[0m\n", num_vertices);
    }
    
    void MaCH::GetQueryConnectivity(){
        std::vector<int> visited(hyper_query->GetNumHyperedges());
        query_edge_adj = std::vector<std::vector<int>>(hyper_query->GetNumHyperedges()); 
        query_edge_adj_idx = std::vector<std::vector<int>>(hyper_query->GetNumHyperedges(), std::vector<int>(hyper_query->GetNumHyperedges())); 
        query_edge_adj_label = std::vector<std::vector<int>>(hyper_query->GetNumHyperedges(), std::vector<int>(hyper_query->GetNumHyperedges()));
        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){
            for(int v: hyper_query->GetHyperedge(e)){
                for (int e_ : hyper_query->GetIncidentHyperedges(v)){
                    if(visited[e_] || e == e_) continue;
                    query_edge_adj[e].push_back(e_);
                    visited[e_] = 1;
                    query_edge_adj_label[e][e_] = hyper_query->GetVertexLabel(v);
                }
            }
            for(int i = 0; i < query_edge_adj[e].size(); i++){
                int e_ = query_edge_adj[e][i];
                query_edge_adj_idx[e][e_] = i;
                visited[e_] = 0;
            }
        }
    }

    void MaCH::Initialize() {

        GetQueryConnectivity();

        MaxStage = 1;
        Vg_64 = (hyper_data->GetNumHyperedges()-1)/64+1;

        matched = std::vector<int>(hyper_query->GetNumHyperedges(), -1);
        unmapped_deg = std::vector<int>(hyper_query->GetNumHyperedges());
        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){
            unmapped_deg[e] = query_edge_adj[e].size();
        }

        unmapped_edges = (int *) calloc(hyper_query->GetNumHyperedges(), sizeof(int));
        unmapped_edges_idx = (int *) calloc(hyper_query->GetNumHyperedges(), sizeof(int));

        unmapped_edges_sz = hyper_query->GetNumHyperedges();
        for(int e=0; e<hyper_query->GetNumHyperedges(); e++){
            unmapped_edges[e]=e;
            unmapped_edges_idx[e]=e;
        }


        for(int e=0; e < hyper_query->GetNumHyperedges();e++){
            int e_label = hyper_query->GetHyperedgeLabel(e);
            auto &Hyperedges_with_same_label = hyper_data->GetHyperedgesByLabel(e_label);
            int num_same_label = Hyperedges_with_same_label.size();
            MaxStage += 1;
            VqCu += num_same_label;
            VqVg_64 += Vg_64;
            for(int f = 0; f < hyper_data->GetNumHyperedges(); f++){
                VqVg++;
                VqVgDu += query_edge_adj[e].size();
                VqVgDuVg += query_edge_adj[e].size()*hyper_data->GetNumHyperedges();
                VqVgDuVg_64 += query_edge_adj[e].size()*Vg_64;
            }
        }

        hyperedge_weight = (double*) calloc(hyper_query->GetNumHyperedges(), sizeof(double));
        cand_sz = (int **) calloc(MaxStage, sizeof(int *));
        int *cand_sz_1 = (int *) calloc(MaxStage * hyper_query->GetNumHyperedges(), sizeof(int));

        for(int s=0;s<MaxStage;s++){
            cand_sz[s] = cand_sz_1;
            cand_sz_1 += hyper_query->GetNumHyperedges();
        }

        cand = (int **) calloc(hyper_query->GetNumHyperedges(), sizeof(int *));
        int *cand_1 = (int *) calloc(VqCu, sizeof(int));

        cand_idx = (int **) calloc(hyper_query->GetNumHyperedges(), sizeof(int *));
        int *cand_idx_1 = (int *) calloc(VqVg, sizeof(int));
        memset(cand_idx_1, -1, VqVg*sizeof(int));


        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            int e_label = hyper_query->GetHyperedgeLabel(e);
            int max_num_of_cand = hyper_data->GetHyperedgesByLabel(e_label).size();
            cand[e] = cand_1; cand_1 += max_num_of_cand;
            cand_idx[e] = cand_idx_1; cand_idx_1 += hyper_data->GetNumHyperedges();
        }

        global_match_of_v = (int *) calloc(hyper_data->GetNumHyperedges(), sizeof(int));
        memset(global_match_of_v, -1, sizeof(int)*hyper_data->GetNumHyperedges());


        in_queue_removed = (bool **) calloc(hyper_query->GetNumHyperedges(), sizeof(bool *));
        bool *in_queue_removed_1 = (bool *) calloc(hyper_query->GetNumHyperedges()*hyper_query->GetNumHyperedges(), sizeof(bool));
        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){
            in_queue_removed[e] = in_queue_removed_1;
            in_queue_removed_1 += hyper_query->GetNumHyperedges();
        }
    }

    inline bool MaCH::BuildHyperCandidateSpace() {
        BuildInitialHCS();
        bool ret = true;
        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){
            if(cand_sz[0][e] == 0) return false;
        }
        InitializeMemory();
        BuildConnection();
        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){
            if(cand_sz[0][e] == 0) return false;
        }
        return FilterHCS();
    }


    inline void MaCH::PushCand(int *cand, int *idx, int &sz, int v){
        idx[v] = sz;
        cand[sz++]=v;
    }


    inline void MaCH::PushCand(int *cand, int &sz, int v){
        cand[sz++]=v;
    }


    inline int MaCH::SwapAndPop(int *cand, int *idx, int ***cand_nbr,  int **cand_nbr_sz, int &sz, int pos){
        sz--;
        if(sz==pos) return sz;
        else{
            std::swap(cand[sz], cand[pos]);
            std::swap(cand_nbr[sz], cand_nbr[pos]);
            std::swap(cand_nbr_sz[sz], cand_nbr_sz[pos]);
            std::swap(idx[cand[sz]], idx[cand[pos]]);
            return sz;
        }
    }

    inline int MaCH::SwapAndPop(int *cand, int *idx, int &sz, int pos){
        sz--;
        if(sz==pos) return sz;
        else{
            std::swap(cand[sz], cand[pos]);
            std::swap(idx[cand[sz]], idx[cand[pos]]);
            return sz;
        }
    }

    inline int MaCH::SwapAndPop(int *cand, int &sz, int pos){
        sz--;
        if(sz==pos) return sz;
        else{
            std::swap(cand[sz], cand[pos]);
            return sz;
        }
    }
    
    inline void MaCH::EmptyQueue(){
        while(!queue_removed.empty()){
            auto [u, e_idx] = queue_removed.front();
            queue_removed.pop();
            in_queue_removed[u][e_idx] = 0;
        }
    }

    void MaCH::BuildInitialHCS() {
        
        int min_e = -1;
        int min_sz = __INT_MAX__;
        fprintf(log_to, "# hyperedges with the same signature\n");
        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){
            int e_label = hyper_query->GetHyperedgeLabel(e);
            int num_same_label = hyper_data->GetHyperedgesByLabel(e_label).size();
            fprintf(log_to, "Hyperedge %d: %d\n", e, num_same_label);
            if(min_sz > num_same_label){
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

        std::vector<int> visited(hyper_query->GetNumHyperedges());
        visited[min_e] = 1;
        data_visited = std::vector<int>(hyper_data->GetNumHyperedges());
        while(!pq.empty()){
            auto [_sz, e] = pq.top(); pq.pop();
            for(int i = 0; i < cand_sz[0][e]; i++){
                int f = cand[e][i];
                cell_signature->AddMapping(e, f);
                for(int en_idx = 0; en_idx < query_edge_adj[e].size(); en_idx++){
                    int en = query_edge_adj[e][en_idx];
                    if(visited[en] == 1){
                        continue;
                    }
                    cell_signature->AddQueryMapping(en);
                    int en_label = hyper_query->GetHyperedgeLabel(en);
                    for(int v: hyper_data->GetHyperedge(f)){
                        if(hyper_data->GetVertexLabel(v) != query_edge_adj_label[e][en]) continue;
                        for(int fn: hyper_data->GetIncidentHyperedges(v)){
                            if(hyper_data->GetHyperedgeLabel(fn) != en_label) continue;
                            if(data_visited[fn]) continue;
                            if(isCandidate(en, fn, 0)) continue;
                            data_visited[fn] = true;
                            adj_data.push(fn);
                            if(cell_signature->CheckMapping(en, fn)){
                                PushCand(cand[en], cand_idx[en], cand_sz[0][en], fn);
                            }
                        }
                    }
                    cell_signature->RemoveQueryMapping(en);
                    while(!adj_data.empty()){
                        int fn = adj_data.top(); adj_data.pop();
                        data_visited[fn] = false;
                    }
                }
                cell_signature->RemoveMapping(e, f);
            }
            for(int en_idx=0; en_idx < query_edge_adj[e].size(); en_idx++){
                int en = query_edge_adj[e][en_idx];
                if(visited[en] == 1) continue;
                visited[en] = 1;
                pq.emplace(-cand_sz[0][en], en);
            }
        }
        
        PrintCSStatistics(1, 0);

        fflush(stderr);

    }

    
    void MaCH::InitializeMemory() {

        VqCu = 0;
        for(int e=0; e< hyper_query->GetNumHyperedges(); e++){
            VqCu += cand_sz[0][e];
        }
        
        for(int e=0; e< hyper_query->GetNumHyperedges(); e++){
            int e_label = hyper_query->GetHyperedgeLabel(e);
            for(int en_idx = 0; en_idx < query_edge_adj[e].size(); en_idx++){
                int en = query_edge_adj[e][en_idx];
                VqCuDu += cand_sz[0][e];
            }
        }

        cand_nbr_sz = (int***) calloc(hyper_query->GetNumHyperedges(), sizeof(int**));
        for (int e = 0; e < hyper_query->GetNumHyperedges(); e++) {
            cand_nbr_sz[e] = (int**) calloc(cand_sz[0][e], sizeof(int*));
            for (int f = 0; f < cand_sz[0][e]; f++) {
                cand_nbr_sz[e][f] = (int*) calloc(query_edge_adj[e].size(), sizeof(int));
            }
        }
    }

    
    void MaCH::BuildConnection() {

        int min_e = -1;
        int min_sz = __INT_MAX__;
        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){
            if(cand_sz[0][e] < min_sz){
                min_e = e;
                min_sz = cand_sz[0][e];
            }
        }
        std::vector<int> data_visited(hyper_data->GetNumHyperedges());
        std::stack<int> adj_data;
        int maximum_adj = 0;
        cand_nbr = (int ****) calloc(hyper_query->GetNumHyperedges(), sizeof(int ***));

        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){
            cand_nbr[e] = (int ***) calloc(cand_sz[0][e], sizeof(int**));
            for(int f_idx=0; f_idx < cand_sz[0][e]; f_idx++){
                cand_nbr[e][f_idx] = (int **) malloc(query_edge_adj[e].size() * sizeof(int *));
            }
        }

        int *temp_nbr = (int *) malloc(hyper_data->GetNumHyperedges() * sizeof(int));
        int sz = 0;
        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){

            for(int f_idx = 0; f_idx < cand_sz[0][e]; f_idx++){
                // int now_adj = 0;
                int f = cand[e][f_idx];
                int sz = 0;
                cell_signature->AddMapping(e, f);
                for(int en_idx = 0; en_idx < query_edge_adj[e].size(); en_idx++){
                    int sz = 0;
                    int en = query_edge_adj[e][en_idx];
                    cell_signature->AddQueryMapping(en);
                    for(int v: hyper_data->GetHyperedge(f)){
                        if(hyper_data->GetVertexLabel(v) != query_edge_adj_label[e][en]) continue;
                        for(int fn: hyper_data->GetIncidentHyperedges(v)){
                            if(!isCandidate(en, fn, 0)) continue;
                            int fn_idx = cand_idx[en][fn];
                            if(data_visited[fn]) continue;
                            data_visited[fn] = true;
                            adj_data.push(fn);
                            if(cell_signature->CheckMapping(en, fn)){
                                temp_nbr[sz++] = fn;
                            }
                        }
                    }
                    while(!adj_data.empty()){
                        int fn = adj_data.top(); adj_data.pop();
                        data_visited[fn] = false;
                    }
                    cell_signature->RemoveQueryMapping(en);
                    cand_nbr[e][f_idx][en_idx] = (int *) malloc(sz * sizeof(int));
                    memcpy(cand_nbr[e][f_idx][en_idx], temp_nbr, sz * sizeof(int));
                    cand_nbr_sz[e][f_idx][en_idx] = sz;
                }
                cell_signature->RemoveMapping(e, f);
            }
        }

        for(int e = 0; e < hyper_query->GetNumHyperedges(); e++){
            for(int en: query_edge_adj[e]){
                int e_idx = query_edge_adj_idx[en][e];
                if(!in_queue_removed[en][e_idx]){
                    in_queue_removed[en][e_idx] = true;
                    queue_removed.emplace(en, e_idx);
                }
            }
        }
        fflush(stderr);

    }


    bool MaCH::RestoreLocalMatch(int e, int i, int en_idx, int stage){
        int en = query_edge_adj[e][en_idx];
        int f = cand[e][i];
        int &sz = cand_nbr_sz[e][i][en_idx];
        while(sz){
            int fn = cand_nbr[e][i][en_idx][sz-1];
            if(isCandidate(en, fn, stage)){
                break;
            } 
            else sz--;
        }
        return sz;
    }

    bool MaCH::FilterHCS(int stage) {
        while (!queue_removed.empty()) {
            auto [e, en_idx] = queue_removed.front(); queue_removed.pop();
            int en = query_edge_adj[e][en_idx];
            if(matched[e] != -1) continue;
            bool filtered = false;
            for(int i = cand_sz[stage][e] - 1; i >= 0; i--){
                int f = cand[e][i];
                if(!RestoreLocalMatch(e, i, en_idx, stage)){
                    conflict++;
                    if(SwapAndPop(cand[e], cand_idx[e], cand_nbr[e], cand_nbr_sz[e], cand_sz[stage][e], i)==0){
                        EmptyQueue();
                        return false;
                    }
                    filtered = true;
                }
            }
            if(filtered){
                for(int en2: query_edge_adj[e]){
                    if(en2 == en) continue;
                    int e_idx = query_edge_adj_idx[en2][e];
                    if(!in_queue_removed[en2][e_idx]){
                        in_queue_removed[en2][e_idx] = true;
                        queue_removed.emplace(en2, e_idx);
                    }
                }
            }
            in_queue_removed[e][en_idx] = false;
        }
        return true;
    }


    bool MaCH::backtracking(int e_, int f_){


        std::vector<std::vector<int>> adj_list(hyper_query->GetNumHyperedges(), std::vector<int>(hyper_query->GetNumHyperedges(), -1));
        for(int e=0; e<hyper_query->GetNumHyperedges(); e++){
            for(int j = 0; j < query_edge_adj[e].size(); j++){
                int en = query_edge_adj[e][j];
                adj_list[e][en] = j;
            } 
        }
        int initial_stage = 0;
        int initial_d = 0;

        std::stack<std::array<int, 4>> something;
        something.push({initial_stage, 0, 0, hyper_query->GetNumHyperedges() - unmapped_edges_sz});
        
        while(!something.empty()){
            if(++traversed_node%1'000'000==0){
                fprintf(log_to, "Traversed Nodes: %lld, Embeddings: %lld\n",
                    traversed_node, num_embedding);
            }
            auto [stage, divided_e, sz, num_matched] = something.top(); something.pop();
            while(matched_order.size() > num_matched){
                unmapping_last();
            }
            bool valid = true;
            if(!(stage==0 && sz==0)){
                memcpy(cand_sz[stage], cand_sz[stage-1], sizeof(int) * hyper_query->GetNumHyperedges());
                
                cand_sz[stage][divided_e] = 1;
                int f = cand[divided_e][0];
                int f_idx = cand_idx[divided_e][f];

                for(int en_idx = 0; en_idx < query_edge_adj[divided_e].size(); en_idx++){
                    int en = query_edge_adj[divided_e][en_idx];
                    if(matched[en] != -1) continue;
                    
                    if(false){
                        for(int i = cand_sz[stage][en]-1; i>=0; i--){
                            num_check_cleaning++;
                            int fn = cand[en][i];
                            int e_idx = query_edge_adj_idx[en][divided_e];
                            if(!isCandidate(en, fn, stage)){
                                SwapAndPop(cand[en], cand_idx[en], cand_nbr[en], cand_nbr_sz[en], cand_sz[stage][en], i);
                            }
                        }
                    }
                    else{
                        int sz1 = 0;
                        for(int i = 0; i < cand_nbr_sz[divided_e][f_idx][en_idx]; i++){
                            num_check_cleaning++;
                            int fn = cand_nbr[divided_e][f_idx][en_idx][i];
                            if(isCandidate(en, fn, stage)){
                                int idx = cand_idx[en][fn];
                                if(idx != sz1){
                                    std::swap(cand[en][sz1], cand[en][idx]);
                                    std::swap(cand_nbr[en][sz1], cand_nbr[en][idx]);
                                    std::swap(cand_nbr_sz[en][sz1], cand_nbr_sz[en][idx]);
                                    std::swap(cand_idx[en][cand[en][sz1]], cand_idx[en][cand[en][idx]]);
                                }
                                sz1++;
                            }
                        }
                        cand_sz[stage][en] = sz1;
                    }
                    if(cand_sz[stage][en]  == 0){
                        valid = false;
                        break;
                    }
                }
            }
            if(valid){
                if(stage!=0){
                    mapping(divided_e, cand[divided_e][0]);
                    for(int i=0; i<unmapped_edges_sz; i++){
                        int en = unmapped_edges[i];
                        cell_signature->AddQueryMapping(en);
                        for(int j = cand_sz[stage][en]-1; j >= 0; j--){
                            int fn = cand[en][j];
                            bool poss = cell_signature->CheckMapping(en, fn);
                            num_check_intersection++;
                            if(!poss){
                                conflict2++;
                                int n = SwapAndPop(cand[en], cand_idx[en], cand_nbr[en], cand_nbr_sz[en], cand_sz[stage][en], cand_idx[en][fn]);
                                if(n == 0){
                                    valid = false;
                                    break;
                                }
                            }
                        }
                        cell_signature->RemoveQueryMapping(en);
                        if(!valid){
                            break;
                        }
                    }
                
                    if(!valid){
                        if(stage!=0){
                            int f = cand[divided_e][0];
                            SwapAndPop(cand[divided_e], cand_idx[divided_e], cand_nbr[divided_e], cand_nbr_sz[divided_e], cand_sz[stage-1][divided_e], 0);
                        }
                        continue;
                    }
                }
                
                long long now_sz = __LONG_LONG_MAX__;
                int now_e = -1;
                long long now_deg = -1;

                for(int i=0; i < unmapped_edges_sz; i++){
                    int e = unmapped_edges[i];
                    if(std::pair<int, long long>{-unmapped_deg[e], cand_sz[stage][e]} < std::pair<int, long long>{-now_deg, now_sz}){
                        now_sz = cand_sz[stage][e];
                        now_deg = unmapped_deg[e];
                        now_e = e;
                    }
                }
                for(int i=0; i < unmapped_edges_sz; i++){
                    int e = unmapped_edges[i];
                    if(cand_sz[stage][e] == 1){
                        now_e = e;
                    }
                }

                if(unmapped_edges_sz == 0){
                    num_embedding++;
                }
                else if(unmapped_edges_sz == 1){
                    num_embedding += cand_sz[stage][now_e];
                }
                else{
                    for(int i=0; i < cand_sz[stage][now_e]; i++){
                        something.push({stage+1, now_e, i, hyper_query->GetNumHyperedges() - unmapped_edges_sz});
                    }
                }
            }
            if(stage!=0) {
                int f = cand[divided_e][0];
                SwapAndPop(cand[divided_e], cand_idx[divided_e], cand_nbr[divided_e], cand_nbr_sz[divided_e], cand_sz[stage-1][divided_e], 0);
            }
        }
        return true;
    }
}