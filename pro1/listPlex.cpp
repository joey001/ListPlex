int k;
thread_local long cntT=0;
int workers,workers1t2;
int taskCnt = 0,overCnt=0;
//statistics for schedule
#include <cmath>
#include "listPlex.h"
#if defined(TIMER)
double decomposeTime(0);
double buildTime(0);
double listTime(0);
#endif

namespace ListPlex
{
    struct PlexEmitor{
        //only count the number
        void emitPlex(){
            cntT++;
        }
    };

    struct DecomposeMaxChecker
    {
        const graph<int> g;
        const graph<int> subg; //check in the subgraph(contains all possible nodes and edges) for cache-efficiency
        const vector<int> &blk;
        const vector<uint8_t> &inBlk;
        vector<uint8_t> &visited;
        vector<int> &uni;

        DecomposeMaxChecker(const graph<int> &_g,const graph<int> &_subg,const vector<int>& _blk,vector<uint8_t> &_visited,vector<int> &_uni,const vector<uint8_t> &_inBlk)
            : g(_g),subg(_subg),blk(_blk),visited(_visited),uni(_uni),inBlk(_inBlk){}
        
        bool isAdjBinD2K(const int ele,const int u){
            return binary_search(subg.V[ele].NeighborsLeft,subg.V[ele].NeighborsLeft+subg.V[ele].degreeLeft,u);
        }

        //Count the number of adjacent vertices of u in the plex
        bool canGlobalAddD2K(int u,int* plex,int szplex){
            int cnt = 0;
            for (int j = 0; j < szplex;++j)
            {
                if(isAdjBinD2K(plex[j],u)){
                    cnt++;
                }
            }
            if (cnt + k >= szplex + 1)
                return true; //can add
            return false;	//can not add
        }

        bool isAdjBin(const int ele,const int u){
            return g.isAdj(ele,u);
        }

        //Count the number of adjacent vertices of u in the plex
        bool canGlobalAdd(int u,int* plex,int szplex){
            int cnt = 0;
            for (int j = 0; j < szplex;++j)
            {
                const int ele=blk[plex[j]];
                if(isAdjBin(ele,u)){
                    cnt++;
                }
            }
            if (cnt + k >= szplex + 1)
                return true; //can add
            return false;	//can not add
        }
        
        //Check if the vertices in list in a globaly maximal k-plex
        //i.e., maximality : all vertices in Preds cannot form a larger k-plex with list.
        bool isMaximal(int *list, int szplex, int* neiInP)
        {   
            if(szplex<k)return false;
            // saturated vertices;
            assert(szplex < PSIZE_LIMIT);
            int sat[PSIZE_LIMIT]; //saturated vertices in plex
            int nsat = 0;
            bool max = true;
            for (int i = 0; i < szplex; i++){
                const int v = list[i];
                if (neiInP[v] + k == szplex){
                    sat[nsat++] = v;
                }
            }
            const bool isD2K=(szplex+1>2*k-2);
            if(isD2K){
                if (nsat){
                    const int v0 = sat[0];
                    int * const commons = new int[subg.V[v0].degreeLeft];                             
                    copy(subg.V[v0].NeighborsLeft, subg.V[v0].NeighborsLeft + subg.V[v0].degreeLeft, commons);
                    int szcom = subg.V[v0].degreeLeft;
                    for (int i = 1; i < nsat; i++){
                        const int u = sat[i];
                        szcom = utils::interSection(subg.V[u].NeighborsLeft, subg.V[u].degreeLeft, commons, szcom, commons);                    
                    }
                    // (commons,szcom) contains all common left-adajcent vertices of saturated vertices                
                    for (int i = 0; i < szcom; i++){
                        const int u = commons[i];
                        if(canGlobalAddD2K(u,list,szplex)){ // not maximal because we can add u to the solution
                            max = false;
                            break;
                        }
                    }
                    delete[] commons;
                }
                else{//no saturated vertex
                    for (int i = 0; i < k;++i){
                        const int u = list[i];
                        for (int j = 0; j < subg.V[u].degreeLeft;++j){
                            const int ele = subg.V[u].NeighborsLeft[j];
                            if(!visited[ele]){
                                uni.push_back(ele);
                                visited[ele]=true;
                                if (canGlobalAddD2K(ele,list,szplex)){
                                    max = false;
                                    break;
                                }
                            }
                        }
                        if(!max)break;
                    }
                    for(const auto  u : uni){
                        visited[u]=false;
                    }
                    uni.clear();
                }
            }
            else{
                if (nsat){
                    const int v0 = blk[sat[0]];
                    int * const commons = new int[g.V[v0].degree];                             
                    copy(g.V[v0].Neighbors, g.V[v0].Neighbors + g.V[v0].degree, commons);
                    int szcom = g.V[v0].degree;
                    for (int i = 1; i < nsat; i++){
                        const int u = blk[sat[i]];
                        szcom = utils::interSection(g.V[u].Neighbors, g.V[u].degree, commons, szcom, commons);                    
                    }
                    for (int i = 0; i < szcom; i++){
                        const int u = commons[i];
                        if(inBlk[u])continue;
                        if(canGlobalAdd(u,list,szplex)){ // not maximal because we can add u to the solution
                            max = false;
                            break;
                        }
                    }
                    delete[] commons;
                }
                else{//no saturated vertex
                    for (int i = 0; i < k;++i){
                        const int u = blk[list[i]];
                        for (int j = 0; j < g.V[u].degree;++j){
                            const int ele = g.V[u].Neighbors[j];
                            if(inBlk[ele])continue;
                            if(!visited[ele]){
                                uni.push_back(ele);
                                visited[ele]=true;
                                if (canGlobalAdd(ele,list,szplex)){
                                    max = false;
                                    break;
                                }
                            }
                        }
                        if(!max)break;
                    }
                    for(const auto u : uni){
                        visited[u]=false;
                    }
                    uni.clear();
                }
            }
            return max;
		}
    };

    //Decompose a local graph contains vstart and its 1a2hop-neighbors from the global graph
    template<typename Pred>
    int decompose(const Graph &g,const int vstart,vector<int> &neibors_1a2hop,vector<int> &neibors_1a2hopLeft,Pred &&pred,vector<uint8_t> &visited){
        //collect 1 hop neigbors
        visited[vstart]  = 1;
        for(int j = 0; j < g.V[vstart].degree; j++){
            const int nei = g.V[vstart].Neighbors[j];                    
            if (!visited[nei]){
                visited[nei] = 1;
                if(pred(vstart,nei))//1-hop right-neighbor
                    neibors_1a2hop.emplace_back(nei);
                else//1-hop left-neighbor
                    neibors_1a2hopLeft.emplace_back(nei);
            }
        }
        const int sz=neibors_1a2hop.size();
        neibors_1a2hop.emplace_back(vstart);
        swap(neibors_1a2hop[0], neibors_1a2hop[sz]);
        //In fact, neibors_1a2hop contains the vstart at the first
        //collect 2-hop neigbors
        for(int i = 1; i <= sz; ++i){
            const int nei1=neibors_1a2hop[i];
            for (int j = 0; j < g.V[nei1].degree; j++)
            {
                const int nei2 = g.V[nei1].Neighbors[j];
                if (!visited[nei2]){
                    visited[nei2] = 1;
                    if(pred(vstart, nei2))
                        neibors_1a2hop.emplace_back(nei2);
                    else
                        neibors_1a2hopLeft.emplace_back(nei2);
                }
            }
        }
        //restore the global variables
        for (auto v : neibors_1a2hop ){
            visited[v] = 0;
        }
        for (auto v : neibors_1a2hopLeft){
            visited[v] = 0;
        }
        return sz+1;
    }

    //Build edges among the local nodes 
    Graph buildSubgraph(const Graph &g,const vector<int>& blk,const vector<int>& left,const vector<uint8_t> &inBlk, const vector<uint8_t> &inLeft,const int hopSz,vector<uint8_t> &isNei,vector<uint8_t> &commonMtx){  
        const int blkSz = blk.size();
        const int leftSz = left.size();
        vertex<int> *vertices = newA(vertex<int>, blkSz);
        for (int i = 0; i < blkSz; i++){
            const int origin = blk[i];
            int ne = 0, neLeft = 0;
            for (int j = 0; j < g.V[origin].degree; j++){
                const int nei = g.V[origin].Neighbors[j];
                if (inBlk[nei]){ // in block
                    ne++;
                    isNei[nei]=true;
                }
                if (inLeft[nei]){
                    neLeft++;
                    isNei[nei]=true;
                }
            }
            vertices[i].Neighbors = newA(int,ne);
            vertices[i].NeighborsLeft = newA(int,neLeft);
            vertices[i].degree=ne;
            vertices[i].degreeLeft=neLeft;

            //add right edges      
            int cnt = 0;  
            for (int j = 0; j < blkSz; j++){
                const int nei = blk[j];
                commonMtx[i*blkSz+j]=isNei[nei];
                if (isNei[nei]){
                    vertices[i].Neighbors[cnt++] = j;
                    isNei[nei]=false;
                }
                if(j==hopSz-1)vertices[i].degreeHop=cnt;
            }
            //add left edges
            cnt = 0;
            for (int j = 0; j < leftSz; j++){
                const int nei = left[j];
                if (isNei[nei]){
                    vertices[i].NeighborsLeft[cnt++] = j;
                    isNei[nei]=false;
                }
            }
        }
        return Graph(vertices, blkSz, 0, nullptr);
    }

    //Local decompose --> Tree Search
    void decomposableSearch(const Graph &g){
        int *dseq = new int[g.n];
        int *forwardNei = new int[g.n];
        int *dpos = new int[g.n];
        auto pred = [dpos](int v1 ,int v2){
            return dpos[v1] < dpos[v2];
        };
        int validblk=0;
        long cntMaxPlex=0;
        Timer watch;
        watch.start();
        degeneracyOrder<int>(g, dseq, forwardNei, dpos);        
        #pragma omp parallel
        {
            thread_local vector<int> uni;
            thread_local vector<uint8_t> visited(g.n,0);
            thread_local Stack<int> plex(PSIZE_LIMIT);
            thread_local VtxSet cand1(g.n);
            thread_local VtxSet cand2(g.n);
            thread_local VtxSet candO(g.n);
            thread_local VtxSet excl(g.n);
            thread_local VtxSet exclBK(g.n);
            thread_local VtxSet exclSM(g.n);
            thread_local Stack<int> exclStack(g.n);
            thread_local vector<int> neiInP(g.n,0);
            thread_local vector<int> neiInG(g.n);
            thread_local vector<uint8_t> commonMtx(BLK_LIMIT*BLK_LIMIT);
            thread_local vector<int> blk;//vstart and its right 1a2hop-neighbors
            thread_local vector<int> left;//vstart's left 1a2hop-neighbors
            thread_local vector<uint8_t> inBlk(g.n,0);
            thread_local vector<uint8_t> inLeft(g.n,0);

            #pragma omp for schedule(dynamic) reduction(+:validblk)
            for(int idx = 0; idx < g.n; idx++){
                //updata task statistics
                #pragma omp atomic
                taskCnt++;
                // Search with leading vertex vstart;
				const int vstart = dseq[idx];
 
                // decompose the graph with leading vertex vstart
                blk.clear();
                left.clear();
                const int hopSz=decompose(g,vstart,blk,left,pred,visited);
                const int blkSz=blk.size();
                const int hop2Sz=blkSz-hopSz;
                const int leftSz=left.size();
                validblk++;
                //build subgraph 
                for(const auto u: blk)
                    inBlk[u] = true;
                for(const auto u: left)
                    inLeft[u] = true;
                Graph subg = buildSubgraph(g,blk,left,inBlk,inLeft,hopSz,visited,commonMtx);

                //build search context
                plex.resize(hopSz+k);
                cand1.resize(hopSz);
                cand2.resize(hop2Sz);
                excl.resize(blkSz);
                exclBK.resize(blkSz);
                for (int j = 0; j < subg.n; j++){
                    neiInG[j] = subg.V[j].degreeHop;
                    if(j<hopSz)cand1.add(j);
                    else cand2.add(j);
                }
                DecomposeMaxChecker maxchecker(g,subg,blk,visited,uni,inBlk);
                PlexEmitor emitor;
                KplexListor<DecomposeMaxChecker, PlexEmitor> listor(
                    subg, &maxchecker, &emitor, &plex, &cand1, &cand2, &excl, &exclBK, &exclStack,
                    neiInP.data(), neiInG.data(), hopSz, hop2Sz, commonMtx);

                listor.cand1ToPlex(0)->kSearch(k-1);//start from the first vertex
                listor.plexToCand1();
                cand1.clear();
                cand2.clear();
                subg.del();

                plex.resize(2*k-2);
                cand1.resize(g.n);
                cand2.resize(g.n);
                candO.resize(g.n);
                SmallKplexListor<PlexEmitor> smallListor(
                    g, &emitor, &plex, &cand1, &cand2, &candO, &exclSM, &exclStack,
                    neiInP.data());
                const int beg=exclSM.sz;
                for(int j=beg;j<idx;++j){
                    const int ele=dseq[j];
                    exclSM.add(ele);
                }
                for(int j=0;j<hopSz;++j){
                    const int ele=blk[j];
                    cand1.add(ele);
                }
                for(int j=hopSz;j<blkSz;++j){
                    const int ele=blk[j];
                    cand2.add(ele);
                }
                for(int j=idx;j<g.n;++j){
                    const int ele=dseq[j];
                    if(!inBlk[ele])candO.add(ele);
                }
                smallListor.cand1ToPlex(vstart)->listO(k-1);//start from the first vertex
                smallListor.plexToCand1();
                cand1.clear();
                cand2.clear();
                candO.clear();
                excl.clear();
                for(auto u: blk)
                    inBlk[u] = false;
                for(auto u: left)
                    inLeft[u] = false;
            }
            //updata task statistics
            #pragma omp atomic
            overCnt++;
        }
        #pragma omp parallel reduction(+:cntMaxPlex)
        {
            cntMaxPlex+=cntT;
        }
        double durat = watch.stop();
        printf("Total valid block %d\n", validblk);
        printf("Number of %d-plex: %ld\n", k, cntMaxPlex);
        printf("Running time %.3f\n",durat);
        printf("Task Cnt: %d\n", taskCnt);
        delete[] dseq;
        delete[] forwardNei;
        delete[] dpos;
    }
}
int main(int argc, char **argv)
{
    if (argc == 3 || argc == 4){       
        k = atoi(argv[2]);
        if(argc==4)
            setWorkers(workers = atoi(argv[3]));
        else 
            workers = getWorkers();
        workers1t2 = workers >> 1;
        graph<intT> g = ListPlex::readBinaryGraph(argv[1]); 
        ListPlex::decomposableSearch(g);
        g.del();
    }else {
        fprintf(stderr, "usage: listPlex <filename> <k> [worker-num]\n");
    }
    return 0;
}