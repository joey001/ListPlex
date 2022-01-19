int k,lb,bd;
thread_local long cntT=0;
int workers,workers1t2;
int taskCnt = 0,overCnt=0;
//statistics for schedule
#include <cmath>
#include "listPlex.h"
#include <queue>
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
        const graph<int> subg; //check in the subgraph(contains all possible nodes and edges) for cache-efficiency

        DecomposeMaxChecker(const graph<int> &_subg)
            : subg(_subg){}

        //Count the number of adjacent vertices of u in the plex
        bool canGlobalAdd(int u, int szplex,int* plex){
            int cnt = 0;
            for (int j = 0; j < szplex;++j)
            {
                const int ele = plex[j];
                if(binary_search(subg.V[ele].NeighborsLeft,subg.V[ele].NeighborsLeft+subg.V[ele].degreeLeft,u))
                    cnt++;
            }
            if (cnt + k >= szplex + 1)
                return true; //can add
            return false;	//can not add
        }
        
        //Check if the vertices in list in a globaly maximal k-plex
        //i.e., maximality : all vertices in Preds cannot form a larger k-plex with list.
        bool isMaximal(int *list, int szplex, int* neiInP)
        {   
            // saturated vertices;
            assert(szplex < PSIZE_LIMIT);
            int sat[PSIZE_LIMIT]; //saturated vertices in plex
            int nsat = 0;
            bool max = true;
            thread_local RandSet<int> uni(CAP_LIMIT2);
            for (int i = 0; i < szplex; i++){
                const int v = list[i];
                if (neiInP[v] + k == szplex){
                    sat[nsat++] = v;
                }
            }
            if (nsat){
                int v0 = sat[0];
                int *commons = new int[subg.V[v0].degreeLeft];                             
                copy(subg.V[v0].NeighborsLeft, subg.V[v0].NeighborsLeft + subg.V[v0].degreeLeft, commons);
                int szcom = subg.V[v0].degreeLeft;
                for (int i = 1; i < nsat; i++){
                    int u = sat[i];
                    szcom = utils::interSection(subg.V[u].NeighborsLeft, subg.V[u].degreeLeft, commons, szcom, commons);                    
                }
                //(commons,szcom) contains all common left-adajcent vertices of saturated vertices                
                if(nsat==1){
                    for (int i = 0; i < k;++i){
                        int u = list[i];
                        for (int j = 0; j < subg.V[u].degreeLeft;++j){
                            const int ele = subg.V[u].NeighborsLeft[j];
                            if(!uni.contains(ele))uni.add(ele);
                        }
                    }
                    int sz = szcom;
                    szcom = 0;
                    for (int i = 0; i < sz;++i)
                    {
                        if(uni.contains(commons[i]))
                            commons[szcom++]=commons[i];
                    }
                    uni.clearByEle();
                }
                for (int i = 0; i < szcom; i++){
                    int u = commons[i];
                    if (canGlobalAdd(u, szplex, list)){ // not maximal because we can add u to the solution
                        max = false;
                        break;
                    }
                }   
                delete[] commons;
            }             
            else{//no saturated vertex
                for (int i = 0; i < k;++i){
                    int u = list[i];
                    for (int j = 0; j < subg.V[u].degreeLeft;++j){
                        const int ele = subg.V[u].NeighborsLeft[j];
                        if(!uni.contains(ele))uni.add(ele);
                    }
                }
                for (int i = 0; i < uni.sz;++i)
                {
                    const int ele = uni.members[i];
                    //local->global->degeneracy->status
                    if (canGlobalAdd(ele, szplex,list)){
                        max = false;
                        break;
                    }
                }
                uni.clearByEle();
	    	}
	    	return max;
		}      
    };

    //Decompose a local graph contains vstart and its 1a2hop-neighbors from the global graph
    int decompose(const Graph &g,const int vstart,vector<int> &neibors_1a2hop,vector<int> &neibors_1a2hopLeft,const int* dpos,vector<uint8_t> &visited,vector<int> &count,vector<int> &neibors,vector<int> &neiborsLeft){
        //collect 1 hop neigbors
        visited[vstart]  = 1;
        for(int j = 0; j < g.V[vstart].degree; j++){
            const int nei = g.V[vstart].Neighbors[j];                    
            if (!visited[nei]){
                visited[nei] = 1;
                if(dpos[vstart]<dpos[nei])//1-hop right-neighbor
                    neibors_1a2hop.emplace_back(nei);
                else//1-hop left-neighbor
                    neibors_1a2hopLeft.emplace_back(nei);
            }
        }
        auto filter = [&](const std::vector<int>& origin,const std::vector<int>& neibors_1hop,const int thres) {
            std::vector<int> filtered;
            const int sz2 = neibors_1hop.size();
            for (const auto u : origin) {
                size_t cnt = 0;
                const int sz1 = g.V[u].degree;
                //choose proper algorithm by time complexity estimation
                if(sz2 * log2(sz1) < (sz2 + sz1)){
                    for (auto v : neibors_1hop) {
                        if (g.isAdj(u,v))cnt++;
                    }
                }
                else{
                    const int *const lst1 = g.V[u].Neighbors;
                    const int *const lst2 = neibors_1hop.data();
                    int j1 = 0, j2 = 0;
                    while (j1 < sz1 && j2 < sz2) {
                        if (lst1[j1] < lst2[j2])
                            ++j1;
                        else if (lst1[j1] > lst2[j2])++j2;
                        else {
                            cnt++,j1++, j2++;
                        }
                    }
                }
                if (cnt >= thres) filtered.push_back(u);
                else
                    visited[u] = 0;
            }
            return filtered;
        };
        size_t sz;
        const int lb_2k = lb - 2 * k;
        do {
            sz = neibors_1a2hop.size();
            if(sz<bd)goto CLEAN;
            //the number of common neighbors with vstart >= lb-2*k
            neibors_1a2hop = filter(neibors_1a2hop, neibors_1a2hop, lb_2k);
        } while (neibors_1a2hop.size() < sz);
        //the number of common neighbors with vstart >= lb+1-2*k
        neibors_1a2hopLeft = filter(neibors_1a2hopLeft, neibors_1a2hop,lb_2k+1);

        //collect 2-hop neigbors
        for(const auto nei1 : neibors_1a2hop){
            for (int j = 0; j < g.V[nei1].degree; j++)
            {
                const int nei2 = g.V[nei1].Neighbors[j];
                if (!visited[nei2]){
                    if(!count[nei2]){
                        if(dpos[vstart]<dpos[nei2])
                            neibors.emplace_back(nei2);
                        else
                            neiborsLeft.emplace_back(nei2);
                    }
                    count[nei2]++;
                }
            }
        }
        neibors_1a2hop.emplace_back(vstart);
        swap(neibors_1a2hop[0], neibors_1a2hop[neibors_1a2hop.size() - 1]);
        //In fact, neibors_1a2hop contains the vstart at the first

        //restore the global variables
    CLEAN:
        for (auto v : neibors_1a2hop ){
            visited[v] = 0;
        }
        for (auto v : neibors_1a2hopLeft){
            visited[v] = 0;
        }
        for (int v : neibors) {
            //the number of common neighbors with vstart >= lb-2*k+2
            if (count[v] >= lb_2k + 2) {
                neibors_1a2hop.emplace_back(v);
            }
            count[v] = 0;
        }
        for (int v : neiborsLeft) {
            //the number of common neighbors with vstart >= lb+1-2*k+2
            if (count[v] >= lb_2k+3) {
                neibors_1a2hopLeft.emplace_back(v);
            }
            count[v] = 0;
        }
        neibors.clear();
        neiborsLeft.clear();
        return sz+1;
    }

    //Build edges among the local nodes 
    Graph buildSubgraph(const Graph &g,const vector<int>& blk,const vector<int>& left,const vector<uint8_t> &inBlk, const vector<uint8_t> &inLeft,const int hopSz,vector<uint8_t> &isNei,vector<uint8_t> &commonMtx){  
        const int blkSz = blk.size();
        const int leftSz = left.size();
        vertex<int> *vertices = newA(vertex<int>, blkSz);
        int edgesHop=0,edges=0;
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
            if(i<hopSz)edgesHop+=vertices[i].degreeHop;
            edges+=vertices[i].degree;
        }
        const bool proper=(double(edgesHop)/(hopSz*(hopSz-1)))<0.9 && (double(edges)/(blkSz*(blkSz-1)))<0.85;
        return Graph(vertices, blkSz, 0, nullptr, proper);
    }

    Graph peelGraph(const Graph &g,bool* const mark,int * const resNei){
        const int n=g.n;
        #pragma omp parallel
        {
            thread_local queue<int> Q;
            #pragma omp for
            for(int i=0;i<n;++i){
                if(g.V[i].degree<bd){
                    mark[i]=false;
                    Q.push(i);
                }
                else{
                    mark[i]=true;
                    resNei[i]=g.V[i].degree;
                }
            }
            while(Q.size()){
                const int ele=Q.front();
                Q.pop();
                for(int i=0;i<g.V[ele].degree;++i){
                    const int nei=g.V[ele].Neighbors[i];
                    if(mark[nei]){
                        int old=resNei[nei];
                        while(!utils::CAS(&resNei[nei],old,old-1)){
                            old=resNei[nei];
                        }
                        if(old==bd){
                            mark[nei]=false;
                            Q.push(nei);
                        }
                    }
                }
            }
        }

        int* const map=new int[n];
        #pragma omp parallel for
        for(int i=0;i<n;++i)map[i]=i;
        _seq<int> leadList = sequence::pack(map, mark, n);
        const int pn=leadList.n;
        #pragma omp parallel for
        for(int i=0;i<pn;++i)map[leadList.A[i]]=i;
        vertex<int> *vertices = newA(vertex<int>, pn);
        #pragma omp parallel for
        for (int i = 0; i < pn; i++){
            const int ori = leadList.A[i];
            const int en=resNei[ori];
            vertices[i].Neighbors = newA(int,en);
            vertices[i].NeighborsLeft = nullptr;
            vertices[i].degree=en;
            int cursor=0;
            for(int j=0;j<g.V[ori].degree;++j){
                const int nei=g.V[ori].Neighbors[j];
                if(mark[nei])vertices[i].Neighbors[cursor++]=map[nei];
            }
        }
        delete[] map;
        leadList.del();
        return Graph(vertices, pn, 0, nullptr);
    }

    void buildCommonMtx(const Graph &subg,const int blkSz,const int hopSz,vector<uint8_t> &commonMtx){
        const int thresPP1=lb-2*k-max(k-3,0),thresPP2=lb-2*k+2-max(k-3,0);
        const int thresPC1=lb-2*k-max(k-2,0),thresPC2=lb-2*k+2-max(k-2,0);
        const int thresCC1=lb-2*k-(k-1),thresCC2=lb-2*k+2-(k-1);
        for(int i=0;i<hopSz;++i){
            for(int j=0;j<i;++j){
                const int common=utils::commonEle(
                    subg.V[i].Neighbors,
                    subg.V[i].degreeHop,
                    subg.V[j].Neighbors,
                    subg.V[j].degreeHop
                );
                if(commonMtx[i*blkSz+j]){
                    if(common>thresCC1)commonMtx[i*blkSz+j]=LINK2MORE;
                    else if(common==thresCC1)commonMtx[i*blkSz+j]=LINK2EQUAL;
                    else commonMtx[i*blkSz+j]=LINK2LESS;
                }
                else{
                    if(common>thresCC2)commonMtx[i*blkSz+j]=UNLINK2MORE;
                    else if(common==thresCC2)commonMtx[i*blkSz+j]=UNLINK2EQUAL;
                    else commonMtx[i*blkSz+j]=UNLINK2LESS;
                }
                commonMtx[j*blkSz+i]=commonMtx[i*blkSz+j];
            }
        }
        for(int i=hopSz;i<blkSz;++i){
            for(int j=0;j<hopSz;++j){
                const int common=utils::commonEle(
                    subg.V[i].Neighbors,
                    subg.V[i].degreeHop,
                    subg.V[j].Neighbors,
                    subg.V[j].degreeHop
                );
                if(commonMtx[i*blkSz+j]){
                    if(common>thresPC1)commonMtx[i*blkSz+j]=LINK2MORE;
                    else if(common==thresPC1)commonMtx[i*blkSz+j]=LINK2EQUAL;
                    else commonMtx[i*blkSz+j]=LINK2LESS;
                }
                else{
                    if(common>thresPC2)commonMtx[i*blkSz+j]=UNLINK2MORE;
                    else if(common==thresPC2)commonMtx[i*blkSz+j]=UNLINK2EQUAL;
                    else commonMtx[i*blkSz+j]=UNLINK2LESS;
                }
                commonMtx[j*blkSz+i]=commonMtx[i*blkSz+j];
            }
            if(k==2)continue;
            for(int j=hopSz;j<i;++j){
                const int common=utils::commonEle(
                    subg.V[i].Neighbors,
                    subg.V[i].degreeHop,
                    subg.V[j].Neighbors,
                    subg.V[j].degreeHop
                );
                if(commonMtx[i*blkSz+j]){
                    if(common>thresPP1)commonMtx[i*blkSz+j]=LINK2MORE;
                    else if(common==thresPP1)commonMtx[i*blkSz+j]=LINK2EQUAL;
                    else commonMtx[i*blkSz+j]=LINK2LESS;
                }
                else{
                    if(common>thresPP2)commonMtx[i*blkSz+j]=UNLINK2MORE;
                    else if(common==thresPP2)commonMtx[i*blkSz+j]=UNLINK2EQUAL;
                    else commonMtx[i*blkSz+j]=UNLINK2LESS;
                }
                commonMtx[j*blkSz+i]=commonMtx[i*blkSz+j];
            }
        }
    }

    //Local decompose --> Tree Search
    void decomposableSearch(const Graph &g){
        int *dpos = new int[g.n];
        int *dseq = new int[g.n];
        bool *mark = new bool[g.n];
        int *resNei = new int[g.n];
        #pragma omp parallel for
        for(int i=0;i<g.n;++i)dpos[i]=INT_MAX;
        int validblk=0;
        long cntMaxPlex=0;

        Timer watch;
        watch.start();
        
        Graph peelG=peelGraph(g,mark,resNei);
        const int pn=peelG.n;
        volatile bool* const ready=(volatile bool*)mark;
        #pragma omp parallel
        {
            thread_local vector<uint8_t> visited(pn,0);
            thread_local vector<int> count(pn,0);
            thread_local vector<int> neibors;
            thread_local vector<int> neiborsLeft;
            thread_local Stack<int> plex(PSIZE_LIMIT);
            thread_local VtxSet cand1(CAP_LIMIT1);
            thread_local VtxSet cand2(CAP_LIMIT2);
            thread_local VtxSet excl(CAP_LIMIT2);
            thread_local VtxSet exclBK(CAP_LIMIT2);
            thread_local Stack<int> exclStack(CAP_LIMIT2);
            thread_local vector<int> neiInP(CAP_LIMIT2,0);
            thread_local vector<int> neiInG(CAP_LIMIT2);
            thread_local vector<uint8_t> commonMtx(CAP_LIMIT2*CAP_LIMIT2);
            thread_local vector<int> blk;//vstart and its right 1a2hop-neighbors
            thread_local vector<int> left;//vstart's left 1a2hop-neighbors
            thread_local vector<uint8_t> inBlk(pn,0);
            thread_local vector<uint8_t> inLeft(pn,0);

            #pragma omp for
            for(int i=0;i<pn;++i)mark[i]=false;
            #pragma omp master
            {
                ListLinearHeap<int> *linear_heap = new ListLinearHeap<int>(pn, pn-1);
                linear_heap->init(pn, pn-1);
                for(int i=0;i<pn;++i){
                    linear_heap->insert(i,peelG.V[i].degree);
                }
                for (int i = 0; i < pn; i++) {
                    int u, key;
                    linear_heap->pop_min(u, key);
                    dpos[u] = i;
                    dseq[i] = u;
                    ready[i] = true;
                    for (int j = 0; j < peelG.V[u].degree; j++){
                        const int nei = peelG.V[u].Neighbors[j];
                        if(dpos[nei]==INT_MAX) {
                            linear_heap->decrement(nei);
                        }
                    }
                }
                delete linear_heap;
            }

            #pragma omp for schedule(dynamic) reduction(+:validblk)
            for(int idx = 0; idx < pn; idx++){
                while(!ready[idx]);
                // Search with leading vertex vstart;
                const int vstart=dseq[idx];

                // decompose the graph with leading vertex vstart
                blk.clear();
                left.clear();
                const int hopSz=decompose(peelG,vstart,blk,left,dpos,visited,count,neibors,neiborsLeft);
                const int blkSz=blk.size();
                if (blkSz>=lb) {
                    const int hop2Sz=blkSz-hopSz;
                    const int leftSz=left.size();
                    validblk++;
                    //updata task statistics
                    #pragma omp atomic
                    taskCnt++;
                    //build subgraph 
                    for(const auto u: blk)
                        inBlk[u] = true;
                    for(const auto u: left)
                        inLeft[u] = true;
                    Graph subg = buildSubgraph(peelG,blk,left,inBlk,inLeft,hopSz,visited,commonMtx);
                    if(subg.proper)buildCommonMtx(subg,blkSz,hopSz,commonMtx);

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

                    DecomposeMaxChecker maxchecker(subg);
                    PlexEmitor emitor;
                    KplexListor<DecomposeMaxChecker, PlexEmitor> listor(
                        subg, &maxchecker, &emitor, &plex, &cand1, &cand2, &excl, &exclBK, &exclStack,
                        neiInP.data(), neiInG.data(), hopSz, hop2Sz, commonMtx);

                    #pragma omp taskgroup
                    listor.cand1ToPlex(0)->kSearch(k-1);//start from the first vertex
                    
                    for(const auto u: blk)
                        inBlk[u] = false;
                    for(const auto u: left)
                        inLeft[u] = false;
                    listor.plexToCand1();
                    cand1.clear();
                    cand2.clear();
                    subg.del();
                    //updata task statistics
                    #pragma omp atomic
                    overCnt++;
                }
            }
        }
        #pragma omp parallel reduction(+:cntMaxPlex)
        cntMaxPlex+=cntT;

        double durat = watch.stop();
        printf("Total valid block %d\n", validblk);
        printf("Number of %d-plex: %ld\n", k, cntMaxPlex);
        printf("Running time %.3f\n",durat);
        printf("Task Cnt: %d\n", taskCnt);
        
        delete[] dseq;
        delete[] resNei;
        delete[] dpos;
        delete[] mark;
        peelG.del();
    }
}
int main(int argc, char **argv)
{
    if (argc == 4 || argc == 5){       
        k = atoi(argv[2]);
        lb = atoi(argv[3]);
        bd=lb-k;
        if(argc==5)
            setWorkers(workers = atoi(argv[4]));
        else 
            workers = getWorkers();
        workers1t2 = workers >> 1;
        graph<intT> g = ListPlex::readBinaryGraph(argv[1]); 
        ListPlex::decomposableSearch(g);
        g.del();
    }else {
        fprintf(stderr, "usage: listPlex <filename> <k> <lb> [worker-num]\n");
    }
    return 0;
}