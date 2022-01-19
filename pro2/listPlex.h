#include <cassert>
#include <functional>
#include "utils.h"
#include "timer.h"
#include "MRandSet.h"
#include "MStack.h"
#include "degeneracy.h"
#include "graph.h"
#include "graphIO.h"
#include "sequence.h"
namespace ListPlex{
    using VtxSet = RandSet<int>;
    using Graph = graph<int>;
    enum : uint8_t{
        //odd for linked
        UNLINK2LESS=0,
        LINK2LESS=1,
        UNLINK2EQUAL=2,
        LINK2EQUAL=3,
        UNLINK2MORE=4,
        LINK2MORE=5
    };

    struct PlexEmitor;
    struct DecomposeMaxChecker;

    template <typename MaximalityChecker, typename Emitor>
    struct KplexListor
    {
        const Graph subg;
        MaximalityChecker *plexMaxChecker;
        Emitor *emitor;
        
        int * neiInP;
        int * const neiInG;
        Stack<int> * const plex;
        VtxSet * const cand1;
        VtxSet * const cand2;
        VtxSet *excl;
        VtxSet * const exclBK;
        Stack<int> * const exclStack;
        const int hopSz, hop2Sz;
        const vector<uint8_t>& commonMtx;

        KplexListor(Graph &_subg,MaximalityChecker *_plexMaxChecker, Emitor *_emitor,
                    Stack<int> *_plex, VtxSet *_cand1, VtxSet *_cand2, VtxSet *_excl, VtxSet *_exclBK, Stack<int> *_exclStack,
                    int *_neiInP, int *_neiInG, int _hopSz, int _hop2Sz, const vector<uint8_t>& _commonMtx)
        :subg(_subg),plexMaxChecker(_plexMaxChecker),emitor(_emitor)
        ,plex(_plex),cand1(_cand1),cand2(_cand2),excl(_excl),exclBK(_exclBK),exclStack(_exclStack)
        ,neiInP(_neiInP),neiInG(_neiInG),hopSz(_hopSz),hop2Sz(_hop2Sz),commonMtx(_commonMtx)
        {
        }

        KplexListor(const KplexListor& o)
        :subg(o.subg),plexMaxChecker(o.plexMaxChecker),emitor(o.emitor)
        ,plex(new Stack<int>(*(o.plex))),cand1(new VtxSet(*(o.cand1))),cand2(nullptr),excl(new VtxSet(*(o.excl))),exclBK(nullptr),exclStack(new Stack<int>(o.exclStack->cap))
        ,neiInP(new int[subg.n]),neiInG(new int[subg.n]),hopSz(o.hopSz),hop2Sz(o.hop2Sz),commonMtx(o.commonMtx)
        {
            memcpy(neiInP, o.neiInP,sizeof(int)*subg.n);
            memcpy(neiInG, o.neiInG, sizeof(int)*subg.n);
        }

        void del(){
            delete plex;
            delete cand1;
            delete excl;
            delete exclStack;
            delete[] neiInP;
            delete[] neiInG;
        }

        int getIdx(const int v1,const int v2){
            return v1*subg.n+v2;
        }

        bool isAdjMtx(const int v1,const int v2)
        {
            return commonMtx[getIdx(v1,v2)]&1;
        }

        void addG(const int u){
            for (int i = 0; i < subg.V[u].degree; i++)
            {
                const int nei = subg.V[u].Neighbors[i];
                neiInG[nei]++;
            }
        }
        void subG(const int u){
            for (int i = 0; i < subg.V[u].degree; i++)
            {
                const int nei = subg.V[u].Neighbors[i];
                neiInG[nei]--;
            }
        }
        void addP(const int u){
            for (int i = 0; i < subg.V[u].degree; i++)
            {
                const int nei = subg.V[u].Neighbors[i];
                neiInP[nei]++;
            }
        }        
        void subP(const int u){
            for (int i = 0; i < subg.V[u].degree; i++)
            {
                const int nei = subg.V[u].Neighbors[i];
                neiInP[nei]--;
            }
        }

        KplexListor* cand1ToPlex(const int v)
        {
            plex->push(v);
            cand1->remove(v);
            addP(v);
            return this;
        }

        KplexListor* plexToCand1()
        {
            assert(plex->sz > 0);
            const int u = plex->top();
            cand1->add(u);
            plex->pop();
            subP(u);
            return this;
        }

        int cand2BackToPlex()
        {
            const int v=cand2->pop_back();
            plex->push(v);
            addP(v);
            return v;
        }

        int cand2BackToPlexInc()
        {
            const int v=cand2->pop_back();
            plex->push(v);
            addP(v);
            addG(v);
            return v;
        }

        KplexListor* plexToCand2()
        {
            assert(plex->sz > 0);
            const int u = plex->top();
            cand2->add(u);
            plex->pop();
            subP(u);
            return this;
        }

        KplexListor* plexToCand2Inc()
        {
            assert(plex->sz > 0);
            const int u = plex->top();
            cand2->add(u);
            plex->pop();
            subP(u);
            subG(u);
            return this;
        }

        KplexListor* cand1ToExcl(const int v)
        {
            cand1->remove(v);
            excl->add(v);
            subG(v);
            return this;
        }

        KplexListor* exclToCand1(const int v)
        {
            excl->remove(v);
            cand1->add(v);
            addG(v);
            return this;
        }

        int cand2BackToExcl()
        {
            const int u=cand2->pop_back();
            excl->add(u);
            subG(u);
            return u;
        }

        int cand2BackToExclInc()
        {
            const int u=cand2->pop_back();
            excl->add(u);
            return u;
        }

        KplexListor* exclToCand2(const int v)
        {
            excl->remove(v);
            cand2->add(v);
            addG(v);
            return this;
        }

        KplexListor* exclToCand2Inc(const int v)
        {
            excl->remove(v);
            cand2->add(v);
            return this;
        }

        KplexListor* updateExcl(int& recExcl,const int v2add){
            recExcl=excl->sz;
            for (int i = 0; i < excl->sz;)
            {
                const int ele = excl->members[i];
                if((subg.proper&&UNLINK2MORE>commonMtx[getIdx(ele,v2add)])||!canFormPlex(ele,1)){
                    excl->remove(ele);
                    exclStack->push(ele);
                }
                else ++i;
            }
            recExcl-=excl->sz;
            return this;
        }
        
        KplexListor* updateExclK(int& recExcl,const int v2add){
            recExcl=excl->sz;
            for (int i = 0; i < excl->sz;)
            {
                int ele = excl->members[i];
                if(subg.proper&&UNLINK2MORE>commonMtx[getIdx(ele,v2add)]){
                    excl->remove(ele);
                    exclStack->push(ele);
                }
                else ++i;
            }
            recExcl-=excl->sz;
            return this;
        }

        KplexListor* updateCand1(int& recCand1,const int v2add){
            recCand1=cand1->sz;
            for (int i = 0; i < cand1->sz;)
            {
                const int ele = cand1->members[i];
                if((subg.proper&&UNLINK2EQUAL>commonMtx[getIdx(ele,v2add)])||!canFormPlex(ele,0)){
                    cand1->remove(ele);
                    subG(ele);
                }
                else ++i;
            }
            recCand1-=cand1->sz;
            return this;
        }

        KplexListor* updateCand1Fake(int& recCand1,const  int v2add){
            recCand1=cand1->sz;
            for (int i = 0; i < cand1->sz;)
            {
                int ele = cand1->members[i];
                if((subg.proper&&UNLINK2EQUAL>commonMtx[getIdx(v2add,ele)])||!canFormPlex(ele,0)){
                    cand1->fakeRemove(ele);
                    subG(ele);
                }
                else ++i;
            }
            recCand1-=cand1->sz;
            return this;
        }

        KplexListor* updateCand1KFake(int& recCand1,const  int v2add){
            recCand1=cand1->sz;
            for (int i = 0; i < cand1->sz;)
            {
                int ele = cand1->members[i];
                if(subg.proper&&UNLINK2EQUAL>commonMtx[getIdx(v2add,ele)]){
                    cand1->fakeRemove(ele);
                    subG(ele);
                }
                else ++i;
            }
            recCand1-=cand1->sz;
            return this;
        }

        KplexListor* updateCand2Fake(int& recCand2,const int v2add){
            recCand2=cand2->sz;
            for (int i = 0; i < cand2->sz;)
            {
                int ele = cand2->members[i];
                if(subg.proper&&UNLINK2EQUAL>commonMtx[getIdx(v2add,ele)]){
                    cand2->fakeRemove(ele);
                    subG(ele);
                }
                else ++i;
            }
            recCand2-=cand2->sz;
            return this;
        }

        KplexListor* updateCand2FakeInc(int& recCand2,const int v2add){
            recCand2=cand2->sz;
            for (int i = 0; i < cand2->sz;)
            {
                int ele = cand2->members[i];
                if(subg.proper&&UNLINK2EQUAL>commonMtx[getIdx(v2add,ele)]){
                    cand2->fakeRemove(ele);
                }
                else ++i;
            }
            recCand2-=cand2->sz;
            return this;
        }

        //Check if plex+u is a k-plex
        bool canFormPlex(const int u,const int extra)
        {
            if (neiInP[u] + k < plex->sz + 1 || neiInG[u] + k < max(lb,plex->sz)+extra)
                return false;
            for (int i = 0; i < plex->sz; i++)
            {
                const int v = plex->members[i];
                assert(v != u);
                if (neiInP[v] + k == plex->sz && !isAdjMtx(v, u))
                { //v is saturated by v,u are not neighbors
                    return false;
                }
            }
            return true;
        }

        void kSearch(const int res){
            int recExcl=0,recExclTmp;
            int recCand1[K_LIMIT],recCand2[K_LIMIT];
            if(cand2->sz==0){
                listByCase();
                return;
            }
            //br0;
            int v2delete=cand2BackToExclInc();
            kSearch(res);
            exclToCand2Inc(v2delete);
            
            //start from the 1st br.
            int br=1;
            for (; br<res; br++)
            {
                const int v2add=cand2BackToPlexInc();
                updateCand1KFake(recCand1[br],v2add);
                updateCand2FakeInc(recCand2[br],v2add);
                updateExclK(recExclTmp,v2add);
                recExcl+=recExclTmp;
                if(cand2->sz){
                    v2delete = cand2BackToExclInc();
                    kSearch(res-br);
                    exclToCand2Inc(v2delete);
                }
                else{
                    kSearch(res-br);
                    break;
                }
            }
            //The last branch.
            if(br==res){
                recCand2[br]=0;
                const int v2add=cand2BackToPlexInc();
                VtxSet* tmp=excl;
                excl=exclBK;
                updateCand1Fake(recCand1[br],v2add);
                listByCase();
                excl=tmp;
            }
            for (int i = br; i >=1 ; i--)
            {
                cand1->fakeRecoverAdd(recCand1[i],this);
                cand2->fakeRecover(recCand2[i]);
                plexToCand2Inc();
            }
            for(int i=0;i<recExcl;++i){
                excl->add(exclStack->top());
                exclStack->pop();
            }
        }

        //Check if plex+cand is a k-plex
        void checkAsMaxPlex()
        {
            memcpy(plex->members + plex->sz, cand1->members, cand1->sz * sizeof(int));
            plex->sz += cand1->sz;
            int *tmp = neiInP;
            neiInP = neiInG;
            //assume plex+cand as a plex

            bool flag = true;
            for (int i = 0; i < excl->sz;++i)
            {
                if (canFormPlex(excl->members[i],1))
                {
                    flag = false;
                    break;
                }
            }
            if (flag&&plexMaxChecker->isMaximal(plex->members, plex->sz, neiInP))
            {
                emitor->emitPlex();
            }
            plex->sz -= cand1->sz;
            neiInP = tmp;
        }

        void branchInCandBase(const int pivot)
        {
            //In the first branch, select pivot
            int recCand1,recExcl;
            cand1ToPlex(pivot)->updateCand1Fake(recCand1,pivot)->updateExcl(recExcl,pivot)->listByCase();

            //Recover
            cand1->fakeRecoverAdd(recCand1,this);
            plexToCand1();
            while(recExcl){
                excl->add(exclStack->top());
                exclStack->pop();
                recExcl--;
            }
            //In the second branch, remove pivot
            cand1ToExcl(pivot)->listByCase();
            exclToCand1(pivot);
        }

        void branchInCand(const int pivot)
        {
            //In the first branch, select pivot
            //update task statistics
            #pragma omp atomic
            taskCnt++;
            KplexListor *listor = new KplexListor(*this);
            #pragma omp task firstprivate(listor,pivot)
            {
                int tmp;
                listor->cand1ToPlex(pivot)->updateCand1Fake(tmp,pivot)->updateExcl(tmp,pivot)->listByCase();
                listor->del();
                delete listor;
                //update task statistics
                #pragma omp atomic
                overCnt++;
            }
            //In the second branch, remove pivot
            cand1ToExcl(pivot)->listByCase();
            exclToCand1(pivot);
        }

        void branchInPlexBase(const int pivot){
            int tot = k + neiInP[pivot] - plex->sz;
            assert(tot > 0);
            int addList[K_LIMIT];
            int recCand1[K_LIMIT]; 
            int recExcl=0,recExclTmp;
            int sz = 0;
            //Filter cand1 and find non-adjacent nodes
            for (VtxSet::TIndex i = 0; i < cand1->sz; i++)
            {
                const int u = cand1->members[i];
                if (!isAdjMtx(pivot, u))
                {
                    addList[sz++] = u;
                }
                if (sz == tot)
                    break;
            }
            assert(sz == tot);
            //There are at most tol+1 branches.
            //br0;
            int v2delete = addList[0];
            cand1ToExcl(v2delete)->listByCase();
            exclToCand1(v2delete);

            //Proceed br1 in [1,tol-1); there are at least 2 branches
            //add in addList[0,...,br), remove addList[br].
            int v2add;
            int movals = 0; // recording the number of vertices pushing into kplex

            //start from the 1st br.
            int br = 1;
            for (; br < tot; br++)
            {
                v2add = addList[br - 1];
                v2delete = addList[br];
                //if v2add has not been reduced by feasibility condition, i.e., v2add in cand1
                //or the graph size is above lb. stop
                if (cand1->contains(v2add) && plex->sz + cand1->sz >= lb) 
                {
                    //move v2 add to plex and update cand1
                    cand1ToPlex(v2add)->updateCand1(recCand1[movals],v2add)->updateExcl(recExclTmp,v2add);
                    recExcl+=recExclTmp;
                    movals +=1; //mark the number of pushed vertices
                    
                    if (cand1->contains(v2delete))
                    {
                        //if v2delete is still in cand1, remove it ;
                        cand1ToExcl(v2delete);
                        listByCase();
                        exclToCand1(v2delete);
                        //recover the remove vertices
                    }
                    else
                    {
                        listByCase();
                    }
                }
                else
                {
                    break;
                }
            }
            //The last (tolth) branch.
            if (br == tot && cand1->contains(addList[br - 1]) && plex->sz + cand1->sz >= lb)
            {
                v2add = addList[tot - 1];
                cand1ToPlex(v2add)->updateCand1(recCand1[movals],v2add)->updateExcl(recExclTmp,v2add); //all vertices in addList has moved to plex        
                recExcl+=recExclTmp;
                movals += 1;
                listByCase();
            }
            //recover cand1 and excl
            //pop the last br vertices
            movals--;
            while(movals>=0){
                cand1->recoverAdd(recCand1[movals],this);
                plexToCand1();
                movals--;
            }
            while(recExcl){
                excl->add(exclStack->top());
                exclStack->pop();
                recExcl--;
            }
        }

        void branchInPlex(const int pivot){
            const int tot = k + neiInP[pivot] - plex->sz;
            assert(tot > 0);
            int addList[K_LIMIT];
            int tmp = 0;
            //Filter cand1 and find non-adjacent nodes
            for (VtxSet::TIndex i = 0; i < cand1->sz; i++)
            {
                const int u = cand1->members[i];
                if (!isAdjMtx(pivot, u))
                {
                    addList[tmp++] = u;
                }
                if (tmp == tot)
                    break;
            }
            assert(tmp == tot);
            //There are at most tol+1 branches.
            //br0
            KplexListor* listorBK = new KplexListor(*this); //stand
            KplexListor* listorNOW= new  KplexListor(*this);
            KplexListor *listorNEXT;
            for (int j = 0; j < tot ;++j){
                #pragma omp atomic
                taskCnt++;
                //1.exclusive the first
                listorNOW->cand1ToExcl(addList[j]);
                if(j!=tot-1){
                    //2. store for the next
                    listorNEXT = new KplexListor(*listorNOW);
                    //3. mark the second
                    const int v2add=addList[j + 1];
                    listorNOW->cand1ToPlex(v2add)->updateCand1(tmp,v2add)->updateExcl(tmp,v2add);
                }
                //4. parallel
                #pragma omp task firstprivate(listorNOW)
                {
                    listorNOW->listByCase();
                    listorNOW->del();
                    delete listorNOW;
                    #pragma omp atomic
                    overCnt++;
                }
                listorNOW = listorNEXT;
            }

            //Proceed br1 in [1,tol-1); there are at least 2 branches
            //add in addList[0,...,br), remove addList[br].
            int v2delete,v2add;
            //start from the 1st br.
            int br = 1;
            for (; br < tot; br++)
            {
                v2add = addList[br - 1];
                v2delete = addList[br];
                //if v2add has not been reduced by feasibility condition, i.e., v2add in cand1
                //or the graph size is above lb. stop
                if (listorBK->cand1->contains(v2add) && listorBK->plex->sz + listorBK->cand1->sz >= lb) 
                {
                    #pragma omp atomic
                    taskCnt++;
                    //move v2 add to plex and update cand1
                    listorBK->cand1ToPlex(v2add)->updateCand1(tmp,v2add)->updateExcl(tmp,v2add);
                    KplexListor* listorBR = new KplexListor(*listorBK);
                    #pragma omp task firstprivate(listorBR)
                    {
                        //if v2delete is still in cand1, remove it ;
                        if (listorBR->cand1->contains(v2delete))listorBR->cand1ToExcl(v2delete);
                        listorBR->listByCase();
                        listorBR->del();
                        delete listorBR;
                        #pragma omp atomic
                        overCnt++;
                    }
                }
                else
                {
                    break;
                }
            }
            //The last (tolth) branch.
            if (br == tot && listorBK->cand1->contains(addList[br - 1]) && listorBK->plex->sz + listorBK->cand1->sz >= lb)
            {
                v2add = addList[tot - 1];
                listorBK->cand1ToPlex(v2add)->updateCand1(tmp,v2add)->updateExcl(tmp,v2add); //all vertices in addList has moved to plex        
                listorBK->listByCase();
            }
            listorBK->del();
            delete listorBK;
        }
        

        void listByCase()
        {
            if (plex->sz + cand1->sz < lb){
                return;
            }
            if (cand1->sz == 0)
            {
                if (excl->sz == 0 && plexMaxChecker->isMaximal(plex->members, plex->sz, neiInP))
                {
                    emitor->emitPlex();
                }
                return;
            }
            int minnei=INT_MAX;int pivot;
            auto maxswap = [&](int u) {
                if (neiInG[u] < minnei)
                {
                    minnei = neiInG[u];
                    pivot = u;
                }
            };
            plex->for_each(maxswap);
            const int minneiPlex=minnei;
            if(minneiPlex + k < max(lb,plex->sz))
                return;
            cand1->for_each(maxswap);
            //TODO: for better branching prediction:change to switch
            
            if (minnei >= plex->sz + cand1->sz - k)
            {
                checkAsMaxPlex();   //P+C is a k-plex
                return;
            }
            const int grainSize = 10;
            if (minneiPlex==minnei) //Pivot is in plex
            {
                if(taskCnt>=workers1t2&&(taskCnt-overCnt<workers)&&(cand1->sz>grainSize))
                    branchInPlex(pivot);
                else branchInPlexBase(pivot);
            }
            else
            {
                if((taskCnt>=workers1t2)&&(taskCnt-overCnt<workers)&&(cand1->sz>grainSize))
                    branchInCand(pivot);
                else branchInCandBase(pivot);
            }
        }

        void printContextInfo(){
            printf("\n\n[DEPTH]: %d\n", plex->sz);
            printf("PLEX:");
            plex->for_each([](int v){printf("%d ",v);});
            printf("\nCAND1:");
            cand1->for_each([](int v){printf("%d ",v);});
            printf("\nCAND2:");
            cand2->for_each([](int v){printf("%d ",v);});
            printf("\nEXCLUSIVE:");
            excl->for_each([](int v){printf("%d ",v);});
            printf("\n");
            
            for (int i = 0; i < subg.n; i++)
                printf("[%d]%d:%d ", i, neiInP[i], neiInG[i]);
            printf("\n");
        }
    };
}
