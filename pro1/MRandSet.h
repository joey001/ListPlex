#pragma once
#include <cassert>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include "MBitSet.h"
#include "config.h"


template <typename T=int> 
class RandSet{    
public:
    using TIndex = uint16_t;
    int sz;
    int cap;
    T *members;
    TIndex *pos; 
    MBitSet32 *flag;
    
public:

    //default constructor        
    RandSet():sz(0),cap(0),members(nullptr),flag(nullptr),pos(nullptr){}

    RandSet(int _maxvalue):cap(_maxvalue),sz(0){       
        members = new T[_maxvalue];
        pos = new TIndex[_maxvalue];
        flag = new MBitSet32(_maxvalue);
    }

    
    RandSet(const RandSet &_set){// copy constructor
        cap = _set.cap;
        sz = _set.sz;
        
        members = new T[cap];
        memcpy(members, _set.members, sizeof(T) * sz);

        pos = new TIndex[cap];
        memcpy(pos, _set.pos, sizeof(TIndex) * cap);

        flag = new MBitSet32(*(_set.flag));        
        
    }
    RandSet& operator=(const RandSet &_set) = delete;

    ~RandSet(){
        delete[] members;
        delete[] pos;
        delete flag;        
    }

    void init(int _maxvalue){       
        cap=_maxvalue;
        sz=0;
        members = new T[_maxvalue];
        pos = new TIndex[_maxvalue];
        flag = new MBitSet32(_maxvalue);
    }

    void resize(const int _cap){
        cap=_cap;
        flag->cap=_cap;
        flag->n=_cap>>5;
    }

    bool contains(const T& v){
        assert(v<cap);        
        return flag->test(v);
    }

    void add(const T& v){
        pos[v] = sz; 
        members[sz++]=v;
        flag->set(v); 
    }

    void fakeRemove(const T& v){
        //swap
        TIndex idx = pos[v];            
        pos[members[sz-1]] = idx;
        pos[v] = sz-1;
        std::swap(members[idx],members[--sz]);
    }

    void remove(const T& v){
        fakeRemove(v);
        flag->set(v);
    }

    template<typename Type>
    void recoverAdd(const int len,Type* listor){
        T* cursor=members+sz;
        for(int i=0;i<len;++i){
            const int u=*cursor;
            flag->set(u);
            listor->addG(u);
            cursor++;
        }
        sz+=len;
    }

    template<typename Type>
    void fakeRecoverAdd(const int len,Type* listor){
        T* cursor=members+sz;
        for(int i=0;i<len;++i){
            const int ele=*cursor;
            listor->addG(ele);
            cursor++;
        }
        sz+=len;
    }

    void fakeRecover(const int len){
        sz+=len;
    }
     
    int pop_back(){
        int v=members[--sz];
        flag->set(v);
        return v;
    }

    RandSet& operator+(const T& v){
        add(v);
        return *this;
    }
    
    void clear(){
        flag->clear();
        sz = 0;                
    }

    void clearByEle(){
        for(int i=0;i<sz;++i){
            flag->set(members[i]);
        }
        sz = 0;                
    }

    template<typename Fnc>
    void for_each(const Fnc &f){        
        for (TIndex i = 0; i < sz; i++){
            f(members[i]);
        }
    }
};
