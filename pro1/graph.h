#ifndef _GRAPH_INCLUDED
#define _GRAPH_INCLUDED

#include <iostream>
#include <algorithm>
#include "utils.h"
#include "config.h"

template <class intT>
struct vertex
{
  intT *Neighbors;
  intT *NeighborsLeft;
  intT degree,degreeHop,degreeLeft;
  void del() {
    free(Neighbors);
    free(NeighborsLeft);
  }
};

template <class intT>
struct graph
{
  vertex<intT> *V;
  const intT n;
  intT m;   
  intT *allocatedInplace;
  graph():n(0),m(0),V(nullptr),allocatedInplace(nullptr){};
  graph(vertex<intT> *VV, intT nn, uintT mm)
      : V(VV), n(nn), m(mm), allocatedInplace(NULL) {}
  graph(vertex<intT> *VV, intT nn, uintT mm, intT *ai)
      : V(VV), n(nn), m(mm), allocatedInplace(ai) {}
  graph copy()
  {
    vertex<intT> *VN = newA(vertex<intT>, n);
    intT *_allocatedInplace = newA(intT, n + m + 2);
    _allocatedInplace[0] = n;
    _allocatedInplace[1] = m;
    intT *Edges = _allocatedInplace + n + 2;
    intT k = 0;
    for (intT i = 0; i < n; i++) //for each vertex
    {
      _allocatedInplace[i + 2] = allocatedInplace[i + 2]; //copy 
      VN[i] = V[i];
      VN[i].Neighbors = Edges + k;
      for (intT j = 0; j < V[i].degree; j++)
        Edges[k++] = V[i].Neighbors[j];
    }
    return graph(VN, n, m, _allocatedInplace);
  }
  void del()
  {
    if (allocatedInplace == nullptr)
      for (intT i = 0; i < n; i++)
        V[i].del();
    else
      free(allocatedInplace);
    free(V);
  }
  bool isAdj(intT v1, intT v2) const
  {
      if(V[v1].degree>V[v2].degree)
          std::swap(v1, v2);
      return std::binary_search(V[v1].Neighbors, V[v1].Neighbors + V[v1].degree, v2);
  }
  bool isAdj11(intT v1, intT v2) const
  {
      if(V[v1].degreeHop>V[v2].degreeHop)
          std::swap(v1, v2);
      return std::binary_search(V[v1].Neighbors, V[v1].Neighbors + V[v1].degreeHop, v2);
  }
};

#endif