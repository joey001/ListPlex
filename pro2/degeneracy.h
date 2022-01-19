#pragma once

#include <vector>
#include "graph.h"

template<typename intT>
class ListLinearHeap {
private:
	intT n; // number vertices
	intT key_cap; // the maximum allowed key value

	intT max_key; // possible max key
	intT min_key; // possible min key

	intT *keys; // keys of vertices
			  // keys[i] > key_cap if vertex i is not in the data structure

	intT *heads; // head of doubly-linked list for a specific weight
	intT *pres; // pre for doubly-linked list
	intT *nexts; // next for doubly-linked list

public:
	ListLinearHeap(intT _n, intT _key_cap) {
		n = _n;
		key_cap = _key_cap;

		min_key = key_cap;
		max_key = 0;

		heads = keys = pres = nexts = nullptr;
	}
	~ListLinearHeap() {
		if(heads != nullptr) {
			delete[] heads;
			heads = nullptr;
		}
		if(pres != nullptr) {
			delete[] pres;
			pres = nullptr;
		}
		if(nexts != nullptr) {
			delete[] nexts;
			nexts = nullptr;
		}
		if(keys != nullptr) {
			delete[] keys;
			keys = nullptr;
		}
	}

	// initialize the data structure by (id, key) pairs
	// _n is the number of pairs, _key_cap is the maximum possible key value
	void init(intT _n, intT _key_cap) {
		keys = new intT[n];
		pres = new intT[n];
		nexts = new intT[n];
		heads = new intT[key_cap+1];
		min_key = _key_cap; max_key = 0;
		#pragma omp parallel for
		for(intT i = 0;i <= _key_cap;i ++) heads[i] = n;
	}

	// insert (id, key) pair into the data structure
	void insert(intT id, intT key) {
		assert(id < n); assert(key <= key_cap);
		//assert(keys[id] > key_cap);

		keys[id] = key; pres[id] = n; nexts[id] = heads[key];
		if(heads[key] != n) pres[heads[key]] = id;
		heads[key] = id;

		if(key < min_key) min_key = key;
		if(key > max_key) max_key = key;
	}

	// remove a vertex from the data structure
	intT remove(intT id) {
		assert(keys[id] <= max_key);
		if(pres[id] == n) {
			assert(heads[keys[id]] == id);
			heads[keys[id]] = nexts[id];
			if(nexts[id] != n) pres[nexts[id]] = n;
		}
		else {
			intT pid = pres[id];
			nexts[pid] = nexts[id];
			if(nexts[id] != n) pres[nexts[id]] = pid;
		}

		return keys[id];
	}

	intT get_n() { return n; }
	intT get_key_cap() { return key_cap; }
	intT get_key(intT id) { return keys[id]; }

	void get_ids(std::vector<intT> &ids) {
		ids.clear();
		tighten();
		for(intT i = min_key;i <= max_key;i ++) {
			for(intT id = heads[i];id != n;id = nexts[id]) {
				ids.pb(id);
			}
		}
	}

	void get_ids_keys(std::vector<intT> &ids, std::vector<intT> &_keys) {
		ids.clear(); _keys.clear();
		tighten();
		for(intT i = min_key;i <= max_key;i ++) {
			for(intT id = heads[i];id != n;id = nexts[id]) {
				ids.pb(id); _keys.pb(id);
			}
		}
	}

	bool empty() {
		tighten();
		return min_key > max_key;
	}

	intT size() {
		tighten();
		intT res = 0;
		for(intT i = min_key;i <= max_key;i ++) for(intT id = heads[i];id != n;id = nexts[id]) ++ res;
		return res;
	}

	// get the (id,key) pair with the maximum key value; return true if success, return false otherwise
	bool get_max(intT &id, intT &key) {
		if(empty()) return false;

		id = heads[max_key];
		key = max_key;
		assert(keys[id] == key);
		return true;
	}

	// pop the (id,key) pair with the maximum key value; return true if success, return false otherwise
	bool pop_max(intT &id, intT &key) {
		if(empty()) return false;

		id = heads[max_key];
		key = max_key;
		assert(keys[id] == key);

		heads[max_key] = nexts[id];
		if(heads[max_key] != n) pres[heads[max_key]] = n;
		return true;
	}

	// get the (id,key) pair with the minimum key value; return true if success, return false otherwise
	bool get_min(intT &id, intT &key) {
		if(empty()) return false;

		id = heads[min_key];
		key = min_key;
		assert(keys[id] == key);

		return true;
	}

	// pop the (id,key) pair with the minimum key value; return true if success, return false otherwise
	void pop_min(intT &id, intT &key) {
		tighten();

		id = heads[min_key];
		key = min_key;
		assert(keys[id] == key);

		heads[min_key] = nexts[id];
		if(heads[min_key] != n) pres[heads[min_key]] = n;
	}

	// increment the key of vertex id by inc
	intT increment(intT id, intT inc = 1) {
		assert(keys[id]+inc <= key_cap);

		intT new_key = keys[id] + inc;

		remove(id);
		insert(id, new_key);

		return new_key;
	}

	// decrement the key of vertex id by dec
	intT decrement(const intT id,const intT dec = 1) {
		assert(keys[id] >= dec);

		const intT new_key = keys[id] - dec;

		remove(id);
		insert(id, new_key);

		return new_key;
	}

private:
	void tighten() {
		while(heads[min_key] == n) ++ min_key;
		// while(min_key <= max_key&&heads[max_key] == n) -- max_key;
	}
};

void degeneracyOrder(const graph<int> &g, intT *dseq, int *resNei, int *dpos) {
	const int n=g.n;
	ListLinearHeap<int> *linear_heap = new ListLinearHeap<int>(n, n-1);
	linear_heap->init(n, n-1);
	for(int i=0;i<n;++i){
		linear_heap->insert(i,g.V[i].degree);
	}
	for (int i = 0; i < n; i++) {
		int u, key;
		linear_heap->pop_min(u, key);
		dpos[u] = i;
		dseq[i] = u;
		resNei[u]=key;
		for (int j = 0; j < g.V[u].degree; j++){
			const int nei = g.V[u].Neighbors[j];
			if(dpos[nei]==INT_MAX) {
				linear_heap->decrement(nei);
			}
		}
	}
	delete linear_heap;
}
