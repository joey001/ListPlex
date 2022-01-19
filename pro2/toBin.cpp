/*
This file transform text graph file to binary graph file.
Format used in SNAP and DIMACS10 are supported.
Format in SNAP are identified by the file which has suffix "txt".
Format in DIMACS10 are identified by the file which has suffix "graph".
All direct, weighted are transformed to undirected, unweighted.
Self-loops are eliminated in the transformed.
Singleton are eliminated in the transformed file. The vertex id are also re-assigned.
*/

#define _CRT_SECURE_NO_WARNINGS


#include <cassert>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <map>
#include <cstring>
//using namespace std;


using ui = unsigned int;
const ui FILELEN = 1024;

static int fileSuffixPos(char* filepath) {
	int j = strlen(filepath) - 1;
	while (filepath[j] != '.')
		j--;
	return j + 1;
}

static FILE *open_file(const char *file_name, const char *mode) {
	FILE *f = fopen(file_name, mode);
	if (f == nullptr) {
		printf("Can not open file: %s\n", file_name);
		exit(1);
	}

	return f;
}

static std::string integer_to_string(long long number) {
	std::vector<ui> sequence;
	if (number == 0) sequence.push_back(0);
	while (number > 0) {
		sequence.push_back(number % 1000);
		number /= 1000;
	}

	char buf[5];
	std::string res;
	for (unsigned int i = sequence.size(); i > 0; i--) {
		if (i == sequence.size()) sprintf(buf, "%u", sequence[i - 1]);
		else sprintf(buf, ",%03u", sequence[i - 1]);
		res += std::string(buf);
	}
	return res;
}
//----------------------------------------------------------------------
ui n, m;
ui* pstart;
ui* edges;
ui* reverse;
void readDIMACS2Text(const char* filepath) {
	std::ifstream infile;
	char buf[1024];
	char tmps1[1024], tmps2[1024];
	std::vector<std::pair<ui, ui> > epairs;
	std::vector<ui> nodes;
	infile.open(filepath, std::ios::in);
	if (!infile.is_open()) {
		fprintf(stderr, "can not find file %s\n", filepath);
		exit(1);
	}
	int from, to;
	infile.getline(buf, 1024);
	//BUG:buff can be potentially empty string
	while (buf[0] != 'p') infile.getline(buf, 1024);
	sscanf(buf, "%s %s %d %d", tmps1, tmps2, &n, &m);
	int lno = 1;
	while (infile.getline(buf, 1024)) {
		if (strlen(buf) == 0) continue;
		if (buf[0] != 'e') {
			fprintf(stderr, "ERROR in line %d\n", lno );
			continue;
		}
		sscanf(buf, "%s %d %d", tmps1, &from, &to);
		if (from != to) {
			epairs.push_back(std::make_pair(from, to));
			epairs.push_back(std::make_pair(to, from));
			nodes.push_back(from);
			nodes.push_back(to);
		}
		lno++;
	}
	sort(nodes.begin(), nodes.end());
	nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());

	sort(epairs.begin(), epairs.end());
	epairs.erase(unique(epairs.begin(), epairs.end()), epairs.end());

	ui contn = 1;
	std::map<ui, ui> idmp;
	for (ui i = 0; i < nodes.size(); i++) {
		idmp[nodes[i]] = i;
		if (nodes[i] != i) {
			contn = 0;
		}
	}
	if (contn == 0) printf("Node ids are not preserved! \n");

	n = nodes.size();
	m = epairs.size();
	printf("n = %s, (undirected) m = %s\n",
		integer_to_string(n).c_str(),
		integer_to_string(m / 2).c_str());

	pstart = new ui[n + 1];
	edges = new ui[m];
	reverse = new ui[m];
	ui j = 0;
	for (ui i = 0; i < n; i++) {
		pstart[i] = j;
		while (j < m && epairs[j].first == nodes[i]) {
			edges[j] = idmp[epairs[j].second];
			reverse[j] = i;
			++j;
		}
		std::sort(edges + pstart[i], edges + j);
	}
	pstart[n] = j;
}

void readRawSNAPText(const char* filepath) {
	std::ifstream infile;
	char buf[1024];
	std::vector<std::pair<ui, ui> > epairs;
	std::vector<ui> nodes;
	//FILE *f = Utility::open_file(filepath, "r");
	infile.open(filepath, std::ios::in);
	if (!infile.is_open()) {
		fprintf(stderr, "can not find file %s\n", filepath);
		exit(1);
	}
	int max_id = 0;
	int from, to;
	while (infile.getline(buf, 1024)) {
		char *p = buf;
		while (*p == ' ' && *p != '\0') p++;
		if (*p == '#' || *p == '\0') continue;
		std::stringstream ss(buf);
		ss >> from >> to;
		if (from != to) {
			epairs.push_back(std::make_pair(from, to));
			epairs.push_back(std::make_pair(to, from));
			nodes.push_back(from);
			nodes.push_back(to);
		}
	}
	infile.close();

	sort(nodes.begin(), nodes.end());
	nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());

	sort(epairs.begin(), epairs.end());
	epairs.erase(unique(epairs.begin(), epairs.end()), epairs.end());

	ui contn = 1;
	std::map<ui, ui> idmp;
	for (ui i = 0; i < nodes.size(); i++) {
		idmp[nodes[i]] = i;
		if (nodes[i] != i) {
			contn = 0;
		}
	}
	if (contn == 0) printf("Node ids are not preserved! \n");

	n = nodes.size();
	m = epairs.size();
	printf("n = %s, (undirected) m = %s\n",
		integer_to_string(n).c_str(),
		integer_to_string(m / 2).c_str());

	pstart = new ui[n + 1];
	edges = new ui[m];
	reverse = new ui[m];
	ui j = 0;
	for (ui i = 0; i < n; i++) {
		pstart[i] = j;
		while (j < m && epairs[j].first == nodes[i]) {
			edges[j] = idmp[epairs[j].second];
			reverse[j] = i;
			++j;
		}
	}
	pstart[n] = j;
}


int readRawDIM10Text(const char* filepath) {
	std::ifstream infile;
	const int SZBUF = 99999999;
	char *buf = new char[SZBUF];
	std::vector<std::pair<ui, ui> > epairs;
	std::vector<ui> nodes;
	//FILE *f = Utility::open_file(filepath, "r");
	infile.open(filepath, std::ios::in);
	if (!infile.is_open()) {
		fprintf(stderr, "can not find file %s\n", filepath);
		exit(1);
	}

	infile.getline(buf, SZBUF);
	while (buf[0] == '%') infile.getline(buf, SZBUF);

	std::stringstream ss(buf);
	int fmt = 0;
	ss >> n >> m >> fmt;
	if (fmt != 0){
		printf("Format of %s is not supported yet\n", filepath);
		exit(0);
	}
	m *= 2;
	pstart = new ui[n + 1];
	edges = new ui[m];
	reverse = new ui[m];
	ui j = 0;
	for (ui u = 0; u < n; u++) {
		pstart[u] = j;
		infile.getline(buf, SZBUF);
		std::stringstream ss(buf);
		int nei;
		while (ss >> nei) {
			//printf("%d ", nei);
			if ((nei - 1) != u) {
				edges[j] = nei - 1;
				reverse[j] = u;
				j++;
				//if (j==745)
				//	printf("pause\n");
			}
		}
		//printf("\n");
		std::sort(edges + pstart[u], edges + j);
	}
	pstart[n] = j;	
	assert(j == m);
	printf("n:%d m:%d\n",n,m/2);
	return 0;
}

int writeBinaryGraph(const char* filepath) {
	FILE *f = open_file(filepath, "wb");
	ui tt = sizeof(ui);
	fwrite(&tt, sizeof(ui), 1, f); //length of ui
	fwrite(&n, sizeof(ui), 1, f);
	fwrite(&m, sizeof(ui), 1, f);
	ui *degree = new ui[n];
	for (ui i = 0; i < n; i++)
		degree[i] = pstart[i + 1] - pstart[i];
	fwrite(degree, sizeof(ui), n, f);
	fwrite(edges, sizeof(ui), m, f);
	fclose(f);
	return 0;
}

void promot() {
	printf("text2bin textfile [binfile]\n");
}

int main(int argc, char** argv) {
	char txtFilePath[FILELEN] = "\0";
	char binFilePath[FILELEN] = "\0";
	if (argc < 2) {
		promot();
		exit(1);;
	}
	else if (argc == 2) {
		strncpy(txtFilePath, argv[1], FILELEN);
	}
	else if (argc == 3) {
		strncpy(txtFilePath, argv[1], FILELEN);
		strncpy(binFilePath, argv[2], FILELEN);
	}
	else {
		printf("Redundant args!\n");
	}
	//Identify format
	int pos = fileSuffixPos(txtFilePath);
	if (strcmp(txtFilePath + pos, "graph") == 0) {//dimacs 10
		readRawDIM10Text(txtFilePath);
	}
	else if (strcmp(txtFilePath + pos, "txt") == 0) { //snap
		readRawSNAPText(txtFilePath);
	}
	else if (strcmp(txtFilePath + pos, "clq") == 0) { //DIMACS2, clq file
		readDIMACS2Text(txtFilePath);
	}
	else {
		printf("Graph file not support\n");
		exit(2);
	}
	//Check ouput name
	if (strlen(binFilePath) == 0) {
		strncpy(binFilePath, txtFilePath, FILELEN);
		binFilePath[pos++] = 'b';
		binFilePath[pos++] = 'i';
		binFilePath[pos++] = 'n';
		binFilePath[pos++] = '\0';
	}
	//Write
	writeBinaryGraph(binFilePath);
	return 0;
}