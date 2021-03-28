#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <deque>
#include <algorithm>
#include <omp.h>
#include <time.h>
#define CacheLineSize 128

int n = 0;
unsigned int my = 0;
unsigned int mx = 0;
const int MaxSize = 280000;
using namespace std;
static struct Account {
    unsigned int x;
    unsigned int y;
} num[MaxSize];

static deque<unsigned int> Q[3000000];

void SplitString(string &s, vector<string> &v, const string &c) {
    int pos2;
    s += c;
    int size = s.size();
    for (int i = 0; i < size; i++) {
        pos2 = s.find(c, i);
        if (pos2 < size) {
            v.push_back(s.substr(i, pos2 - i));
            i = pos2 + c.size() - 1;
        }
    }
}

void readTxt(string file) {

    ifstream in(file, ios::in);
    if (!in.is_open()) {
        cout << "Error: opening file fail" << endl;
        exit(1);
    }
    string s;

    while (getline(in, s)) {
        vector<string> v;
        SplitString(s, v, ",");
        num[n].x = strtoul((const char *) v[0].c_str(), NULL, 10);
        num[n].y = strtoul((const char *) v[1].c_str(), NULL, 10);

        if (num[n].x > mx) {
            mx = num[n].x;
        }
        if (num[n].y > my) {
            my = num[n].y;
        }
        n++;
    }
    in.close();
}

struct ArcNode {
    unsigned int adjvex;
    ArcNode *next;
};

struct VertexNode {
    unsigned int vextex;
    ArcNode *firstedge;
};
int items = CacheLineSize / sizeof(ArcNode);
int visited[MaxSize] = {0};

class AlGraph {
public:
    AlGraph(Account nx[], int m, int e);

    ~AlGraph() {}

    void CheckCircle();

    void DFS(unsigned int x, int *visited, unsigned int *stack, int &top, bool *instack, int &count);

    void deleteArc(unsigned int x);

    void DeOrder(deque<unsigned int> *q, int count);

private:
    VertexNode adjlist[MaxSize];
    int vertexNum, arcNum;
};

AlGraph::AlGraph(Account nx[], int m, int e) {
    vertexNum = m;
    arcNum = e;
#pragma omp parallel for
    for (int i = 0; i <= vertexNum; i+=items) {
        for (int j = 0; j <items ; ++j) {
            adjlist[i+j].vextex = i+j;
            adjlist[i+j].firstedge = NULL;
        }
    }
#pragma omp parallel for num_threads(5)
    for (int k = 0; k < arcNum; k+=items) {
        for (int i = 0; i < items; ++i) {
            if(nx[k+i].x>m || nx[k+i].y>m){
                continue;
            }else{
                ArcNode *s = new ArcNode;
                s->adjvex = nx[k+i].y;

                ArcNode *p = adjlist[nx[k+i].x].firstedge;
                if (p == NULL || p->adjvex > s->adjvex) {
                    s->next = adjlist[nx[k+i].x].firstedge;
                    adjlist[nx[k+i].x].firstedge = s;
                } else {
                    while (p->next != NULL && p->next->adjvex < s->adjvex) {
                        p = p->next;
                    }
                    s->next = p->next;
                    p->next = s;
                }
            }
        }
    }
}

void AlGraph::deleteArc(unsigned int x) {
    ArcNode *ds = adjlist[x].firstedge;
    adjlist[x].firstedge = ds->next;
    delete[] ds;
}

void AlGraph::CheckCircle() {
    int count = 0;
    int top = -1;
    unsigned int stack[MaxSize];
    bool instack[MaxSize] = {false};
#pragma omp parallel for schedule(guided)
    //int items=CacheLineSize/sizeof(VertexNode);
    for (int i = 0; i <= vertexNum; i+=items) {
        for (int j = 0; j <items ; ++j) {
            if (adjlist[i+j].firstedge != NULL) {
                if (!visited[adjlist[i+j].firstedge->adjvex]) {
                    DFS(i+j, visited, stack, top, instack, count);
                    while (adjlist[i+j].firstedge->next != NULL) {
                        deleteArc(i+j);
                        if (!visited[adjlist[i+j].firstedge->adjvex])
                            DFS(i+j, visited, stack, top, instack, count);
                    }
                } else {
                    while (adjlist[i+j].firstedge->next != NULL) {
                        deleteArc(i+j);
                        if (!visited[adjlist[i+j].firstedge->adjvex])
                            DFS(i+j, visited, stack, top, instack, count);
                    }
                }
            }
            visited[i+j] = 1;
        }

    }
    DeOrder(Q, count);
}

void AlGraph::DFS(unsigned int x, int *visited, unsigned int *stack, int &top, bool *instack, int &count) {

    stack[++top] = x;
    instack[x] = true;
    if (top < 7) {
        if (adjlist[x].firstedge != NULL && !visited[adjlist[x].firstedge->adjvex]) {
            if (!instack[adjlist[x].firstedge->adjvex]) {
                DFS(adjlist[x].firstedge->adjvex, visited, stack, top, instack, count);
            } else if (stack[0] == adjlist[x].firstedge->adjvex && top >= 2) {
                count++;
                int t = 0;
                for (int j = t; j <= top; j++) {
                    Q[count].push_back(adjlist[stack[j]].vextex);

                }
            }
        }
    }
    if (top < 7 && top > 0) {
        if (adjlist[stack[top]].firstedge != NULL) {
            ArcNode *p = new ArcNode;
            p = adjlist[stack[top]].firstedge->next;
            while (p != NULL) {
                if (!visited[p->adjvex]) {
                    if (!instack[p->adjvex]) {
                        DFS(p->adjvex, visited, stack, top, instack, count);
                    } else if (stack[0] == p->adjvex && top >= 2) {
                        count++;
                        int t = 0;
                        int pos = top;
                        for (int j = t; j <= pos; j++) {
                            Q[count].push_back(adjlist[stack[j]].vextex);
                        }
                    }
                }
                p = p->next;
            }
            delete p;
        }
    }
    instack[stack[top--]] = false;
}

void WriteFile(ofstream &ofile, deque<unsigned int> q, int i) {
    for (int j = 0; j < i - 1; ++j) {
        ofile << q[j] << ",";
    }
    ofile << q[i - 1] << endl;
}


void AlGraph::DeOrder(deque<unsigned int> q[], int count) {
    string filename = "/projects/student/result.txt";
    ofstream ofile;
    ofile.open(filename);
    ofile << count << endl;
#pragma omp parallel for schedule(static)
    for (int i = 3; i <= 7; ++i) {
        for (int j = 1; j <= count; j+=items) {
            for (int k = 0; k <items ; ++k) {
                if (q[j+k].size() == i) {
                    WriteFile(ofile, q[j+k], i);
                }
            }
        }
    }
    ofile.close();
}

int main() {
    string file = "/data/test_data.txt";
    readTxt(file);
    int max = (mx > my ? my : mx);
    static AlGraph G(num, max, n);
    G.CheckCircle();
    return 0;
}
