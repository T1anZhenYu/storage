#define _CRT_SECURE_NO_WARNINGS
#include <cassert>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <sstream>
#include <stack>
#include <stdlib.h>
#include <string>
#include <vector>
using namespace std;

#define CACHESIZE 3
#define FILENUM 40
vector<float> filePopular(FILENUM, 0);
void Zipf(vector<float>& filePopular, int gamma)
{
    float sum = 0;
    for (int i = 1; i < FILENUM + 1; i++) {
        sum += pow(i, -gamma);
        // cout<<"sum:"<<sum<<endl;
    }
    // cout<<"sum:"<<sum<<endl;
    for (int i = 1; i < FILENUM + 1; i++) {
        filePopular[i - 1] = pow(i, -gamma) / sum;
    }
    for (int i = 0; i < FILENUM; i++) {
        cout << filePopular[i] << " ";
    }
    cout << endl;
}

class Vertex {
private:
    vector<int> cacheLRU;
    void cachefileLRU()
    {
        float randnum, temp;
        int pos;
        int flag = 0;
        // cout<<"vertex num:"<<this->val<<endl;
        for (int i = 0; i < CACHESIZE; i++) {
            temp = 0;
            flag = 0;
            randnum = (rand() % 100);

            pos = 0;
            while (temp < randnum) {
                temp += filePopular[pos] * 100;
                pos++;
            }
            cout << "randnum:" << randnum << " temp:" << temp << " pos:" << pos - 1 << endl;
            pos--;
            for (int j = 0; j < i; j++) {
                if (cacheLRU[j] == pos) {
                    flag = 1;
                }
            }
            if (flag == 1) {
                i--;
            } else {
                cacheLRU[i] = pos;
            }
        }
    }

public:
    int val;

    Vertex()
    {

        cacheLRU.resize(CACHESIZE);
        cachefileLRU();
    }
    Vertex(int val)
    {
        this->val = val;

        cacheLRU.resize(CACHESIZE);
        cachefileLRU();
    }
    ~Vertex() {}
    void showCache()
    {
        cout << "Vertex:" << this->val << " cache:";
        for (int i = 0; i < CACHESIZE; i++) {
            cout << this->cacheLRU[i] << " ";
        }
        // cout << endl;
    }
};
vector<Vertex> vertexs;
template <typename Weight>
class Edge {
private:
    Vertex *a, *b;
    Weight weight;

public:
    Edge(int a, int b, Weight weight)
    {
        this->a = &vertexs[a];
        this->b = &vertexs[b];
        this->a->val = a;
        this->b->val = b;
        this->weight = weight;
    }
    Edge() {}
    ~Edge() {}
    int v() { return a->val; }
    int w() { return b->val; }
    Weight wt() { return weight; }
    int other(int x)
    {
        assert(x == a->val || x == b->val);
        return x == a->val ? b->val : a->val;
    }
    friend ostream& operator<<(ostream& os, const Edge& e)
    {
        os << e.a->val << "-" << e.b->val << ": " << e.weight << "\n";
        return os;
    }
    bool operator<(Edge<Weight>& e)
    {
        return weight < e.wt();
    }
    bool operator<=(Edge<Weight>& e)
    {
        return weight <= e.wt();
    }
    bool operator>(Edge<Weight>& e)
    {
        return weight > e.wt();
    }
    bool operator>=(Edge<Weight>& e)
    {
        return weight >= e.wt();
    }
    bool operator==(Edge<Weight>& e)
    {
        return weight == e.wt();
    }
};

template <typename Item>
class MinHeap {

private:
    Item* data;
    int count;
    int capacity;
    void shiftUp(int k)
    {
        while (k > 1 && data[k / 2] > data[k]) {
            swap(data[k / 2], data[k]);
            k /= 2;
        }
    }
    void shiftDown(int k)
    {
        while (2 * k <= count) {
            int j = 2 * k;
            if (j + 1 <= count && data[j + 1] < data[j])
                j++;
            if (data[k] <= data[j])
                break;
            swap(data[k], data[j]);
            k = j;
        }
    }

public:
    MinHeap(int capacity)
    {
        data = new Item[capacity + 1];
        count = 0;
        this->capacity = capacity;
    }
    MinHeap(Item arr[], int n)
    {
        data = new Item[n + 1];
        capacity = n;
        for (int i = 0; i < n; i++)
            data[i + 1] = arr[i];
        count = n;
        for (int i = count / 2; i >= 1; i--)
            shiftDown(i);
    }
    ~MinHeap()
    {
        delete[] data;
    }
    int size()
    {
        return count;
    }
    bool isEmpty()
    {
        return count == 0;
    }
    void insert(Item item)
    {
        assert(count + 1 <= capacity);
        data[count + 1] = item;
        shiftUp(count + 1);
        count++;
    }
    Item extractMin()
    {
        assert(count > 0);
        Item ret = data[1];
        swap(data[1], data[count]);
        count--;
        shiftDown(1);
        return ret;
    }
    Item getMin()
    {
        assert(count > 0);
        return data[1];
    }
    void show()
    {
        cout << "| ";
        for (int i = 1; i <= count; i++)
            cout << data[i]->wt() /*<<","<<data[i]*/ << " | ";
        cout << endl;
    }
};
// 稠密图 - 邻接矩阵
template <typename Weight>
class DenseGraph {
private:
    int n, m;
    bool directed;
    vector<vector<Edge<Weight>*>> g;

public:
    DenseGraph(int n, bool directed)
    {
        this->n = n;
        this->m = 0;
        this->directed = directed;
        //构建有n个点的图，每个点的边的个数现在是0
        for (int i = 0; i < n; i++) {
            g.push_back(vector<Edge<Weight>*>(n, NULL));
        }
    }
    ~DenseGraph()
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (g[i][j] != NULL)
                    delete g[i][j];
    }
    int V() { return n; }
    int E() { return m; }
    void addEdge(int v, int w, Weight weight)
    {
        assert(v >= 0 && v < n);
        assert(w >= 0 && w < n);
        if (hasEdge(v, w)) {
            delete g[v][w];
            if (!directed)
                delete g[w][v];
            m--;
        }
        g[v][w] = new Edge<Weight>(v, w, weight);
        if (!directed)
            g[w][v] = new Edge<Weight>(w, v, weight);
        m++;
    }
    bool hasEdge(int v, int w)
    {
        assert(v >= 0 && v < n);
        assert(w >= 0 && w < n);
        return g[v][w] != NULL;
    }
    void show()
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                if (g[i][j])
                    cout << g[i][j]->wt() << "\t";
                else
                    cout << "NULL\t";
            cout << endl;
        }
    }
    class adjIterator {
    private:
        DenseGraph& G;
        int v;
        int index;

    public:
        //这个构造函数没看懂
        adjIterator(DenseGraph& graph, int v)
            : G(graph)
        {
            this->v = v;
            this->index = -1;
        }
        Edge<Weight>* begin()
        {
            index = -1;
            return next();
        }
        Edge<Weight>* next()
        {
            for (index += 1; index < G.V(); index++)
                if (G.g[v][index])
                    return G.g[v][index];
            return NULL;
        }
        bool end()
        {
            return index >= G.V();
        }
    };
};

template <typename Weight>
class SparseGraph {
private:
    int n, m;
    bool directed;
    vector<vector<Edge<Weight>*>> g;

public:
    SparseGraph(int n, bool directed)
    {
        this->n = n;
        this->m = 0;
        this->directed = directed;
        for (int i = 0; i < n; i++)
            g.push_back(vector<Edge<Weight>*>());
    }
    ~SparseGraph()
    {

        for (int i = 0; i < n; i++)
            for (int j = 0; j < g[i].size(); j++)
                delete g[i][j];
    }
    int V() { return n; }
    int E() { return m; }
    void addEdge(int v, int w, Weight weight)
    {
        assert(v >= 0 && v < n);
        assert(w >= 0 && w < n);

        g[v].push_back(new Edge<Weight>(v, w, weight));
        if (v != w && !directed)
            g[w].push_back(new Edge<Weight>(w, v, weight));
        m++;
    }
    bool hasEdge(int v, int w)
    {
        assert(v >= 0 && v < n);
        assert(w >= 0 && w < n);
        for (int i = 0; i < g[v].size(); i++)
            if (g[v][i]->other(v) == w)
                return true;
        return false;
    }
    void show()
    {

        for (int i = 0; i < n; i++) {
            cout << "vertex " << i << ":\t";
            for (int j = 0; j < g[i].size(); j++)
                cout << "( to:" << g[i][j]->w() << ",wt:" << g[i][j]->wt() << ")\t";
            cout << endl;
        }
    }
    class adjIterator {
    private:
        SparseGraph& G;
        int v;
        int index;

    public:
        adjIterator(SparseGraph& graph, int v)
            : G(graph)
        {
            this->v = v;
            this->index = 0;
        }
        Edge<Weight>* begin()
        {
            index = 0;
            if (G.g[v].size())
                return G.g[v][index];
            return NULL;
        }
        Edge<Weight>* next()
        {
            index += 1;
            if (index < G.g[v].size())
                return G.g[v][index];
            return NULL;
        }
        bool end()
        {
            return index >= G.g[v].size();
        }
    };
};

template <typename Graph, typename Weight>
class ReadGraph {
public:
    ReadGraph(Graph& graph, const string& filename)
    {
        ifstream file(filename);
        string line;
        int V, E;
        assert(file.is_open());
        assert(getline(file, line));
        stringstream ss(line);
        ss >> V >> E;
        assert(graph.V() == V);
        for (int i = 0; i < E; i++) {
            assert(getline(file, line));
            stringstream ss(line);
            int a, b;
            Weight w;
            ss >> a >> b >> w;
            assert(a >= 0 && a < V);
            assert(b >= 0 && b < V);
            graph.addEdge(a, b, w);
        }
    }
};

template <typename Graph, typename Weight>
class LazyPrimMST {
private:
    Graph& G;
    MinHeap<Edge<Weight>> pq;
    bool* marked;
    vector<Edge<Weight>> mst;
    Weight mstWeight;
    void visit(int v)
    {
        assert(!marked[v]);
        marked[v] = true;
        typename Graph::adjIterator adj(G, v);
        for (Edge<Weight>* e = adj.begin(); !adj.end(); e = adj.next())
            if (!marked[e->other(v)])
                pq.insert(*e);
    }

public:
    LazyPrimMST(Graph& graph)
        : G(graph)
        , pq(MinHeap<Edge<Weight>>(graph.E()))
    {
        marked = new bool[G.V()];
        for (int i = 0; i < G.V(); i++)
            marked[i] = false;
        mst.clear();
        // Lazy Prim
        visit(0);
        while (!pq.isEmpty()) {
            Edge<Weight> e = pq.extractMin();
            if (marked[e.v()] == marked[e.w()])
                continue;
            mst.push_back(e);
            if (!marked[e.v()])
                visit(e.v());
            else
                visit(e.w());
        }
        mstWeight = mst[0].wt();
        for (int i = 1; i < mst.size(); i++)
            mstWeight += mst[i].wt();
    }
    ~LazyPrimMST()
    {
        delete[] marked;
    }
    vector<Edge<Weight>> mstEdges()
    {
        return mst;
    };
    Weight result()
    {
        return mstWeight;
    };
};

int main(void)
{
    string filename = "testG1.txt";
    int V = 8;
    Zipf(filePopular, 1);
    vertexs.resize(V);
    // cout<<"heer"<<endl;
    for (int i = 0; i < V; i++) {
        vertexs[i].val = i;
        vertexs[i].showCache();
        cout<<endl;
    }
    SparseGraph<double> g = SparseGraph<double>(V, false);
    ReadGraph<SparseGraph<double>, double> readGraph(g, filename);
    // cout << "Test Lazy Prim MST:" << endl;
    // LazyPrimMST<SparseGraph<double>, double> lazyPrimMST(g);
    // vector<Edge<double>> mst = lazyPrimMST.mstEdges();
    // for (int i = 0; i < mst.size(); i++)
    //     cout << mst[i] << endl;
    // cout << "The MST weight is: " << lazyPrimMST.result() << endl;
    // cout << endl;
    // cin.get();
    
    return 0;
}
