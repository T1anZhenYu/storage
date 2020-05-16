#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <cmath>
#include <errno.h>
#include <iostream>
#include <iterator>
#include <list>
#include <queue>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <vector>
using namespace std;
#define MAX_LINE_LEN 4000
#define MAX_VERTEX_NUM 600 //图中最大顶点个数
#define FILENUM 40 //文件数目
#define CACHESIZE 2 //顶点缓存空间
#define MYINFINITY 1000000 //将MYINFINITY定义为无穷大的值
#define BACKHUAL 20 //回传时延
vector<float> filePopular(FILENUM, 0); //各个文件的流行度
vector<float> soft2TightTimeLim(FILENUM, 0); //从松到紧的时延要求
vector<float> tight2SoftTimeLim(FILENUM, 0); //从紧到送的时延要求
vector<float> equalTimeLim(FILENUM, 1 / FILENUM); //处处相等的时延要求
vector<vector<int>> Distance; //图中各点到其他点的最短距离

/*******************************************************
 * 函数名称 timelimtGenerate()
 * 作用说明 生成三种时延要求
 * 参数列表 void
 * 返回值 void
********************************************************/
void timelimtGenerate()
{
    int dispart = int(FILENUM / (BACKHUAL / 2));
    int dis = BACKHUAL;
    for (int i = 0; i < FILENUM; i++) {
        soft2TightTimeLim[i] = dis;
        tight2SoftTimeLim[FILENUM - 1 - i] = dis;
        equalTimeLim[i] = BACKHUAL / 2;
        if (i % dispart == 0) {
            dis--;
        }
    }
    cout << "soft2Tight" << endl;
    for (int i = 0; i < FILENUM; i++) {
        cout << soft2TightTimeLim[i] << " ";
    }
    cout << endl;

    cout << "tight2Soft" << endl;
    for (int i = 0; i < FILENUM; i++) {
        cout << tight2SoftTimeLim[i] << " ";
    }
    cout << endl;
}
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
}
//保存每个顶点信息的数据结构
class GraphNode {
private:
    vector<int> cache;
    // LRU 的缓存方式，在cache中缓存文件
    void cachefileLRU();

public:
    bool known; //当前顶点距离起点的距离是否确定
    int dist; //当前顶点到起点的最短距离
    int path; //当前顶点距离起点的最短路径的前一个顶点
    GraphNode()
    {
        cache.resize(CACHESIZE);
        cachefileLRU();
    }
    void showCache()
    {
        // cout << "Vertex:" << this->val << " cache:";
        for (int i = 0; i < CACHESIZE; i++) {
            cout << this->cache[i] << " ";
        }
        cout << endl;
    }
    int cacheHit(int fileid)
    {
        for (int i = 0; i < CACHESIZE; i++) {
            // cout<<"cache "<<this->cache[i]<<" fileid:"<<fileid<<endl;
            if (this->cache[i] == fileid) {
                return 1;
            }
        }
        return 0;
    }
};
/*************************************************
 * 函数名称：cachefileLRU();
 * 功能描述：按LRU的方式缓存
 * 参数列表：void
 * 返回结果：void
**************************************************/
void GraphNode::cachefileLRU()
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
        pos--;
        for (int j = 0; j < i; j++) {
            if (cache[j] == pos) {
                flag = 1;
            }
        }
        if (flag == 1) {
            i--;
        } else {
            cache[i] = pos;
        }
    }
}

vector<GraphNode> nodeArr; //保存每个顶点信息的数组
//图节点信息
typedef struct Node {
    int edge_num; //边号
    int src; //源点
    int vertex; //自身
    int weight; //边的权重
} Node;

/*******************************************************
*  类名称： 邻接表图
********************************************************/
class Graph {
private:
    int edge_num; //图边的个数
    int vertex_num; //图的顶点数目
    list<Node>* graph_list; //邻接表

public:
    Graph() {}
    Graph(char* graph[], int edgenum);
    ~Graph();
    void print();
    void unwightedShorestPath(int src); //算法1求最短距离
    void unwightedShorestPathAdv(int src); //算法2求最短距离
    void printShorestPath(int src); //输出顶点src到各个顶点最短距离的信息
    int getVertexNum(); //获取顶点个数
    void cacheInDistance(float& HitS2T, float& HitT2S, float& HitEqual); //计算在规定范围内，缓存命中的概率

private:
    vector<int> get_graph_value(char* graph[], int columns);
    void addEdge(char* graph[], int columns);
};
/*************************************************
 * 函数名称：cacheInDistance();
 * 功能描述：计算在规定范围内，缓存命中的概率
 * 参数列表：void
 * 返回结果：命中概率
**************************************************/
void Graph::cacheInDistance(float& HitS2T, float& HitT2S, float& HitEqual)
{
    /*算法描述：
        for(file in allFIle):
            for(server in allServer):
                for(neibor in allNeibor):
                    if(neiborDis < fileTimelim):
                        for(cache in allCache):
                            if(cache == file):
                                hit++;
    */
    float hitS2T, hitT2S, hitEqual;
    hitS2T = 0;
    hitT2S = 0;
    hitEqual = 0;
    for (int f = 0; f < FILENUM; f++) {
        int fhitS2T = 0;
        int fhitT2S = 0;
        int fhitEqual = 0;

        for (int ser = 1; ser < this->getVertexNum(); ser++) {
            int flagS2T = 0;
            int flagT2S = 0;
            int flagEqual = 0;
            for (int nei = 1; nei < this->getVertexNum(); nei++) {
                int flag = nodeArr[nei].cacheHit(f);
                if (Distance[ser][nei] < soft2TightTimeLim[f] && flag) {

                    if (!flagS2T) {
                        fhitS2T++;
                        flagS2T = 1;
                    }
                }
                if (Distance[ser][nei] < tight2SoftTimeLim[f] && flag) {
                    if (!flagT2S) {
                        fhitT2S++;
                        flagT2S = 1;
                    }
                }
                if (Distance[ser][nei] < tight2SoftTimeLim[f] && flag) {
                    if (!flagEqual) {
                        fhitEqual++;
                        flagEqual = 1;
                    }
                }
            }
        }
        hitS2T += (float(fhitS2T) / float(this->getVertexNum() - 1));
        hitT2S += (float(fhitT2S) / float(this->getVertexNum() - 1));
        hitEqual += (float(fhitEqual) / float(this->getVertexNum() - 1));
        // cout<<"f:"<<f<<" hitS2t:"<<hitS2T<<" hitT2S"<<hitT2S<<endl;
    }

    HitEqual = hitEqual / FILENUM;
    HitS2T = hitS2T / FILENUM;
    HitT2S = hitT2S / FILENUM;
}

/*************************************************
 * 函数名称：getVertexNum();
 * 功能描述：获取图的顶点个数
 * 参数列表：void
 * 返回结果：图的顶点个数
**************************************************/
int Graph::getVertexNum()
{
    return vertex_num;
}

/*************************************************
*  函数名称：unwightedShorestPathAdv(int src)
*  功能描述：求无权图的任意点到其它顶点的距离,
*            该算法比unwightedShorestPathAdv要好
*  参数列表：src是起点
*  返回结果：void 
*************************************************/
void Graph::unwightedShorestPathAdv(int src)
{
    queue<int> que;

    //步骤1，将各个顶点的dist设置为无穷大
    for (int i = 0; i < vertex_num; ++i)
        nodeArr[i].dist = MYINFINITY;

    //步骤2，将顶点src的dist字段设为0
    nodeArr[src].dist = 0;

    //步骤3，将顶点src压入队列que中
    que.push(src);

    //步骤5
    while (!que.empty()) {
        //步骤4
        int top = que.front(); //获得队列的首元素
        que.pop(); //弹出队列的首元素

        //遍历与顶点相邻的所有顶点
        for (list<Node>::iterator it = graph_list[top].begin(); it != graph_list[top].end(); ++it) {
            if (nodeArr[(*it).vertex].dist == MYINFINITY) {
                nodeArr[(*it).vertex].dist = nodeArr[top].dist + 1;
                nodeArr[(*it).vertex].path = top;

                que.push((*it).vertex);
            }
        }
    }
}

/*************************************************
*  函数名称：unwightedShorestPath(int src)
*  功能描述：求无权图的任意点到其它顶点的距离
*  参数列表：src是起点
*  返回结果：void 
*************************************************/
void Graph::unwightedShorestPath(int src)
{
    //步骤1，初始化配置表
    for (int i = 0; i < vertex_num; ++i) {
        nodeArr[i].known = false;
        nodeArr[i].dist = MYINFINITY;
        nodeArr[i].path = 0;
    }

    //步骤2，先把距离为0的顶点的dist设置为0
    nodeArr[src].dist = 0;

    //步骤4
    for (int currentDist = 0; currentDist < vertex_num; ++currentDist) {
        //步骤3
        for (int i = 0; i < vertex_num; ++i) {
            if ((!nodeArr[i].known) && (nodeArr[i].dist == currentDist)) {
                nodeArr[i].known = true;

                //遍历与顶点i相邻的所有顶点
                for (list<Node>::iterator it = graph_list[i].begin(); it != graph_list[i].end(); ++it) {
                    if (nodeArr[(*it).vertex].dist == MYINFINITY) {
                        nodeArr[(*it).vertex].dist = currentDist + 1;
                        nodeArr[(*it).vertex].path = i;
                    }
                }
            }
        }
    }
}

/*************************************************
*  函数名称：printShorestPath()
*  功能描述：将获得的src顶点到其它顶点的最短路径输出
*  参数列表：无
*  返回结果：无
*************************************************/
void Graph::printShorestPath(int src)
{
    for (int i = 0; i < vertex_num; ++i) {
        if (nodeArr[i].known) {
            // cout << i << "\t" << nodeArr[i].known << "\t" << nodeArr[i].dist << "\t" << nodeArr[i].path << "\t";
            Distance[src][i] = nodeArr[i].dist;
            // nodeArr[i].showCache();
        }
    }
}

/*************************************************
*  函数名称：print
*  功能描述：将图的信息以邻接表的形式输出到标准输出
*  参数列表：无
*  返回结果：无
*************************************************/
void Graph::print()
{
    cout << "******************************************************************" << endl;
    //for(int i = 0 ; i < MAX_VERTEX_NUM; ++i){
    for (int i = 0; i < vertex_num; ++i) {
        if (graph_list[i].begin() != graph_list[i].end()) {
            cout << i << "-->";
            for (list<Node>::iterator it = graph_list[i].begin(); it != graph_list[i].end(); ++it) {
                cout << (*it).vertex << "(边号:" << (*it).edge_num << ",权重:" << (*it).weight << ")-->";
            }
            cout << "NULL" << endl;
        }
    }

    cout << "******************************************************************" << endl;
}

/*************************************************
*  函数名称：get_graph_value
*  功能描述：将图的每一条边的信息保存到一个数组中
*  参数列表: graph：指向图信息的二维数组
             columns:图的第几条边
*  返回结果：无
*************************************************/
vector<int> Graph::get_graph_value(char* graph[], int columns)
{
    vector<int> v;
    char buff[20];
    int i = 0, j = 0, val;
    memset(buff, 0, 20);

    while ((graph[columns][i] != '\n') && (graph[columns][i] != '\0')) {
        if (graph[columns][i] != ',') {
            buff[j] = graph[columns][i];
            j++;
        } else {
            j = 0;
            val = atoi(buff);
            v.push_back(val);
            memset(buff, 0, 20);
        }
        i++;
    }
    val = atoi(buff);
    v.push_back(val);

    return v;
}

/*************************************************
*  函数名称：addEdge
*  功能描述：将图的每一条边的信息加入图的邻接表中
*  参数列表：graph：指向图信息的二维数组
             columns:图的第几条边
*  返回结果：无
*************************************************/
void Graph::addEdge(char* graph[], int columns)
{
    Node node;
    vector<int> v = get_graph_value(graph, columns);

    node.edge_num = v[0];
    node.src = v[1];
    node.vertex = v[2];
    node.weight = v[3];

    //根据顶点的标号，求的总的顶点数目
    if (node.vertex > vertex_num)
        vertex_num = node.vertex;

    //要考虑重复的边，但是边的权重不一样
    for (list<Node>::iterator it = graph_list[node.src].begin(); it != graph_list[node.src].end(); ++it) {
        if ((*it).vertex == node.vertex) {
            if ((*it).weight > node.weight) {
                (*it).weight = node.weight;
            }
            return;
        }
    }

    graph_list[node.src].push_back(node);
}

/*************************************************
*  函数名称：构造函数
*  功能描述：以邻接表的形式保存图的信息,并保存必须经过的顶点
*  参数列表：graph：指向图信息的二维数组
             edgenum:图的边的个数
*  返回结果：无
*************************************************/
Graph::Graph(char* graph[], int edgenum)
{
    edge_num = edgenum;
    vertex_num = 0;
    graph_list = new list<Node>[MAX_VERTEX_NUM + 1];

    for (int i = 0; i < edgenum; ++i) {
        addEdge(graph, i);
    }

    vertex_num++;
}

/*************************************************
*  函数名称：析构函数
*  功能描述：释放动态分配的内存
*  参数列表：无
*  返回结果：无
*************************************************/
Graph::~Graph()
{
    delete[] graph_list;
}

#endif

/****************************************************************
*   函数名称：read_file
*   功能描述: 读取文件中的图的数据信息
*   参数列表: buff是将文件读取的图信息保存到buff指向的二维数组中 
*             spec是文件中图最大允许的边的个数
*             filename是要打开的图文件
*   返回结果：无
*****************************************************************/
int read_file(char** const buff, const unsigned int spec, const char* const filename)
{
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Fail to open file %s, %s.\n", filename, strerror(errno));
        return 0;
    }
    printf("Open file %s OK.\n", filename);

    char line[MAX_LINE_LEN + 2];
    unsigned int cnt = 0;
    while ((cnt < spec) && !feof(fp)) {
        line[0] = 0;
        fgets(line, MAX_LINE_LEN + 2, fp);
        if (line[0] == 0)
            continue;
        buff[cnt] = (char*)malloc(MAX_LINE_LEN + 2);
        strncpy(buff[cnt], line, MAX_LINE_LEN + 2 - 1);
        buff[cnt][4001] = 0;
        cnt++;
    }
    fclose(fp);
    printf("There are %d lines in file %s.\n", cnt, filename);

    return cnt;
}

/****************************************************************
*   函数名称：release_buff
*   功能描述: 释放刚才读取的文件中的图的数据信息
*   参数列表: buff是指向文件读取的图信息
*             valid_item_num是指图中边的个数
*   返回结果：void
*****************************************************************/
void release_buff(char** const buff, const int valid_item_num)
{
    for (int i = 0; i < valid_item_num; i++)
        free(buff[i]);
}

int main(int argc, char* argv[])
{
    Zipf(filePopular, 1);
    timelimtGenerate();
    nodeArr.resize(MAX_VERTEX_NUM);
    char* topo[5000];
    int edge_num;
    char* demand;
    int demand_num;

    char* topo_file = "undirectedGraph.txt";

    edge_num = read_file(topo, 5000, topo_file);
    if (edge_num == 0) {
        printf("Please input valid topo file.\n");
        return -1;
    }

    int src = 1;

    Graph G(topo, edge_num);
    G.print();
    int vertex_num = G.getVertexNum();
    Distance.resize(vertex_num);
    for (int i = 0; i < vertex_num; i++) {
        Distance[i].resize(vertex_num);
    }

    // cout << "算法1: " << endl;
    for (int s = 1; s < vertex_num; s++) {
        G.unwightedShorestPath(s);
        G.printShorestPath(s);
    }
    cout << endl
         << "邻接矩阵" << endl;
    cout << "#######################################" << endl;
    for (int i = 1; i < vertex_num; i++) {
        cout << "vertex:" << i << "| ";
        for (int j = 1; j < vertex_num; j++) {
            cout << Distance[i][j] << " ";
        }
        cout << " cached file: ";
        nodeArr[i].showCache();

        // cout << endl;
    }
    float hitS2T, hitT2S, hitEqual;
    G.cacheInDistance(hitS2T, hitT2S, hitEqual);
    cout << endl
         << "Hit in Distance" << endl;
    cout << "#######################################" << endl;
    cout << hitS2T << "\t" << hitT2S << "\t" << hitEqual << endl;
    release_buff(topo, edge_num);

    return 0;
}