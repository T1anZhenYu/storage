#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <cmath>
#include <errno.h>
#include <iomanip>
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
#define FILENUM 400 //文件数目
#define CACHESIZE 2 //顶点缓存空间
#define MYINFINITY 1000000 //将MYINFINITY定义为无穷大的值
#define BACKHUAL 20 //回传时延
#define LAMBDA 0.05 //服务器密度

vector<float> QF(FILENUM, 0); //各个文件的流行度
vector<float> soft2TightTimeLim(FILENUM, 0); //从松到紧的时延要求
vector<float> tight2SoftTimeLim(FILENUM, 0); //从紧到送的时延要求
vector<float> equalTimeLim(FILENUM, 1 / FILENUM); //处处相等的时延要求
vector<int> facNum(10, 1); //存储10个阶乘结果，防止反复计算
vector<float> LRUFileProb(FILENUM, 0); //LRUCache
vector<float> tempTLACFileProbS2T(FILENUM, 0); //PfS2T
vector<float> tempTLACFileProbT2S(FILENUM, 0); //PfT2S
vector<float> tempTLACFileProbEqual(FILENUM, 0); //PfEqual
vector<float> maxTLACFileProbS2T(FILENUM, 0); //max PfS2T
vector<float> maxTLACFileProbT2S(FILENUM, 0); //max PfT2S
vector<float> maxTLACFileProbEqual(FILENUM, 0); //max PfEqual
vector<vector<int>> Distance; //图中各点到其他点的最短距离
/*******************************************************
 * 函数名称 fac()
 * 作用说明 提前计算好阶乘，防止反复计算
 * 参数列表 void
 * 返回值 void
********************************************************/
void fac()
{
    int res = 1;
    for (int i = 1; i < facNum.size(); i++) {
        facNum[i] = facNum[i - 1] * (i + 1);
    }
    cout << "fac :\n";
    for (int i = 0; i < facNum.size(); i++) {
        cout << facNum[i] << " ";
    }
    cout << endl;
}
/*******************************************************
 * 函数名称 timelimtGenerate()
 * 作用说明 生成三种时延要求
 * 参数列表 void
 * 返回值 void
*******************************************************/
void timelimtGenerate()
{
    int dispart = int(FILENUM / 20);
    int dis = BACKHUAL + 2;
    for (int i = 0; i < FILENUM; i++) {
        soft2TightTimeLim[i] = dis;
        tight2SoftTimeLim[FILENUM - 1 - i] = soft2TightTimeLim[i];
        equalTimeLim[i] = BACKHUAL / 3;
        if (i % 20 == 19) {
            dis = dis - 1;
        }
    }
    cout << "soft2Tight" << endl;
    for (int i = 0; i < FILENUM; i++) {
        cout << setw(4)<<soft2TightTimeLim[i] << " ";
    }
    cout << endl;

    cout << "tight2Soft" << endl;
    for (int i = 0; i < FILENUM; i++) {
        cout<< setw(4)<<tight2SoftTimeLim[i] << " ";
    }
    cout << endl;
}
void Zipf(vector<float>& QF, int gamma)
{
    float sum = 0;
    float temp = 0;
    for (int i = 1; i < FILENUM + 1; i++) {

        if (i % 20 == 1) {
            temp = pow(i, -gamma);
        }
        sum += temp;
        QF[i - 1] = temp;
        // cout<<"sum:"<<sum<<endl;
    }
    cout << "QF:" << sum << endl;
    for (int i = 1; i < FILENUM + 1; i++) {
        float temp = QF[i - 1];
        QF[i - 1] = temp / sum;
        cout <<setw(4)<< " " << QF[i - 1];
    }
    cout << endl;
}

/*******************************************************
 * 函数名称 taylorLW()
 * 作用说明 朗伯W函数在x=0处的泰勒展开 展开级数暂定为7级
 * 参数列表 需要估计的位置
 * 返回值 估计值
********************************************************/
float taylorLW(float x, float maxNum)
{
    float res = 0;
    // cout << "in taylorLw x:"<<x << endl;
    /************************************************************************
     * 下面的估算方法只在x>3有用
    ************************************************************************/
    // float L1 = log(x);
    // float L2 = log(log(x));

    // float A = L1 - L2;
    // float B = L2 / L1;
    // float C = L2 * (-2 + L2) / (2 * pow(L1, 2));
    // float D = L2 * (6 - 9 * L2 + 2 * pow(L2, 2)) / (6 * pow(L1, 3));
    // float E = L2 * (-12 + 36 * L2 - 22 * pow(L2, 2) + 3 * pow(L2, 3)) / (12 * pow(L1, 4));
    // float F = L2 * (60 - 300 * L2 + 350 * pow(L2, 2) - 125 * pow(L2, 3) + 12 * pow(L2, 4))
    // / (60 * pow(L1, 5));
    // res = A + B + C + D + E + F;
    /************************************************************************
     * 下面的估算方法是一个比较粗略的估计
    ************************************************************************/
    float minNum = 1.0;
    res = log(x + 1) / log(3.6);
    // res = min(res, maxNum);
    // res = max(res, minNum);
    return res;
}

//保存每个顶点信息的数据结构
class GraphNode {
private:
    vector<int> cacheLRU; //LRU 或 RCS的各个文件的缓存概率
    vector<int> cacheTLACS2T; //TLAC 的各个文件的缓存概率 在S2T下
    vector<int> cacheTLACT2S; //TLAC 的各个文件的缓存概率 在T2S下
    vector<int> cacheTLACEqual; //TLAC 的各个文件的缓存概率 在Equal下
    // LRU 的缓存方式，在cache中缓存文件

public:
    bool known; //当前顶点距离起点的距离是否确定
    int dist; //当前顶点到起点的最短距离
    int path; //当前顶点距离起点的最短路径的前一个顶点
    void cachefileLRU(); //根据LRU下各个文件的缓存概率，缓存文件
    void cachefileTLAC(); //根据TLAC下各个文件的缓存概率，缓存文件
    GraphNode()
    {
        cacheLRU.resize(CACHESIZE);
        cacheTLACS2T.resize(CACHESIZE);
        cacheTLACT2S.resize(CACHESIZE);
        cacheTLACEqual.resize(CACHESIZE);
    }
    void showCacheLRU()
    {
        // cout << "Vertex:" << this->val << " cache:";
        for (int i = 0; i < CACHESIZE; i++) {
            cout << setw(4) << this->cacheLRU[i] << " ";
        }
        // cout << endl;
    }
    void showCacheTLAC()
    {
        cout << "\t|TLAC S2T: ";
        for (int i = 0; i < CACHESIZE; i++) {
            cout << setw(4) << this->cacheTLACS2T[i] << " ";
        }
        cout << "\t|T2S:";
        for (int i = 0; i < CACHESIZE; i++) {
            cout << setw(4) << this->cacheTLACT2S[i] << " ";
        }
        cout << "\t|Equal:";
        for (int i = 0; i < CACHESIZE; i++) {
            cout << setw(4) << this->cacheTLACEqual[i] << " ";
        }
    }
    int cacheHit(int id, int filed) //0:LRU,1:RCS,2:TLAC S2T,3:TLAC T2S,4:TLAC Equal
    {
        if (filed == 0) {
            for (int i = 0; i < CACHESIZE; i++) {
                // cout<<"cache "<<this->cache[i]<<" id:"<<id<<endl;
                if (this->cacheLRU[i] == id) {
                    return 1;
                }
            }
        } else if (filed == 2) {
            for (int i = 0; i < CACHESIZE; i++) {
                // cout<<"cache "<<this->cache[i]<<" id:"<<id<<endl;
                if (this->cacheTLACS2T[i] == id) {
                    return 1;
                }
            }
        } else if (filed == 3) {
            for (int i = 0; i < CACHESIZE; i++) {
                // cout<<"cache "<<this->cache[i]<<" id:"<<id<<endl;
                if (this->cacheTLACT2S[i] == id) {
                    return 1;
                }
            }
        } else if (filed == 4) {
            for (int i = 0; i < CACHESIZE; i++) {
                // cout<<"cache "<<this->cache[i]<<" id:"<<id<<endl;
                if (this->cacheTLACEqual[i] == id) {
                    return 1;
                }
            }
        }
        return 0;
    }
};
/*************************************************
 * 函数名称：targetFunc();
 * 功能描述：给定Pf，求出目标函数的值
 * 参数列表：三种时延下的目标函数值的引用
 * 返回结果：void
**************************************************/
void targetFunc(float& DS2T, float& DT2S, float& DEqual)
{
    float sumS2T = 0;
    float sumT2S = 0;
    float sumEqual = 0;
    for (int i = 0; i < FILENUM; i++) {
        sumS2T += QF[i] * (1 - exp(-1 * M_PI * pow(soft2TightTimeLim[i], 2) * LAMBDA * tempTLACFileProbS2T[i]));
        sumT2S += QF[i] * (1 - exp(-1 * M_PI * pow(tight2SoftTimeLim[i], 2) * LAMBDA * tempTLACFileProbT2S[i]));
        sumEqual += QF[i] * (1 - exp(-1 * M_PI * pow(equalTimeLim[i], 2) * LAMBDA * tempTLACFileProbEqual[i]));
    }
    DS2T = sumS2T;
    DT2S = sumT2S;
    DEqual = sumEqual;
}
/*************************************************
 * 函数名称：cachefileTLAC();
 * 功能描述：按TLAC的方式缓存
 * 参数列表：void
 * 返回结果：void
**************************************************/
void GraphNode::cachefileTLAC()
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
            temp += maxTLACFileProbS2T[pos] * 100;
            pos++;
        }
        pos--;
        for (int j = 0; j < i; j++) {
            if (cacheTLACS2T[j] == pos) {
                flag = 1;
            }
        }
        if (flag == 1) {
            i--;
        } else {
            cacheTLACS2T[i] = pos;
        }
    }
    // cout<<"vertex num:"<<this->val<<endl;
    for (int i = 0; i < CACHESIZE; i++) {
        temp = 0;
        flag = 0;
        randnum = (rand() % 100);

        pos = 0;
        while (temp < randnum) {
            temp += maxTLACFileProbT2S[pos] * 100;
            pos++;
        }
        pos--;
        for (int j = 0; j < i; j++) {
            if (cacheTLACT2S[j] == pos) {
                flag = 1;
            }
        }
        if (flag == 1) {
            i--;
        } else {
            cacheTLACT2S[i] = pos;
        }
    }
    for (int i = 0; i < CACHESIZE; i++) {
        temp = 0;
        flag = 0;
        randnum = (rand() % 100);

        pos = 0;
        while (temp < randnum) {
            temp += maxTLACFileProbEqual[pos] * 100;
            pos++;
        }
        pos--;
        for (int j = 0; j < i; j++) {
            if (cacheTLACEqual[j] == pos) {
                flag = 1;
            }
        }
        if (flag == 1) {
            i--;
        } else {
            cacheTLACEqual[i] = pos;
        }
    }
}
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
            temp += LRUFileProb[pos] * 100;
            pos++;
        }
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
    void cacheInDistance(float& HitS2T, float& HitT2S, float& HitEqual, int filed); //计算在规定范围内，缓存命中的概率
    void timeLimAwareCacheProb(); //计算TLAC算法的各个文件的缓存概率
    void LRUCacheProb(); //计算LRU中各个文件的缓存概率
private:
    vector<int> get_graph_value(char* graph[], int columns);
    void addEdge(char* graph[], int columns);
};
/*************************************************
 * 函数名称：LRUCacheProb();
 * 功能描述：按TLAC的方式缓存
 * 参数列表：void
 * 返回结果：void
**************************************************/
void Graph::LRUCacheProb()
{
    float tempSum = 0;
    for (int i = 0; i < FILENUM; i++) {
        tempSum += QF[i];
    }
    for (int i = 0; i < FILENUM; i++) {
        LRUFileProb[i] = QF[i] / tempSum;
    }
    cout << "LRUCacheProb\n";
    for (int i = 0; i < FILENUM; i++) {
        if (i % 20 == 0) {
            cout << LRUFileProb[i] << " ";
        }
    }
    cout << endl;
}
/*************************************************
 * 函数名称：timeLimAwareCacheProb();
 * 功能描述：按TLAC的方式缓存
 * 参数列表：void
 * 返回结果：void
**************************************************/
void Graph::timeLimAwareCacheProb()
{
    float tempDS2T, tempDT2S, tempDEqual, maxDS2T, maxDT2S, maxDEqual;
    maxDS2T = 0;
    maxDT2S = 0;
    maxDEqual = 0;
    int B = 200;
    // for (int B = 0; B < 400; B++) {
    tempDS2T = 0;
    tempDT2S = 0;
    tempDEqual = 0;
    float sumfpS2T = 0;
    float sumfpT2S = 0;
    float sumfpEqual = 0;
    float minfpS2T = 0;
    float minfpT2S = 0;
    float minfpEqual = 0;
    for (int i = 0; i < FILENUM; i++) {
        float AS2T = M_PI * LAMBDA * pow(soft2TightTimeLim[i], 2);
        float AT2S = M_PI * LAMBDA * pow(tight2SoftTimeLim[i], 2);
        float AEqual = M_PI * LAMBDA * pow(equalTimeLim[i], 2);

        float a1S2T = log(QF[i]*float(B));
        // float a2S2T = taylorLW(float(B) / 1000.0 / CACHESIZE * a1S2T, AS2T + 1);
        float pfS2T = a1S2T / AS2T;

        float a1T2S = log(QF[i]*float(B));
        // float a2T2S = taylorLW(float(B) / 1000.0 / CACHESIZE * a1T2S, AT2S + 1);
        float pfT2S = a1T2S / AT2S;

        float a1Equal = log(QF[i]*float(B));
        // float a2Equal = taylorLW(float(B) / 1000.0 / CACHESIZE * a1Equal, AEqual + 1);
        float pfEqual = a1Equal / AEqual;
        // float a1T2S = exp((AT2S + 1) )/ QF[i];
        // float a2T2S = taylorLW(float(B) / 1000.0 / CACHESIZE * a1T2S, AT2S + 1);
        // float pfT2S = (-a2T2S + AT2S + 1) / AT2S;

        // float a1Equal = exp((AEqual + 1)) / QF[i];
        // float a2Equal = taylorLW(float(B) / 1000.0 / CACHESIZE * a1Equal, AEqual + 1);
        // float pfEqual = (-a2Equal + AEqual + 1) / AEqual;

        tempTLACFileProbS2T[i] = pfS2T;
        tempTLACFileProbT2S[i] = pfT2S;
        tempTLACFileProbEqual[i] = pfEqual;
        sumfpS2T += pfS2T;
        sumfpT2S += pfT2S;
        sumfpEqual += pfEqual;
        minfpS2T = min(minfpS2T,pfS2T);
        minfpT2S = min(minfpT2S, pfT2S);
        minfpEqual = min(minfpEqual, pfEqual);
        cout << "AEqual:" << AEqual << " a1Equal:" << a1Equal << " pfEqual" << pfEqual << " AT2S:" << AT2S << " a1T2S:" << a1T2S  << " pfT2S:" << pfT2S << endl;
    }
    sumfpS2T += minfpS2T*FILENUM;
    sumfpT2S += minfpT2S* FILENUM;
    sumfpEqual += minfpEqual*FILENUM;
    for (int i = 0; i < FILENUM; i++) {
        tempTLACFileProbS2T[i] = (tempTLACFileProbS2T[i]+minfpS2T) / sumfpS2T;
        tempTLACFileProbT2S[i] = (tempTLACFileProbT2S[i] +minfpT2S)/ sumfpT2S;
        tempTLACFileProbEqual[i] =( tempTLACFileProbEqual[i]+minfpEqual) / sumfpEqual;
    }
    targetFunc(tempDS2T, tempDT2S, tempDEqual);
    cout << "B:" << float(B) / 100.0 << "tempDS2T:" << tempDS2T << " tempDT2S:" << tempDT2S << " tempDEqual:" << tempDEqual << endl;
    if (tempDS2T > maxDS2T) {
        maxDS2T = tempDS2T;
        for (int i = 0; i < FILENUM; i++) {
            maxTLACFileProbS2T[i] = tempTLACFileProbS2T[i];
        }
    }
    if (tempDT2S > maxDT2S) {
        maxDT2S = tempDT2S;
        for (int i = 0; i < FILENUM; i++) {
            maxTLACFileProbT2S[i] = tempTLACFileProbT2S[i];
        }
    }
    if (tempDEqual > maxDEqual) {
        maxDEqual = tempDEqual;
        for (int i = 0; i < FILENUM; i++) {
            maxTLACFileProbEqual[i] = tempTLACFileProbEqual[i];
        }
    }
    // }

    cout << "maxTLACFileProbS2T" << endl;
    for (int i = 0; i < FILENUM; i++) {
        if (i % 20 == 0) {
            cout << maxTLACFileProbS2T[i] << " ";
        }
    }
    cout << endl;
    cout << "maxTLACFileProbT2S" << endl;
    for (int i = 0; i < FILENUM; i++) {
        if (i % 20 == 0) {
            cout << maxTLACFileProbT2S[i] << " ";
        }
    }
    cout << endl;
    cout << "maxTLACFileProbEqual" << endl;
    for (int i = 0; i < FILENUM; i++) {
        if (i % 20 == 0) {
            cout << maxTLACFileProbEqual[i] << " ";
        }
    }
    cout << endl;
}
/*************************************************
 * 函数名称：cacheInDistance();
 * 功能描述：计算在规定范围内，缓存命中的概率
 * 参数列表：三种时延下的命中率，算法代号：LRU0，TLAC 2（S2T），3（T2S），4（Equal）
 * 返回结果：命中概率
**************************************************/
void Graph::cacheInDistance(float& HitS2T, float& HitT2S, float& HitEqual, int filed)
{
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
                int flag = nodeArr[nei].cacheHit(f, filed);
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
                if (Distance[ser][nei] < equalTimeLim[f] && flag) {
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
            // nodeArr[i].showCacheLRU();
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
    Zipf(QF, 1);
    timelimtGenerate();
    fac();
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
    G.timeLimAwareCacheProb(); //计算TLAC的各个文件的缓存概率
    G.LRUCacheProb(); //计算LRU的各个文件的缓存概率
    for (int s = 1; s < vertex_num; s++) {
        nodeArr[s].cachefileTLAC(); //各个节点根据TLAC的概率，开始缓存文件。
        nodeArr[s].cachefileLRU(); //各个节点根据LRU的概率，开始缓存文件。
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
        cout << "\tcached file: ";
        nodeArr[i].showCacheLRU();
        // cout<<"\tTLAC S2T:";
        nodeArr[i].showCacheTLAC();
        cout << endl;
    }
    float hitS2T, hitT2S, hitEqual;
    G.cacheInDistance(hitS2T, hitT2S, hitEqual, 0);
    cout << endl
         << "Hit in Distance" << endl;
    cout << "#######################################" << endl;
    cout << hitS2T << "\t" << hitT2S << "\t" << hitEqual << endl;
    release_buff(topo, edge_num);

    return 0;
}