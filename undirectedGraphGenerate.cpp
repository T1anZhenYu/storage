#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
using namespace std;

class generageGraph {
private:
    int v_num, degree;
    vector<string> content;

public:
    generageGraph() {}
    generageGraph(int v_num, int degree)
    {
        this->v_num = v_num;
        this->degree = degree;
    }
    ~generageGraph() {}
    void work(string filename);
};
void generageGraph::work(string filename)
{
    int randnum;
    int src, des;
    vector<int> d(v_num, 0);
    vector<int> flag(v_num, 0);
    string res;
    int count;
    count = 1;
    for (int i = 0; i < v_num - degree; i++) {
        flag.assign(flag.size(), 0);
        for (int j = 0; j < degree; j++) {
            res.clear();
            if (d[i] > degree) {
                break;
            }
            d[i]++;
            src = i;
            randnum = rand() % v_num;
            while (randnum <= src || flag[randnum] == 1) {
                randnum = rand() % v_num;
            }
            des = randnum;
            flag[des] = 1;
            // cout<<"des:"<<des<<" src:"<<src<<" d[src]"<<d[src]<<" d[des]"<<d[des]<<endl;
            res.append(to_string(count) + "," + to_string(src + 1) + "," + to_string(des + 1) + ",1");
            content.push_back(res);
            res.clear();
            count++;

            if (d[des] > degree) {
                continue;
            }
            d[des]++;
            res.append(to_string(count) + "," + to_string(des + 1) + "," + to_string(src + 1) + ",1");
            content.push_back(res);
            res.clear();
            count++;
        }
    }
    cout << "content size:" << content.size() << endl;
    ofstream out("Graph.txt");
    if (out.is_open()) {
        for (int i = 0; i < content.size(); i++) {
            cout << content[i] << endl;
            out << content[i] << "\n";
        }
    }
}
int main()
{
    generageGraph* g;
    g = new generageGraph(10, 2);
    g->work("GraphContent.txt");
}