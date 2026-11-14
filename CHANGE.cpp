#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <tuple>
#include <climits>
#include <numeric>
#include <functional>
#include <filesystem>

using namespace std;
using namespace std::chrono;

// 分割字串（以空白分隔）
vector<string> split(const string &s) {
    istringstream iss(s);
    vector<string> tokens;
    string token;
    while (iss >> token)
        tokens.push_back(token);
    return tokens;
}

struct Poly {
    bool CC;
    int xp;
    int yp;
};

// 計算向量叉積，計算 (B-A) 與 (C-B) 的叉積
int cross(const pair<int, int> &A, const pair<int, int> &B, const pair<int, int> &C) {
    return (B.first - A.first) * (C.second - B.second) - (B.second - A.second) * (C.first - B.first);
}

// 讀取 CSV 檔案，分辨凹凸點
void readCSV(int argc, char* argv[],
             unordered_map<string, vector<pair<int, int>>> &cellSize,
             unordered_map<string, vector<vector<vector<Poly>>>> &cellM1Polys,
             unordered_map<string, vector<string>> &libname) {
    ifstream file;
    string line, cell, name;
    vector<pair<int, int>> onepoly;
    vector<Poly> oneCCpoly;
    vector<vector<Poly>> CCpolys;
    pair<int, int> prev, cur, next;
    int xp, yp, libnum;
    //bool pass = false;
    for (int i = 1; i < argc-2; i++) {
        file.open(argv[i]);
        name = filesystem::path(argv[i]).stem().string();
        while(getline(file, line)) {
            auto tokens = split(line);
            if(tokens[0] == "cell") {
                if (!CCpolys.empty()) {
                    cellM1Polys[cell].push_back(CCpolys); //寫入 (CCpolys非空、讀到下個cell行)
                    CCpolys.clear();
                }
                cell = tokens[1];
                libname[cell].push_back(name);
                cellSize[cell].push_back(make_pair((int)(stof(tokens[2])*1000), (int)(stof(tokens[3])*1000)));
            }
            else {
                for(int i = 1; i < tokens.size(); i+=3) {
                    xp = (int)(stof(tokens[i])*1000);
                    yp = (int)(stof(tokens[i+1])*1000);
                    onepoly.push_back({(int)(stof(tokens[i])*1000), (int)(stof(tokens[i+1])*1000)});
                }
                // 向量外積判斷凹凸點
                for(int i = 0; i < onepoly.size(); i++) {  
                    prev = onepoly[(i-1+onepoly.size())%onepoly.size()];
                    cur = onepoly[i];
                    next = onepoly[(i+1)%onepoly.size()];
                    int cp = cross(prev, cur, next);
                    if(cp > 0) // 凸點
                        oneCCpoly.push_back({true, cur.first, cur.second});
                    else if(cp < 0) // 凹點
                        oneCCpoly.push_back({false, cur.first, cur.second});
                }
                CCpolys.push_back(oneCCpoly);
                if(file.peek() == EOF) {
                    cellM1Polys[cell].push_back(CCpolys); //寫入 (最後一行)
                    CCpolys.clear();
                }
                onepoly.clear();
                oneCCpoly.clear();
            }
        }
        file.close();
    }
}

void readNextsCSV(const string &filename, unordered_map<string, unordered_map<int, vector<int>>> &libnexts){
    ifstream file(filename);
    string line, cell;
    // 讀取 CSV 檔案 (LIB 順序)
    while(getline(file, line)) {
        auto tokens = split(line);
        if(tokens[0] == "cell") {
            cell = tokens[1];
        }
        else {
            for(int i = 1; i < tokens.size(); i++){
                libnexts[cell][stoi(tokens[0])].push_back(stoi(tokens[i]));
            }
                
        }
    }
    file.close();
}

struct Rect {
    int xmin;
	int ymin;
	int xmax;
	int ymax;
};

struct Point {
    int xp;
	int yp;
    int ctype;
};

void addM1Zone(const unordered_map<string, vector<pair<int, int>>> &cellSize,
               const unordered_map<string, vector<vector<vector<Poly>>>> &cellM1Polys,
               unordered_map<string, vector<vector<Rect>>> &cellM1Zone,
               unordered_map<string, vector<vector<vector<vector<Point>>>>> &ExtendingPoint,
               unordered_map<string, vector<vector<vector<vector<Point>>>>> &MatchingPoint) {
    pair<int, int> prev, cur, next;
    int x1 = 0, y1 = 0, x2 = 0, y2 = 0, type, cellXmax = 0, cellYmax = 0;
    map<int, int> MY = { {1, 2}, {2, 1}, {3, 4}, {4, 3} };
    map<int, int> R180 = { {1, 3}, {2, 4}, {3, 1}, {4, 2} };
    map<int, int> MX = { {1, 4}, {2, 3}, {3, 2}, {4, 1} };
    vector<Rect> M1Zone;
    for(auto &cell: cellM1Polys){
        for (int libnum = 0; libnum < cell.second.size(); libnum++) {
            cellXmax = cellSize.at(cell.first)[libnum].first;
            cellYmax = cellSize.at(cell.first)[libnum].second;
            // 初始化
            ExtendingPoint[cell.first].push_back(vector<vector<vector<Point>>>(4, vector<vector<Point>>(4, vector<Point>())));
            MatchingPoint[cell.first].push_back(vector<vector<vector<Point>>>(4, vector<vector<Point>>(4, vector<Point>())));
            for(auto &poly: cell.second[libnum]){
                for(int i = 0; i < poly.size(); i++){            
                    prev = {poly[(i-1+poly.size())%poly.size()].xp, poly[(i-1+poly.size())%poly.size()].yp};
                    cur = {poly[i].xp, poly[i].yp};
                    next = {poly[(i+1)%poly.size()].xp, poly[(i+1)%poly.size()].yp};
                    if(poly[i].CC) { // 凸點
                        if(cur.first == next.first){ // 垂直線
                            y1 = min(cur.second, next.second);
                            y2 = max(cur.second, next.second);
                            if(prev.first > cur.first) { // left
                                if(y2-y1 > 36)
                                    x1 = cur.first - 9;
                                else
                                    x1 = cur.first - 16;
                                x2 = cur.first;
                            } else { // right
                                if(y2-y1 > 36)
                                    x2 = cur.first + 9;
                                else
                                    x2 = cur.first + 16;
                                x1 = cur.first;
                            }
                        }
                        else if(cur.second == next.second){ // 水平線
                            x1 = min(cur.first, next.first);
                            x2 = max(cur.first, next.first);
                            if(prev.second < cur.second) { // up
                                if(x2-x1 > 36)
                                    y2 = cur.second + 9;
                                else
                                    y2 = cur.second + 16;
                                y1 = cur.second;
                            } else { // down
                                if(x2-x1 > 36)
                                    y1 = cur.second - 9;
                                else
                                    y1 = cur.second - 16;
                                y2 = cur.second;
                            }
                        }
                    }
                    else { // 凹點
                        if(cur.first == next.first){ // 垂直線
                            y1 = min(cur.second, next.second);
                            y2 = max(cur.second, next.second);
                            if(prev.first > cur.first) { // left
                                if(y2-y1 > 36)
                                    x2 = cur.first + 9;
                                else
                                    x2 = cur.first + 16;
                                x1 = cur.first;
                            } else { // right
                                if(y2-y1 > 36)
                                    x1 = cur.first - 9;
                                else
                                    x1 = cur.first - 16;
                                x2 = cur.first;
                            }
                        }
                        else if(cur.second == next.second){ // 水平線
                            x1 = min(cur.first, next.first);
                            x2 = max(cur.first, next.first);
                            if(prev.second < cur.second) { // up
                                if(x2-x1 > 36)
                                    y1 = cur.second - 9;
                                else
                                    y1 = cur.second - 16;
                                y2 = cur.second;
                            } else { // down
                                if(x2-x1 > 36)
                                    y2 = cur.second + 9;
                                else
                                    y2 = cur.second + 16;
                                y1 = cur.second;
                            }
                        }
                    }
                    M1Zone.push_back({x1, y1, x2, y2});
                    if(poly[i].CC) { // 凸點
                        int flipX = cellXmax - cur.first, flipY = cellYmax - cur.second;
                        if(cur.first < next.first && cur.second < prev.second) type = 1; // 指向左下角
                        else if(cur.second < next.second && cur.first > prev.first) type = 2; // 指向右下角
                        else if(cur.first > next.first && cur.second > prev.second) type = 3; // 指向右上角
                        else if(cur.second > next.second && cur.first < prev.first) type = 4; // 指向左上角
                        if(cur.first < 47) { // Point 左
                            MatchingPoint[cell.first][libnum][0][0].push_back({cur.first, cur.second, type}); //不翻轉 Point 左
                            MatchingPoint[cell.first][libnum][1][1].push_back({flipX, cur.second, MY[type]}); //MY翻轉後 Point 右
                            MatchingPoint[cell.first][libnum][2][1].push_back({flipX, flipY, R180[type]}); //R180翻轉後 Point 右
                            MatchingPoint[cell.first][libnum][3][0].push_back({cur.first, flipY, MX[type]}); //MX翻轉後 Point 左
                            if(cur.first < 9) {
                                ExtendingPoint[cell.first][libnum][0][0].push_back({cur.first, cur.second, type});
                                ExtendingPoint[cell.first][libnum][1][1].push_back({flipX, cur.second, MY[type]});
                                ExtendingPoint[cell.first][libnum][2][1].push_back({flipX, flipY, R180[type]});
                                ExtendingPoint[cell.first][libnum][3][0].push_back({cur.first, flipY, MX[type]});
                            }
                        }
                        if(cur.first > cellXmax - 47) { // Point 右
                            MatchingPoint[cell.first][libnum][0][1].push_back({cur.first, cur.second, type}); //不翻轉 Point 右
                            MatchingPoint[cell.first][libnum][1][0].push_back({flipX, cur.second, MY[type]}); //MY翻轉後 Point 左
                            MatchingPoint[cell.first][libnum][2][0].push_back({flipX, flipY, R180[type]}); //R180翻轉後 Point 左
                            MatchingPoint[cell.first][libnum][3][1].push_back({cur.first, flipY, MX[type]}); //MX翻轉後 Point 右
                            if(cur.first > cellXmax - 9) {
                                ExtendingPoint[cell.first][libnum][0][1].push_back({cur.first, cur.second, type});
                                ExtendingPoint[cell.first][libnum][1][0].push_back({flipX, cur.second, MY[type]});
                                ExtendingPoint[cell.first][libnum][2][0].push_back({flipX, flipY, R180[type]});
                                ExtendingPoint[cell.first][libnum][3][1].push_back({cur.first, flipY, MX[type]});
                            }
                        }
                        if(cur.second > cellYmax - 47) { // Point 上
                            MatchingPoint[cell.first][libnum][0][2].push_back({cur.first, cur.second, type}); //不翻轉 Point 上
                            MatchingPoint[cell.first][libnum][1][2].push_back({flipX, cur.second, MY[type]}); //MY翻轉後 Point 上
                            MatchingPoint[cell.first][libnum][2][3].push_back({flipX, flipY, R180[type]}); //R180翻轉後 Point 下
                            MatchingPoint[cell.first][libnum][3][3].push_back({cur.first, flipY, MX[type]}); //MX翻轉後 Point 下
                            if(cur.second > cellYmax - 9) {
                                ExtendingPoint[cell.first][libnum][0][2].push_back({cur.first, cur.second, type});
                                ExtendingPoint[cell.first][libnum][1][2].push_back({flipX, cur.second, MY[type]});
                                ExtendingPoint[cell.first][libnum][2][3].push_back({flipX, flipY, R180[type]});
                                ExtendingPoint[cell.first][libnum][3][3].push_back({cur.first, flipY, MX[type]});
                            }
                        }
                        if(cur.second < 47) { // Point 下
                            MatchingPoint[cell.first][libnum][0][3].push_back({cur.first, cur.second, type}); //不翻轉 Point 下
                            MatchingPoint[cell.first][libnum][1][3].push_back({flipX, cur.second, MY[type]}); //MY翻轉後 Point 下
                            MatchingPoint[cell.first][libnum][2][2].push_back({flipX, flipY, R180[type]}); //R180翻轉後 Point 上
                            MatchingPoint[cell.first][libnum][3][2].push_back({cur.first, flipY, MX[type]}); //MX翻轉後 Point 上
                            if(cur.second < 9) {
                                ExtendingPoint[cell.first][libnum][0][3].push_back({cur.first, cur.second, type});
                                ExtendingPoint[cell.first][libnum][1][3].push_back({flipX, cur.second, MY[type]});
                                ExtendingPoint[cell.first][libnum][2][2].push_back({flipX, flipY, R180[type]});
                                ExtendingPoint[cell.first][libnum][3][2].push_back({cur.first, flipY, MX[type]});
                            }
                        }
                    }
                }
            }
            cellM1Zone[cell.first].push_back(M1Zone);
            M1Zone.clear();
        }
		// cout << "HAHAHAHAHAHAHAHAHAH!!!!!" << endl;
    }
}

struct Extend {
    int lmin;
	int lmax;
	int lout;
};

struct Event {
    int pos, weight;
};

vector<Extend> mergeMaxIntervals(const vector<Extend>& ivs) {
    vector<Event> evs;
    evs.reserve(ivs.size()*2);
    bool isPositive = true;
    for (auto& iv : ivs) {
        if(iv.lout > 0){
            evs.push_back({iv.lmin,   +iv.lout}); // 添加權重
            evs.push_back({iv.lmax+1, -iv.lout}); // 刪除權重
        }
        else{  // 權重為負數取絕對值 (lout < 0)
            isPositive = false;
            evs.push_back({iv.lmin,   -iv.lout}); // 刪除權重
            evs.push_back({iv.lmax+1, +iv.lout}); // 刪除權重
        }
    }
    sort(evs.begin(), evs.end(), [](auto &a, auto &b){
        if (a.pos != b.pos) return a.pos < b.pos; // 升序排列
        return a.weight > b.weight;  // 先 +w 再 -w
    });

    multiset<int> ms; // 存取目前活動區間的所有權重
    vector<Extend> res;
    int prevPos = 0, prevMax = 0;
    size_t i = 0, n = evs.size();
    while (i < n) {
        int curPos = evs[i].pos;
        // 先處裡相同位置的事件
        while (i < n && evs[i].pos == curPos) {
            if (evs[i].weight > 0) ms.insert(evs[i].weight);
            else ms.erase(ms.find(-evs[i].weight));
            ++i;
        }
        int newMax = ms.empty()?0:*ms.rbegin(); // 獲取目前最大權重(若為空則為0)
        if (newMax != prevMax) {
            if (prevMax > 0 && isPositive)
                res.push_back({prevPos, curPos-1, prevMax});
            else if (prevMax > 0 && !isPositive)
                res.push_back({prevPos, curPos-1, -prevMax});
            prevPos = curPos;
            prevMax = newMax;
        }
    }
    return res;
}

// 填充區段
void fillSeg(int i, int j, int l, int r, int v, vector<vector<vector<int>>>& Matching) {
    if (j == 0 || j == 3){ // 左、下
        for (int k = l; k <= r; k++)
            if (k >= 0 && k < Matching[i][j].size()) 
                Matching[i][j][k] = min(Matching[i][j][k], v);
    }
    else{ // 右、上
        for (int k = l; k <= r; k++)
            if (k >= 0 && k < Matching[i][j].size()) 
                Matching[i][j][k] = max(Matching[i][j][k], v);
    } 
}

struct SparseTable {
    int N, K;
    vector<int> lg2;
    vector<vector<int>> st;
    function<int(int, int)> op;  // 比較操作：最大值 or 最小值

    // 建構：建立陣列 a，0-based，下標 [0…N-1]
    SparseTable(const vector<int>& a, function<int(int, int)> cmp) {
        N = a.size();
        K = floor(log2(N)) + 1;
        op = cmp;
         // 計算 N 取 log2 後無條件捨棄到整數的數值 
        lg2.resize(N + 1);
        lg2[1] = 0;
        for (int i = 2; i <= N; i++)
            lg2[i] = lg2[i/2] + 1;
        st.assign(K, vector<int>(N));
        // k = 0：長度 1 區間
        for (int i = 0; i < N; i++)
            st[0][i] = a[i];
        // k >= 1：合併兩個長度 2^{k-1} 的子區間
        for (int k = 1; k < K; k++) {
            int len = 1 << (k - 1);
            for (int i = 0; i + (1 << k) <= N; i++) {
                st[k][i] = op(st[k-1][i], st[k-1][i + len]);
            }
        }
    }

    // 查詢 [L, R] 的最大值，從0開始
    int query(int L, int R) {
        int len = R - L + 1;
        int k = lg2[len];
        return op(st[k][L], st[k][R - (1<<k) + 1]);
    }
};

// 執行 : 計算Extending cell 凸出的位置(left/right/up/down、xmin、xmax)、cell 四周 Matcing 矩陣
void Scan(const unordered_map<string, vector<pair<int, int>>> &cellSize,
          const unordered_map<string, vector<vector<Rect>>> &cellM1Zone, 
          unordered_map<string, bool> &isExtending,
          unordered_map<string, vector<vector<vector<vector<Extend>>>>> &ExtendingLength,
          unordered_map<string, vector<vector<vector<SparseTable>>>> &MatchingLength){
    // 翻轉順序: [0]為R0(N)、[1]為MY(FN)、[2]為R180(S)、[3]為MX(FS)
    int cellXmin = 0, cellXmax = 0, cellYmin = 0, cellYmax = 0;
    for(auto &cell: cellM1Zone){ //矩形"邊"
        // 初始化
        isExtending[cell.first] = false;
        for(int libnum = 0; libnum < cell.second.size(); libnum++){
            cellXmax = cellSize.at(cell.first)[libnum].first;
            cellYmax = cellSize.at(cell.first)[libnum].second;
            // 初始化    
            ExtendingLength[cell.first].push_back(vector<vector<vector<Extend>>>(4, vector<vector<Extend>>(4, vector<Extend>())));
            MatchingLength[cell.first].push_back(vector<vector<SparseTable>>(4, vector<SparseTable>(4, SparseTable({0}, [](int a, int b){return max(a, b);}))));
            vector<vector<vector<int>>> Matching(4, vector<vector<int>>(4, vector<int>()));
            for(int i = 0; i < 4; i++){
                Matching[i][0].resize(cellYmax+1, cellXmax); //left
                Matching[i][1].resize(cellYmax+1, 0); //right
                Matching[i][2].resize(cellXmax+1, 0); //up
                Matching[i][3].resize(cellXmax+1, cellYmax); //down
            }
            for(auto &r: cell.second[libnum]){
                int flipXmin = cellXmax - r.xmax, flipXmax = cellXmax - r.xmin, flipYmin = cellYmax - r.ymax, flipYmax = cellYmax - r.ymin;
                
                // 填充區段(Matching)(不紀錄四個角落)
                fillSeg(0, 0, r.ymin, r.ymax, r.xmin, Matching); // R0 左
                fillSeg(0, 1, r.ymin, r.ymax, r.xmax, Matching); // R0 右
                fillSeg(0, 2, r.xmin, r.xmax, r.ymax, Matching); // R0 上
                fillSeg(0, 3, r.xmin, r.xmax, r.ymin, Matching); // R0 下
                fillSeg(1, 0, r.ymin, r.ymax, flipXmin, Matching); // MY 左
                fillSeg(1, 1, r.ymin, r.ymax, flipXmax, Matching); // MY 右
                fillSeg(1, 2, flipXmin, flipXmax, r.ymax, Matching); // MY 上
                fillSeg(1, 3, flipXmin, flipXmax, r.ymin, Matching); // MY 下
                fillSeg(2, 0, flipYmin, flipYmax, flipXmin, Matching); // R180 左
                fillSeg(2, 1, flipYmin, flipYmax, flipXmax, Matching); // R180 右
                fillSeg(2, 2, flipXmin, flipXmax, flipYmax, Matching); // R180 上
                fillSeg(2, 3, flipXmin, flipXmax, flipYmin, Matching); // R180 下
                fillSeg(3, 0, flipYmin, flipYmax, r.xmin, Matching); // MX 左
                fillSeg(3, 1, flipYmin, flipYmax, r.xmax, Matching); // MX 右
                fillSeg(3, 2, r.xmin, r.xmax, flipYmax, Matching); // MX 上
                fillSeg(3, 3, r.xmin, r.xmax, flipYmin, Matching); // MX 下

                if(r.xmin < cellXmin){  //凸出左方
                    isExtending[cell.first] = true;
                    ExtendingLength[cell.first][libnum][0][0].push_back({r.ymin, r.ymax, r.xmin}); //不翻轉凸出左方
                    ExtendingLength[cell.first][libnum][1][1].push_back({r.ymin, r.ymax, flipXmax}); //MY翻轉後凸出右方
                    ExtendingLength[cell.first][libnum][2][1].push_back({flipYmin, flipYmax, flipXmax}); //R180翻轉後凸出右方
                    ExtendingLength[cell.first][libnum][3][0].push_back({flipYmin, flipYmax, r.xmin}); //MX翻轉後凸出左方
                }
                if(r.xmax > cellXmax){  //凸出右方
                    isExtending[cell.first] = true;
                    ExtendingLength[cell.first][libnum][0][1].push_back({r.ymin, r.ymax, r.xmax}); //不翻轉凸出右方
                    ExtendingLength[cell.first][libnum][1][0].push_back({r.ymin, r.ymax, flipXmin}); //MY翻轉後凸出左方
                    ExtendingLength[cell.first][libnum][2][0].push_back({flipYmin, flipYmax, flipXmin}); //R180翻轉後凸出左方
                    ExtendingLength[cell.first][libnum][3][1].push_back({flipYmin, flipYmax, r.xmax}); //MX翻轉後凸出右方
                }
                if(r.ymax > cellYmax){  //凸出上方
                    isExtending[cell.first] = true;
                    ExtendingLength[cell.first][libnum][0][2].push_back({r.xmin, r.xmax, r.ymax}); //不翻轉凸出上方
                    ExtendingLength[cell.first][libnum][1][2].push_back({flipXmin, flipXmax, r.ymax}); //MY翻轉後凸出上方
                    ExtendingLength[cell.first][libnum][2][3].push_back({flipXmin, flipXmax, flipYmin}); //R180翻轉後凸出下方
                    ExtendingLength[cell.first][libnum][3][3].push_back({r.xmin, r.xmax, flipYmin}); //MX翻轉後凸出下方
                }
                if(r.ymin < cellYmin){  //凸出下方
                    isExtending[cell.first] = true;
                    ExtendingLength[cell.first][libnum][0][3].push_back({r.xmin, r.xmax, r.ymin}); //不翻轉凸出下方
                    ExtendingLength[cell.first][libnum][1][3].push_back({flipXmin, flipXmax, r.ymin}); //MY翻轉後凸出下方
                    ExtendingLength[cell.first][libnum][2][2].push_back({flipXmin, flipXmax, flipYmax}); //R180翻轉後凸出上方
                    ExtendingLength[cell.first][libnum][3][2].push_back({r.xmin, r.xmax, flipYmax}); //MX翻轉後凸出上方
                }
            }
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    if (j == 0 || j == 3){ // 左、下(取 min)
                        MatchingLength[cell.first][libnum][i][j] = SparseTable(Matching[i][j], [](int x, int y) { return min(x, y); });}
                    else{ // 右、上(取 max)
                        MatchingLength[cell.first][libnum][i][j] = SparseTable(Matching[i][j], [](int x, int y) { return max(x, y); });}
                    if(!ExtendingLength[cell.first][libnum][i][j].empty()) // 重疊區間取(正)最大、(負)最小值
                        ExtendingLength[cell.first][libnum][i][j] = mergeMaxIntervals(ExtendingLength[cell.first][libnum][i][j]);
                }   
            }
        }
    }
}

struct DEFfile {
    string instance;
    string cell;
    int xpoint;
    int ypoint;
	string type;
};

// 比較函數：先比較 a，不相同則依 a 排序；若 a 相同則比較 b
bool DEFcmp(const DEFfile &d1, const DEFfile &d2) {
    if (d1.ypoint != d2.ypoint)
        return d1.ypoint < d2.ypoint;
    return d1.xpoint < d2.xpoint;
}

// 讀取 DEF 檔案
void readDEF(const string &filename, vector<DEFfile> &DEFdata) {
    ifstream file(filename);
    string line, cell;
    bool start = false;
    // 讀取 DEF 檔案
    while(getline(file, line)) {
        if(line.substr(0, 10) == "COMPONENTS") start = true; // 開始讀取 COMPONENTS 區段
        else if(line == "END COMPONENTS") break; // 結束讀取
        else if(start && line.substr(0, 2) == " -") {
            auto tokens = split(line);            
            DEFdata.push_back({tokens[1], tokens[2], stoi(tokens[6]), stoi(tokens[7]), tokens[9]});
        };
    }
    file.close();
    sort(DEFdata.begin(), DEFdata.end(), DEFcmp); // 利用 sort 進行排序
}

struct ChipArray {
    string instance;
    string cell;
    int xstart;
    int ystart;
    int xend;
    int libnum;
    int type;
};

// 定義 pair<int, int> 的 hash 函數
struct pair_hash {
    size_t operator()(const pair<int, int>& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

// 製作 chip 表
void createTable(const unordered_map<string, vector<pair<int, int>>> &cellSize,
                 vector<DEFfile> &DEFdata, 
                 vector<vector<ChipArray>> &chip, 
                 unordered_map<string, unordered_set<pair<int,int>, pair_hash>> &cellIndex) {
    vector<ChipArray> temp;
    int yindex = 0;
    int xindex = 0;
    int ytemp = 0;
    int type;
    int setEO;
    int even = 0, odd = 0;
    for (const auto& d : DEFdata) {
        if(ytemp== 0) ytemp = d.ypoint; 
        else if(d.ypoint == ytemp) xindex++; // 同一行
        else { // 進位
            chip.push_back(temp);
            temp.clear();
            if (d.ypoint > ytemp + cellSize.at(d.cell)[0].second){  //空白歷史(跨整行) (儲存 ORI cell)
                yindex++;
                chip.push_back(temp);
            }
            ytemp = d.ypoint;
            yindex++;
            xindex = 0;
        }
        if (d.type == "N") type = 0;
        else if (d.type == "FN") type = 1;
        else if (d.type == "S") type = 2;
        else if (d.type == "FS") type = 3;

        if (cellSize.at(d.cell)[0].first % 36 == 18 && (type == 1 || type == 2)){ //FN(MY) or S(R180) even/odd reverse
            if ((d.xpoint-2500) % 36 == 0){ //even -> odd
                setEO = -1;
                odd++;
            }
            if((d.xpoint-2500) % 36 == 18){ //odd -> even
                setEO = -2;
                even++;
            }
        }
        else if ((d.xpoint-2500) % 36 == 0){ //even
            setEO = -2;
            even++;
        }
        else if((d.xpoint-2500) % 36 == 18){ //odd
            setEO = -1;
            odd++;
        }
        else
            cout << "error" << endl;
        temp.push_back({d.instance, d.cell, d.xpoint, d.ypoint, d.xpoint + cellSize.at(d.cell)[0].first, setEO, type});
        cellIndex[d.cell].insert({xindex, yindex});
    }
    cout << "even: " << even << " odd: " << odd << endl;
    chip.push_back(temp); //最後一行加上
}

//struct Change {
//    int libnum;
    //int type;
    //bool rotate;
//};

// 取得元素，若元素存在則傳回 true，不存在則傳回 false
bool getElement(const vector<vector<ChipArray>> &matrix, int Xindex, int Yindex) {
    if (Yindex < matrix.size() && Xindex < matrix[Yindex].size())
        return true;
    return false;
}

// 二分搜尋法，找出 xs 到 xe 之間的 cell index
pair<int, int> binarySearch(const vector<ChipArray> &matrix, int xs, int xe) {
    int l = 0, r = matrix.size() - 1;
    int xa = -1, xb = -1;
    while(l <= r){
        if (l == r) {
            if (matrix[l].xend >= xs) xa = r; // 找到起始cell
            break;
        }
        int m = (l + r) / 2;
        if (matrix[m].xend >= xs) r = m; // 縮小至前半段 (保留右端點)
        else l = m + 1; // 縮小至後半段
    }
    l = 0;
    r = matrix.size() - 1;
    while(l <= r){
        if (l == r) {
            if (matrix[l].xstart <= xe) xb = l; // 找到終止cell
            break;
        }
        int m = (l + r - 1) / 2 + 1; // 向上取整數
        if (matrix[m].xstart <= xe) l = m; // 縮小至後半段 (保留左端點)
        else r = m - 1; // 縮小至前半段
    }
    return make_pair(xa, xb);
}

void Judgment(const unordered_map<string, vector<pair<int, int>>> &cellSize,
              const unordered_map<string, unordered_map<int, vector<int>>> &libnexts,
              const unordered_map<string, bool> &isExtending,
              const unordered_map<string, vector<vector<vector<vector<Extend>>>>> &ExtendingLength,
              unordered_map<string, vector<vector<vector<SparseTable>>>> &MatchingLength,
              const unordered_map<string, vector<vector<vector<vector<Point>>>>> &ExtendingPoint,
              const unordered_map<string, vector<vector<vector<vector<Point>>>>> &MatchingPoint,
              unordered_map<pair<int,int>, int, pair_hash> &tempChange,
              const int Xindex, const int Yindex, bool &canChange,
              const vector<vector<ChipArray>> &chip) {
    int subXindex, subYindex, run;
    pair<int, int> subXindexs;
    //unordered_map<int, int> rotateType = {{0, 1}, {1, 0}, {2, 3}, {3, 2}}; //翻轉型態(左右翻轉)  U66311
    ChipArray target, subtarget;
    //vector<int> sublibnums;
    bool isContact = false;
    int libnum, sublibnum;
    target = chip[Yindex][Xindex];
    //cout << Xindex << " " << Yindex << endl;
    unordered_map<pair<int,int>, int, pair_hash> constempChange = tempChange; // 固定候選
    if (tempChange.find({Xindex, Yindex}) != tempChange.end()) libnum = tempChange.at({Xindex, Yindex}); // target的libnum有紀錄接著用
    else libnum = target.libnum;
    //if(target.libnum == -1){ // total未記錄   totalchange.find({Xindex, Yindex}) == totalchange.end()
    for(int ln : libnexts.at(target.cell).at(libnum)) {
        //if(ln == 0) break; // 只對extending的cell做判斷
        target.xend = target.xstart + cellSize.at(target.cell)[ln].first; // 更新 xend
        canChange = true; 
        tempChange[{Xindex, Yindex}] = ln; // 紀錄候選
        if (!ExtendingLength.at(target.cell)[ln][target.type][0].empty() && canChange){  // 左凸
            //cout << "left" << endl;
            subXindex = Xindex - 1;
            subYindex = Yindex;
            if (getElement(chip, subXindex, subYindex)) { // 左邊有cell
                subtarget = chip[subYindex][subXindex];
                if (tempChange.find({subXindex, subYindex}) != tempChange.end()) sublibnum = tempChange.at({subXindex, subYindex}); // subtarget的libnum有紀錄接著用
                else sublibnum = subtarget.libnum;
                for(int sln : libnexts.at(subtarget.cell).at(sublibnum)) {
                    subtarget.xend = subtarget.xstart + cellSize.at(subtarget.cell)[sln].first; // 更新 xend
                    canChange = true;
                    if (subtarget.xend == target.xstart){ // 左邊cell與本cell接觸
                        int add = target.xstart - subtarget.xstart;
                        // tip、side 檢查
                        for(auto &Ext: ExtendingLength.at(target.cell)[ln][target.type][0]) { // 本 cell 的左 Extending    
                            // 左邊 cell 的右 Matching 與本 cell 的左 Extending 相交判斷
                            if(!(MatchingLength.at(subtarget.cell)[sln][subtarget.type][1].query(Ext.lmin, Ext.lmax) <= (Ext.lout + add))){
                                canChange = false;
                                break;
                            }    
                        }
                        //corner 檢查
                        if(canChange) {
                            for(auto &EP: ExtendingPoint.at(target.cell)[ln][target.type][0]) {
                                for(auto &MP: MatchingPoint.at(subtarget.cell)[sln][subtarget.type][1]) {
                                    int xdist = EP.xp + add - MP.xp, ydist = EP.yp - MP.yp;
                                    if (hypot(xdist, ydist) < 20) {
                                        // 指向左下、右上產生 corner
                                        if(((EP.ctype == 1 && MP.ctype == 3) || (EP.ctype == 3 && MP.ctype == 1)) && (xdist * ydist > 0)) canChange = false;
                                        // 指向右下、左上產生 corner
                                        if(((EP.ctype == 2 && MP.ctype == 4) || (EP.ctype == 4 && MP.ctype == 2)) && (xdist * ydist < 0)) canChange = false;
                                    }
                                    if (!canChange) break;
                                }
                                if (!canChange) break;
                            }
                        }
                        if(canChange) { 
                            if (isExtending.at(subtarget.cell) && tempChange.find({subXindex, subYindex}) == tempChange.end()) { // 左邊cell也有凸出且未記錄
                                tempChange[{subXindex, subYindex}] = sln; // 紀錄候選
                                Judgment(cellSize, libnexts, isExtending, ExtendingLength, MatchingLength, ExtendingPoint, MatchingPoint, tempChange, subXindex, subYindex, canChange, chip); // 針對該cell再判斷
                            }
                            else tempChange[{subXindex, subYindex}] = sln; // 紀錄候選
                            if (canChange) break;
                            else tempChange.erase({subXindex, subYindex}); // 清除候選(不可改)
                        }
                    }
                }
            }
        }

        if (!ExtendingLength.at(target.cell)[ln][target.type][1].empty() && canChange){  // 右凸
            //cout << "right" << endl;
            subXindex = Xindex + 1;
            subYindex = Yindex;
            if (getElement(chip, subXindex, subYindex)) { // 右邊有cell
                subtarget = chip[subYindex][subXindex];
                if (tempChange.find({subXindex, subYindex}) != tempChange.end()) sublibnum = tempChange.at({subXindex, subYindex}); // subtarget的libnum有紀錄接著用
                else sublibnum = subtarget.libnum;
                for(int sln : libnexts.at(subtarget.cell).at(sublibnum)) {
                    subtarget.xend = subtarget.xstart + cellSize.at(subtarget.cell)[sln].first; // 更新 xend
                    canChange = true;
                    if (subtarget.xstart == target.xend){ // 右邊cell與本cell接觸
                        int add = target.xstart - subtarget.xstart;
                        // tip、side 檢查
                        for(auto &Ext: ExtendingLength.at(target.cell)[ln][target.type][1]) { // 本 cell 的右 Extending    
                            // 右邊 cell 的左 Matching 與本 cell 的右 Extending 相交判斷
                            if(!(MatchingLength.at(subtarget.cell)[sln][subtarget.type][0].query(Ext.lmin, Ext.lmax) >= (Ext.lout + add))){
                                canChange = false;
                                break;
                            }    
                        }
                        //corner 檢查
                        if(canChange) {
                            for(auto &EP: ExtendingPoint.at(target.cell)[ln][target.type][1]) {
                                for(auto &MP: MatchingPoint.at(subtarget.cell)[sln][subtarget.type][0]) {
                                    int xdist = EP.xp + add - MP.xp, ydist = EP.yp - MP.yp;
                                    if (hypot(xdist, ydist) < 20) {
                                        // 指向左下、右上產生 corner
                                        if(((EP.ctype == 1 && MP.ctype == 3) || (EP.ctype == 3 && MP.ctype == 1)) && (xdist * ydist > 0)) canChange = false;
                                        // 指向右下、左上產生 corner
                                        if(((EP.ctype == 2 && MP.ctype == 4) || (EP.ctype == 4 && MP.ctype == 2)) && (xdist * ydist < 0)) canChange = false;
                                    }
                                    if (!canChange) break;
                                }
                                if (!canChange) break;
                            }
                        }
                        if(canChange) {
                            if (isExtending.at(subtarget.cell) && tempChange.find({subXindex, subYindex}) == tempChange.end()) { // 右邊cell也有凸出且未記錄
                                tempChange[{subXindex, subYindex}] = sln; // 紀錄候選
                                Judgment(cellSize, libnexts, isExtending, ExtendingLength, MatchingLength, ExtendingPoint, MatchingPoint, tempChange, subXindex, subYindex, canChange, chip); // 針對該cell再判斷
                            } 
                            else tempChange[{subXindex, subYindex}] = sln; // 紀錄候選
                            if (canChange) break;
                            else tempChange.erase({subXindex, subYindex}); // 清除候選(不可改)
                        }
                    }  
                }
            }
        }

        if (!ExtendingLength.at(target.cell)[ln][target.type][2].empty() && canChange){  //上凸
            //cout << "up" << endl;
            subYindex = Yindex + 1;
            if (getElement(chip, 0, subYindex)) { // 上邊有整行
                subXindexs = binarySearch(chip[subYindex], target.xstart, target.xend);
                if (subXindexs.first != -1 && subXindexs.second != -1 && subXindexs.first <= subXindexs.second) { // 上邊有cell
                    for(int i = subXindexs.first; i <= subXindexs.second; i++){
                        subtarget = chip[subYindex][i];
                        if (tempChange.find({i, subYindex}) != tempChange.end()) sublibnum = tempChange.at({i, subYindex}); // subtarget的libnum有紀錄接著用
                        else sublibnum = subtarget.libnum;
                        for(int sln : libnexts.at(subtarget.cell).at(sublibnum)) {
                            subtarget.xend = subtarget.xstart + cellSize.at(subtarget.cell)[sln].first; // 更新 xend
                            canChange = true;
                            int Xadd = target.xstart - subtarget.xstart;
                            int Yadd = target.ystart - subtarget.ystart;
                            // tip、side 檢查
                            for(auto &Ext: ExtendingLength.at(target.cell)[ln][target.type][2]) { // 本 cell 的上 Extending
                                if(!((subtarget.xstart < target.xstart + Ext.lmin && subtarget.xend < target.xstart + Ext.lmin) ||
                                    (subtarget.xstart > target.xstart + Ext.lmax && subtarget.xend > target.xstart + Ext.lmax))) { // 接觸本 cell 的 Extending 範圍 (排除法)
                                    isContact = true;
                                    int fitlmin = max(0,Ext.lmin + Xadd);
                                    int fitlmax = min(cellSize.at(subtarget.cell)[sln].first, Ext.lmax + Xadd);
                                    // 上邊 cell 的下 Matching 與本 cell 的上 Extending 相交判斷
                                    if(!(MatchingLength.at(subtarget.cell)[sln][subtarget.type][3].query(fitlmin, fitlmax) >= (Ext.lout + Yadd))){
                                        canChange = false;
                                        break;
                                    }
                                }
                            }
                            //corner 檢查
                            if(canChange && isContact) { // 可改且有接觸到本 cell 的 Extending 範圍
                                for(auto &EP: ExtendingPoint.at(target.cell)[ln][target.type][2]) {
                                    for(auto &MP: MatchingPoint.at(subtarget.cell)[sln][subtarget.type][3]) {
                                        int xdist = EP.xp + Xadd - MP.xp, ydist = EP.yp + Yadd - MP.yp;
                                        if (hypot(xdist, ydist) < 20) {
                                            // 指向左下、右上產生 corner
                                            if(((EP.ctype == 1 && MP.ctype == 3) || (EP.ctype == 3 && MP.ctype == 1)) && (xdist * ydist > 0)) canChange = false;
                                            // 指向右下、左上產生 corner
                                            if(((EP.ctype == 2 && MP.ctype == 4) || (EP.ctype == 4 && MP.ctype == 2)) && (xdist * ydist < 0)) canChange = false;
                                        }
                                        if (!canChange) break;
                                    }
                                    if (!canChange) break;
                                }
                            }
                            if(canChange && isContact) { // 可改且有接觸到本 cell 的 Extending 範圍
                                if (isExtending.at(subtarget.cell) && tempChange.find({i, subYindex}) == tempChange.end()) { // 上邊cell也有凸出且未記錄
                                    tempChange[{i, subYindex}] = sln; // 紀錄候選
                                    Judgment(cellSize, libnexts, isExtending, ExtendingLength, MatchingLength, ExtendingPoint, MatchingPoint, tempChange, i, subYindex, canChange, chip); // 針對該cell再判斷
                                    if (!canChange) tempChange.erase({i, subYindex}); // 清除候選(不可改)
                                }
                                else tempChange[{i, subYindex}] = sln; // 紀錄候選
                            }
                            isContact = false;
                            if(canChange) break;
                        }
                        if(!canChange) break;
                    }
                }
            }
        }

        if (!ExtendingLength.at(target.cell)[ln][target.type][3].empty() && canChange){  //下凸
            //cout << "down" << endl;
            subYindex = Yindex - 1;
            if (getElement(chip, 0, subYindex)) { // 下邊有整行
                subXindexs = binarySearch(chip[subYindex], target.xstart, target.xend);
                if (subXindexs.first != -1 && subXindexs.second != -1 && subXindexs.first <= subXindexs.second) { // 下邊有cell
                    for(int i = subXindexs.first; i <= subXindexs.second; i++){
                        subtarget = chip[subYindex][i];
                        if (tempChange.find({i, subYindex}) != tempChange.end()) sublibnum = tempChange.at({i, subYindex}); // subtarget的libnum有紀錄接著用
                        else sublibnum = subtarget.libnum;
                        for(int sln : libnexts.at(subtarget.cell).at(sublibnum)) {
                            subtarget.xend = subtarget.xstart + cellSize.at(subtarget.cell)[sln].first; // 更新 xend
                            canChange = true;
                            int Xadd = target.xstart - subtarget.xstart;
                            int Yadd = target.ystart - subtarget.ystart;
                            // tip、side 檢查
                            for(auto &Ext: ExtendingLength.at(target.cell)[ln][target.type][3]) { // 本 cell 的下 Extending
                                if(!((subtarget.xstart < target.xstart + Ext.lmin && subtarget.xend < target.xstart + Ext.lmin) ||
                                    (subtarget.xstart > target.xstart + Ext.lmax && subtarget.xend > target.xstart + Ext.lmax))) { // 接觸本 cell 的 Extending 範圍 (排除法)
                                    isContact = true;                
                                    int fitlmin = max(0,Ext.lmin + Xadd);
                                    int fitlmax = min(cellSize.at(subtarget.cell)[sln].first, Ext.lmax + Xadd);
                                    // 下邊 cell 的上 Matching 與本 cell 的下 Extending 相交判斷
                                    if(!(MatchingLength.at(subtarget.cell)[sln][subtarget.type][2].query(fitlmin, fitlmax) <= (Ext.lout + Yadd))){
                                        canChange = false;
                                        break;
                                    }
                                }
                            }
                            //corner 檢查
                            if(canChange && isContact) { // 可改且有接觸到本 cell 的 Extending 範圍
                                for(auto &EP: ExtendingPoint.at(target.cell)[ln][target.type][3]) {
                                    for(auto &MP: MatchingPoint.at(subtarget.cell)[sln][subtarget.type][2]) {
                                        int xdist = EP.xp + Xadd - MP.xp, ydist = EP.yp + Yadd - MP.yp;
                                        if (hypot(xdist, ydist) < 20) {
                                            // 指向左下、右上產生 corner
                                            if(((EP.ctype == 1 && MP.ctype == 3) || (EP.ctype == 3 && MP.ctype == 1)) && (xdist * ydist > 0)) canChange = false;
                                            // 指向右下、左上產生 corner
                                            if(((EP.ctype == 2 && MP.ctype == 4) || (EP.ctype == 4 && MP.ctype == 2)) && (xdist * ydist < 0)) canChange = false;
                                        }
                                        if (!canChange) break;
                                    }
                                    if (!canChange) break;
                                }
                            }
                            if(canChange && isContact) { // 可改且有接觸到本 cell 的 Extending 範圍
                                if (isExtending.at(subtarget.cell) && tempChange.find({i, subYindex}) == tempChange.end()) {// 下邊cell也有凸出且未記錄
                                    tempChange[{i, subYindex}] = sln; // 紀錄候選
                                    Judgment(cellSize, libnexts, isExtending, ExtendingLength, MatchingLength, ExtendingPoint, MatchingPoint, tempChange, i, subYindex, canChange, chip); // 針對該cell再判斷
                                    if (!canChange) tempChange.erase({i, subYindex}); // 清除候選(不可改)
                                }
                                else tempChange[{i, subYindex}] = sln; // 紀錄候選
                                
                            }
                            isContact = false;
                            if(canChange) break;
                        }
                        if(!canChange) break;
                    }
                }
            }
        }
        /*
        run = 0;
        while(run == 0 || (run < 2 && !canChange)){
            run++;
            
            if (run == 2) {
                target.type = rotateType[target.type]; // 嘗試翻轉
                tempChange[{Xindex, Yindex}].rotate = true;
                tempChange[{Xindex, Yindex}].type = target.type; // 紀錄翻轉
            }
            
        }*/
        if (canChange) {
            //cout << "cc" << tempChange.size() << endl;
            break;
        }
        else tempChange = constempChange; // 候選還原(不可改)
    }
}
// 檢查 Constrain 條件是否滿足
void Policy (const unordered_map<string, vector<pair<int, int>>> &cellSize,
             const unordered_map<string, vector<string>> &libName,
             const unordered_map<string, unordered_map<int, vector<int>>> &libnexts,
             const unordered_map<string, bool> &isExtending,
             const unordered_map<string, vector<vector<vector<vector<Extend>>>>> &ExtendingLength,
             unordered_map<string, vector<vector<vector<SparseTable>>>> &MatchingLength,
             const unordered_map<string, vector<vector<vector<vector<Point>>>>> &ExtendingPoint,
             const unordered_map<string, vector<vector<vector<vector<Point>>>>> &MatchingPoint,
             const unordered_map<string, unordered_set<pair<int,int>, pair_hash>> &cellIndex,
             vector<vector<ChipArray>> &chip) {
    unordered_map<pair<int,int>, int, pair_hash> tempchange;
    unordered_map<int, string> Type = {{0, "N"}, {1, "FN"}, {2, "S"}, {3, "FS"}};
    ofstream output("output.tcl");
    bool canChange;
    int count = 0, warning = 0;
    ChipArray tempcell;
    for (const auto& cell : cellIndex) {
        //if(isExtending.at(cell.first)) {         
        cout << "Cell name: " << cell.first << endl;
        for (const auto& item : cell.second) {
            canChange = true;
            //cout << chip[item.second][item.first].instance << " " << item.first << " " << item.second << endl;
            Judgment(cellSize, libnexts, isExtending, ExtendingLength, MatchingLength, ExtendingPoint, MatchingPoint, tempchange, item.first, item.second, canChange, chip);
            if(canChange) { 
                count++;
                for(auto &change: tempchange) {
                    tempcell = chip[change.first.second][change.first.first];
                    tempcell.libnum = change.second;
                    output << "change_link [get_cells " << tempcell.instance << "] " << libName.at(tempcell.cell)[tempcell.libnum] << "/" << tempcell.cell << endl;
                    //if (change.second.rotate) {
                    //    tempcell.type = change.second.type;
                    //    output << "rotate_cell [get_cells " << tempcell.instance << "] -orient " << Type.at(tempcell.type) << endl;
                    //}
                    chip[change.first.second][change.first.first] = tempcell;
                }
                if(tempchange.empty()) {
                    //cout << chip[item.second][item.first].instance << " no change" << endl;
                    warning++;
                }else output << endl;
                //cout << "ff" << tempchange.size() << endl;
                tempchange.clear();
            }
            
        }
        cout << "change num: " << count << "/" << cell.second.size() << endl;
        //cout << "warning num: " << warning << " matching cell" << endl;
        count = 0;
        warning = 0;
        //}
    }
    output.close();
}

// 主程式
int main(int argc, char* argv[]) {
    if(argc < 5) {
        cerr << "Usage: " << argv[0] << " input.csv input.def\n";
        return 0;
    }
    auto startTime = high_resolution_clock::now();// 取得起始時鐘

    // 存放 LEF 檔案中所有cell的尺寸
    unordered_map<string, vector<pair<int, int>>> cellSize;
    // 存放 lib 名稱
    unordered_map<string, vector<string>> libName;
    // 存放 CSV 檔案中所有cell的 M1Polys (Convex-Concave, x, y)
    unordered_map<string, vector<vector<vector<Poly>>>> cellM1Polys;
    readCSV(argc, argv, cellSize, cellM1Polys, libName);

    // 建立每個libnum節點對應的「下一節點集合」
    unordered_map<string, unordered_map<int, vector<int>>> libnexts;
    readNextsCSV(argv[argc-2], libnexts);

    // 翻轉順序: [0]為R0(N)、[1]為MY(FN)、[2]為R180(S)、[3]為MX(FS)
    // 方向順序: [0]為左、[1]為右、[2]為上、[3]為下

    // 存放 LEF 檔案中所有cell的 M1zone (矩形)
    unordered_map<string, vector<vector<Rect>>> cellM1Zone;
    // 紀錄凸型 cell的 Extending 點(Point)
    unordered_map<string, vector<vector<vector<vector<Point>>>>> ExtendingPoint;
    // 紀錄 cell 的 Matching 點(Point)
    unordered_map<string, vector<vector<vector<vector<Point>>>>> MatchingPoint;
    // 對每個 cell 添加 M1Zone、cornerEM points
    addM1Zone(cellSize, cellM1Polys, cellM1Zone, ExtendingPoint, MatchingPoint);

    // 紀錄是否為凸型 cell
    unordered_map<string, bool> isExtending;
    // 紀錄凸型 cell 的 Extending 邊(lmin、lmax、lout)(可能有多個)
    unordered_map<string, vector<vector<vector<vector<Extend>>>>> ExtendingLength;
    // 紀錄 cell 的 Matching 邊(SparseTable)
    unordered_map<string, vector<vector<vector<SparseTable>>>> MatchingLength;
    // 計算所有 cell 的Extending、Matching
    Scan(cellSize, cellM1Zone, isExtending, ExtendingLength, MatchingLength);

    for(auto &cell : isExtending){
        if(cell.second) cout << cell.first << " is extending cell" << endl;
    }
    // 存放 DEF 檔案中所有cell的類型、名稱、座標、方向
    vector<DEFfile> DEFdata;
    readDEF(argv[argc-1], DEFdata);

    // 依 index 存放 cell 名稱對應的類型、座標、左右端點
    vector<vector<ChipArray>> chip;
    // 存放 cell 類型對應的所有 cell index
    unordered_map<string, unordered_set<pair<int,int>, pair_hash>> cellIndex;
    // 製作表格
    createTable(cellSize, DEFdata, chip, cellIndex);
    Policy(cellSize, libName, libnexts, isExtending, ExtendingLength, MatchingLength, ExtendingPoint, MatchingPoint, cellIndex, chip);

    auto endTime = high_resolution_clock::now(); // 取得結束時鐘
    auto duration = duration_cast<microseconds>(endTime - startTime); // 計算執行時間
    cout << "Execution time: " << duration.count() << " microseconds" << endl;
    return 0;
} 