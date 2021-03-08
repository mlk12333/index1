//
//  main.cpp
//  index1
//
//  Created by 孟令凯 on 2020/12/21.
//

#include <iostream>
#include "string"
#include <fstream>
#include "Graph.hpp"
using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    if(argc<5) return 0;
    clock_t startTime,endTime;
    startTime = clock();//计时开始
    
    std::cout << "Hello, World!\n";
    string str = "/Users/milk/test_data/email-Eu-core2.txt";
    string f = "/Users/milk/test_data/index/";
    Graph graph(str,f);
    
//    fstream file;
//    file.open(str,ios::in);
//    int b,c;
//    vector<vector<int>> num;
//    while (!file.eof())            // 若未到文件结束一直循环
//    {
//        vector<int> a;
//
//        if(file>>b>>c){
//            a.push_back(b);
//            a.push_back(c);
//
//            num.push_back(a);
//        }
//
//
//    }
//    file.close();
//
//    sort(num.begin(),num.end());
//
//    int l = (int)num.size();
//
//    int i2 = 0;
//
//    vector<int> mm;
//
//    ofstream fout("/Users/milk/test_data/twitter_combined2.txt",ios::out);
//
//    for(int i = 0;i<l;i++){
//        if(num[i][0] == i2){
//            if(num[i][1]!=i2)
//               mm.push_back(num[i][1]);
//        }else{
//            sort(mm.begin(),mm.end());
//
//            for(int j = 0;j<mm.size();j++){
//                fout<<i2<<" "<<mm[j]<<"\n";
//            }
//            mm.clear();
//            i2 =num[i][0];
//            if(num[i][1]!=i2)
//               mm.push_back(num[i][1]);
//
//        }
//    }
//
//    fout.close();
//
    
    
//    startTime = clock();//计时开始
    
//    graph.creatIndex();
    graph.creatIndex(atof(argv[2]),atof(argv[3]),atoi(argv[4]));
    
//    graph.query(0.34, 0.4, 4);
    
//    endTime = clock();//计时结束
//    cout << "The creatIndex time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    
//    fstream file;
//    file.open(f+"nei_num.txt",ios::in);
//    if(!file){
//        graph.creatIndex();
//        endTime = clock();//计时结束
//
//        cout << "The creatIndex time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
//    }else{
//        graph.readIndex();
//
//        endTime = clock();//计时结束
//
//        cout << "The readIndex time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
//    }


//
//    while (true) {
//        double p1,p2;
//        int p3;
//
//        cin>>p1;
//        if(p1 == -1) break;
//        cin>>p2>>p3;
//
//        startTime = clock();
//
//        graph.query(p1, p2, p3);
//
//        endTime = clock();//计时结束
//
//        cout << "The query time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
//    }
    return 0;
}
