/*************************************************************************
	> File Name: main.cpp
	> Author: 
	> Mail: 
	> Created Time: 2016年12月14日 星期三 16时44分11秒
 ************************************************************************/

#include<iostream>
#include<vector>
#include<string>
#include "vertex_model.h"
using namespace std;

int main(){
	int nn = 0;
	double lx = 1.0, ly = 1.0;
    vertex_model vertex(30);
    vertex.initial();

	for( int i = 0; i< 90; i++)
	{
		string str0("E:\\work\\codes\\data\\votex_model2d\\");
		ly = ly - 0.01 / (2 * ly);
		lx = 1.0/ly;
		vertex.setRatio(lx, ly);
		nn = vertex.Honda(nn, str0);
		if( i % 3 ==0)
		{
			vertex.LongAxis();
			cout << "========================================== Division  "<<  i<<endl;
			vertex.Division();
		}
	}
}
