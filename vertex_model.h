/*************************************************************************
	> File Name: vertex_model.h
	> Author: 
	> Mail: 
	> Created Time: 2016年12月14日 星期三 18时43分06秒
 ************************************************************************/

#ifndef _VERTEX_MODEL_H
#define _VERTEX_MODEL_H
#endif

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<vector>
#include<cmath>
#include<algorithm>
#include<sstream>
#include<string>
#include <Eigen>
#include<Dense>
#include<Eigenvalues>
#include"vertices.h"
#include"cells.h"
using namespace std;
using Eigen::MatrixXd;
using Eigen::EigenSolver;

const double PI  =3.141592653589793238463;
class vertex_model
{
    public:
        vertex_model(int num):
            lx(1.0), ly(1.0),K_alpha(1.0), Gama_alpha(0.001), Lamda_ij(0.0)
        {
            cell_number = num;
        };
		void setRatio(double a, double b)
		{
			lx = a;
			ly = b;
		}
        void initial();
        void voronoi();
		void Division();
		void LongAxis();
		int Honda(int nn, string str0);
		void dump(const string filename1, const string filename2);
        void get();
		void getArea0(const vector<double> &x, const vector<double> &y);

		void getArea(const vector<double> &x, const vector<double> &y);
    private:
		double force(const vector<double> &pxnew, const vector<double> &pynew\
						,vector<double> &pfxnew, vector<double> &pfynew);
		void checkT1(int &sig, int &Nstep);
		void GetV2V_id();
		void GetV2C_id();
		void GetL2V_id();
		void SetCellLocs();
        double lx, ly;
        vector <vertices> node;
        vector <cells> cell;
        vector<vector <int>> C2V_id, V2V_id, V2C_id, L2V_id;
        vector<vector<double>> C2V_nodex, C2V_nodey, Areadx, Aready, Perdx, Perdy;
        vector<double> xsolve_unique, ysolve_unique, lijdx, lijdy;
		const double K_alpha, Gama_alpha, Lamda_ij;
		// K_alpha: area elastic modulus which maintain the area of cell
		// Gama_alpha: perimeter coefficient which minimizes the perimeter of cell
		// Lamda_ij: bond tension parameter. It can take positive or negative values
		// 			accroding to wether surface tension or adhesion dominates the boundaries
        int cell_number;
};

void vertex_model::LongAxis()
{
	double maxr = 0.0;
	int indexi, indexj, index0, index1, index;
	double x0ij, y0ij, x1ij, y1ij, longx, longy, longlength;
	vector<vector<double>> loop;
	vector<double> temp;
	
	fstream outfile;

	
	double dr = 0.001;
	for(int c = 0; c<C2V_id.size(); c++ )
	{
		loop.clear();
		for(int i = 0; i< C2V_id[c].size(); i++)
		{
			temp.clear();
			index0 = C2V_id[c][i];
			index1 = C2V_id[c][(i+1)%C2V_id[c].size()];
			x0ij = xsolve_unique[index0];
			y0ij = ysolve_unique[index0];
			x1ij = xsolve_unique[index1];
			y1ij = ysolve_unique[index1];
			if(abs(cell[c].getx() - x0ij) > lx/2.0)
				x0ij = x0ij + abs(cell[c].getx())/cell[c].getx() * lx;
			if(abs(cell[c].gety() - y0ij) > ly/2.0)
				y0ij = y0ij + abs(cell[c].gety())/cell[c].gety() * ly;
			if(abs(cell[c].getx() - x1ij) > lx/2.0)
				x1ij = x1ij + abs(cell[c].getx())/cell[c].getx() * lx;
			if(abs(cell[c].gety() - y1ij) > ly/2.0)
				y1ij = y1ij + abs(cell[c].gety())/cell[c].gety() * ly;
			
			double rij = sqrt(pow((x1ij - x0ij), 2) + pow((y1ij- y0ij), 2));
			int npoint = rij/ dr;
			double dx = (x1ij - x0ij) / rij * dr;
			double dy = (y1ij - y0ij) / rij * dr;
			temp.push_back(x0ij); temp.push_back(y0ij);
			loop.push_back(temp);
			temp.clear();

			for(int j = 0; j < npoint; j++)
			{
				temp.push_back(x0ij + (j+1)*dx); temp.push_back(y0ij + (j+1)*dy);
				loop.push_back(temp);
				temp.clear();

			}
		}
		
		
		double xx = 0.0, yy = 0.0, xy = 0.0, yx = 0.0;
		temp.resize(2);
		temp[0] = 0.0; temp[1] = 0.0;
		for(int i = 0; i< loop.size(); i++)
		{
			temp[0] += loop[i][0];
			temp[1] += loop[i][1];
		}
		temp[0] /= loop.size();
		temp[1] /= loop.size();
		
	
		
		for(int i = 0; i< loop.size(); i++)
		{
			x0ij = loop[i][0] - temp[0];
			y0ij = loop[i][1] - temp[1];
			xx += pow(x0ij, 2); yy += pow(y0ij, 2);
			xy += x0ij * y0ij; yx += x0ij * y0ij;
		}

		MatrixXd m = MatrixXd::Zero(2, 2);

		m(0,0) = xx;
		m(1,0) = xy;
		m(0,1) = yx;
		m(1,1) = yy;
		EigenSolver<MatrixXd> es(m);
//		cout << es.pseudoEigenvectors()<<endl;

//		cout <<"eigenvalue matrix "<<endl;cout <<es.pseudoEigenvalueMatrix()<<endl;

		temp.clear();
		temp.push_back(es.eigenvalues().real()[0]); temp.push_back(es.eigenvalues().real()[1]);
		if(temp[0] > temp[1])
			index = 1;
		else
			index = 0;
		
		cell[c].setLongx(es.eigenvectors().col((index + 1) % temp.size()).real()[0]);
		cell[c].setLongy(es.eigenvectors().col((index + 1) % temp.size()).real()[1]);
		cell[c].setLongLength(temp[(index+1)%temp.size()]/temp[index]);
		if(c == 26)
		{
			outfile.open("E:\\work\\codes\\data\\delete.txt", ios::out);
			cout << "loop "<< endl;
			for(auto it1: loop)
			{
				for(auto it2: it1)
				{
					outfile<< it2<< "  ";
				}
				outfile << endl;
			}
			outfile.close();
			cout << endl;
			cout << es.eigenvectors().real() << endl;
			cout << " eigenvalue matrix " <<endl; cout << es.pseudoEigenvalueMatrix()<<endl;
			
			cout << cell[c].getLongx() <<"  "<< cell[c].getLongy()<< "   "<< (index +1)%temp.size()<< "   "<< temp.size()<<endl;	
		}
	}

	
/*	for(int c = 0; c < C2V_id.size(); c++)
	{
		for(int i = 0; i < C2V_id[c].size()-1; i++)
		{
			maxr = 0.0;
			for(int j = i+1; j < C2V_id[c].size(); j++)
			{
				index0 = C2V_id[c][i];
				index1 = C2V_id[c][j];
				x0ij = xsolve_unique[index0];
				y0ij = ysolve_unique[index0];
				x1ij = xsolve_unique[index1];
				y1ij = ysolve_unique[index1];
				if(abs(cell[c].getx() - x0ij) > lx/2.0)
					x0ij = x0ij + abs(cell[c].getx())/cell[c].getx() * lx;
				if(abs(cell[c].gety() - y0ij) > ly/2.0)
					y0ij = y0ij + abs(cell[c].gety())/cell[c].gety() * ly;
				if(abs(cell[c].getx() - x1ij) > lx/2.0)
					x1ij = x1ij + abs(cell[c].getx())/cell[c].getx() * lx;
				if(abs(cell[c].gety() - y1ij) > ly/2.0)
					y1ij = y1ij + abs(cell[c].gety())/cell[c].gety() * ly;
				
				double r = sqrt(pow(x0ij - x1ij, 2) + pow(y0ij - y1ij, 2));
				if(r > maxr)
				{
					maxr = r;
					indexi = i; indexj = j;
					longx = x0ij - x1ij;
					longy = y0ij - y1ij;
					longlength = maxr;
				}
			}
		}
		cell[c].setLongx(longx);
		cell[c].setLongy(longy);
		cell[c].setLongLength(longlength);
	}*/
	
}

void vertex_model::Division()
{
	int Divi_C, index0, index1;
	vector<int> :: iterator it;
	double maxr = 0, x0ij, x1ij, y0ij, y1ij;
	for(int i = 0; i< C2V_id.size(); i++)
	{
		if( maxr < cell[i].getLongLength())
		{
			maxr = cell[i].getLongLength();
			Divi_C = i;
		}
	}
	
	cout << "Divi_C  "<< Divi_C<<endl;

	vector<int> temp_id_insert, temp_id0, temp_id1, index_record;
	temp_id_insert.clear();
	temp_id0.clear();
	temp_id1.clear();
	index_record.clear();
	
	double slop = - cell[Divi_C].getLongx()/cell[Divi_C].getLongy();
	for(int j = 0; j < C2V_id[Divi_C].size(); j++)
	{
		int index0 = C2V_id[Divi_C][j];
		int index1 = C2V_id[Divi_C][(j+1)%C2V_id[Divi_C].size()];
		
		x0ij = xsolve_unique[index0];
		y0ij = ysolve_unique[index0];
		x1ij = xsolve_unique[index1];
		y1ij = ysolve_unique[index1];
		if(abs(cell[Divi_C].getx() - x0ij) > lx/2.0)
			x0ij = x0ij + abs(cell[Divi_C].getx())/cell[Divi_C].getx() * lx;
		if(abs(cell[Divi_C].gety() - y0ij) > ly/2.0)
			y0ij = y0ij + abs(cell[Divi_C].gety())/cell[Divi_C].gety() * ly;
		if(abs(cell[Divi_C].getx() - x1ij) > lx/2.0)
			x1ij = x1ij + abs(cell[Divi_C].getx())/cell[Divi_C].getx() * lx;
		if(abs(cell[Divi_C].gety() - y1ij) > ly/2.0)
			y1ij = y1ij + abs(cell[Divi_C].gety())/cell[Divi_C].gety() * ly;
		
		x0ij = x0ij - cell[Divi_C].getx();
		x1ij = x1ij - cell[Divi_C].getx();
		y0ij = y0ij - cell[Divi_C].gety();
		y1ij = y1ij - cell[Divi_C].gety();
		
		double cross0 = x0ij * slop - y0ij * 1;
		double cross1 = x1ij * slop - y1ij * 1;

		if(cross0 * cross1 < 0)
		{
			double k = (y1ij - y0ij)/(x1ij - x0ij);
			double tempx = (y0ij - k * x0ij)/(slop - k);
			double tempy = slop * tempx;
			tempx = tempx + cell[Divi_C].getx();
			tempy = tempy + cell[Divi_C].gety();
			
			xsolve_unique.push_back(tempx);
			ysolve_unique.push_back(tempy);
			
			temp_id_insert.push_back(j);
			
			index_record.push_back(index0);
			
		}
	}
	
	
	
	for(int i = 0; i< index_record.size(); i++)
	{
		vector<int> :: iterator temp_it;
		temp_it = find(V2C_id[index_record[i]].begin(), V2C_id[index_record[i]].end(), Divi_C);
		int position = distance(V2C_id[index_record[i]].begin(), temp_it ) + V2C_id[index_record[i]].size();
		position = (position -1 ) % V2C_id[index_record[i]].size();	
		int cell_index = V2C_id[index_record[i]][position];
/*		 if(i==0)
		 {
			 cout<<" cell index is " << cell_index<<"   " ;
			 for(auto it1: V2C_id[index_record[i]])
				 cout <<it1<< "   ";
		 }
		 
		 if(i==1)
		 {
			 cout<<" cell index is " << cell_index<<"   " ;
			 for(auto it1: V2C_id[index_record[i]])
				 cout <<it1<< "   ";
		 }
		 cout << endl;*/

		temp_it = find(C2V_id[cell_index].begin(), C2V_id[cell_index].end(), index_record[i]);
		C2V_id[cell_index].insert(temp_it, xsolve_unique.size() - (2 - i));
	}
	
	
//	cout << " size is " << C2V_id[Divi_C].size()<<" position is ";
//	cout << temp_id_insert[0]<< "  "<<temp_id_insert[1]<<endl;
	int solve_id = xsolve_unique.size() - 2;

	for(int j = 0; j < 2; j++)
	{
		it = C2V_id[Divi_C].begin();
		int position = temp_id_insert[j] +j;
		position = (position + 1)%C2V_id[Divi_C].size();
		C2V_id[Divi_C].insert(it + position, solve_id);
		solve_id++;

	}
	
	vector<int> :: iterator it_ini = find(C2V_id[Divi_C].begin(), C2V_id[Divi_C].end(), xsolve_unique.size() - 2);
	vector<int> :: iterator it_end = find(C2V_id[Divi_C].begin(), C2V_id[Divi_C].end(), xsolve_unique.size() - 1);
	index1 = *it_end;
	int displace = distance(C2V_id[Divi_C].begin(), it_ini);
	for(int i = 0; i< C2V_id[Divi_C].size(); i++)
	{
		int position = (i + displace)% C2V_id[Divi_C].size();
//		cout << i << " "<< *(it_ini + position - displace)<<endl;
		index0 = *(it_ini + position - displace);
		temp_id0.push_back(index0);
		if(index0 == index1)
			break;
	}
	it_ini = find(C2V_id[Divi_C].begin(), C2V_id[Divi_C].end(), xsolve_unique.size() - 2);
	it_end = find(C2V_id[Divi_C].begin(), C2V_id[Divi_C].end(), xsolve_unique.size() - 1);
	index0 = *it_ini;
	displace = distance(C2V_id[Divi_C].begin(), it_end);
	for(int i= 0; i< C2V_id[Divi_C].size(); i++)
	{
		int position = (i + displace)% C2V_id[Divi_C].size();
//		cout << i << " "<< *(it_end + position - displace)<<endl;
		index1 = *(it_end + position - displace);
		temp_id1.push_back(index1);
		if(index0 == index1)
			break;
	}
	C2V_id[Divi_C] = temp_id0;	
	C2V_id.push_back(temp_id1);

	cells temp_cell(cell[Divi_C].getx(), cell[Divi_C].gety());
	cell.push_back(temp_cell);
	
	SetCellLocs();
	
//=======================get V2C_id=======================================
	GetV2C_id();
//=======================get V2V_id=======================================
	GetV2V_id();
//=======================get L2V_id=======================================
	GetL2V_id();	
	//refresh the lists.
}

void vertex_model::SetCellLocs()
{
	double x, y;
	for(int i = 0; i< C2V_id.size(); i++)
	{
		x = 0.0;
		y = 0.0;
		double centx = 0.0;
		double centy = 0.0;
		double area = 0.0;
		for(int j = 0; j< C2V_id[i].size(); j++)
		{
			int id = C2V_id[i][j];
			int id1 = C2V_id[i][(j+1)%C2V_id[i].size()];
			
			xsolve_unique[id] -= round(xsolve_unique[id]/lx)*lx;
			ysolve_unique[id] -= round(ysolve_unique[id]/ly)*ly;
			
			
			double x0ij = xsolve_unique[id];
			double y0ij = ysolve_unique[id];
			double x1ij = xsolve_unique[id1];
			double y1ij = ysolve_unique[id1];
			
			
			if(abs(cell[i].getx() - x0ij)>lx/2.0)
				x0ij += abs(cell[i].getx())/cell[i].getx() * lx;
			if(abs(cell[i].gety() - y0ij)>ly/2.0)
				y0ij += abs(cell[i].gety())/cell[i].gety() * ly;
			if(abs(cell[i].getx() - x1ij)>lx/2.0)
				x1ij += abs(cell[i].getx())/cell[i].getx() * lx;
			if(abs(cell[i].gety() - y1ij)>ly/2.0)
				y1ij += abs(cell[i].gety())/cell[i].gety() * ly;
			
			area += (x0ij * y1ij - x1ij * y0ij)/2.0;			
			centx += (x0ij + x1ij)*(x0ij * y1ij - x1ij * y0ij);
			centy += (y0ij + y1ij)*(x0ij * y1ij - x1ij * y0ij);

		}
		x = centx / (6 * area);
		y = centy/(6 * area);
		x -= round(x/lx)*lx;
		y -= round(y/ly)*ly;
//		x = cell[i].getx() - round(cell[i].getx()/lx)*lx;
//		y = cell[i].gety() - round(cell[i].gety()/ly)*ly;
		cell[i].set(x, y);
	}
}

void vertex_model::dump(const string filename1, const string filename2)
{
	fstream outfile;
	double x0ij, y0ij, x1ij, y1ij;
	int index0, index1;
	outfile.open(filename1, ios::out);
	int temp_count = 0;
	for(int i = 0; i< C2V_id.size(); i++)
	{
		for(int j = 0; j< C2V_id[i].size() - 1; j++)
		{
			index0 = C2V_id[i][j];
			index1 = C2V_id[i][j+1];
			x0ij = xsolve_unique[index0];
			y0ij = ysolve_unique[index0];
			x1ij = xsolve_unique[index1];
			y1ij = ysolve_unique[index1];
			if(abs(cell[i].getx() - x0ij) > lx/2.0)
				x0ij = x0ij + abs(cell[i].getx())/cell[i].getx() * lx;
			if(abs(cell[i].gety() - y0ij) > ly/2.0)
				y0ij = y0ij + abs(cell[i].gety())/cell[i].gety() * ly;
			if(abs(cell[i].getx() - x1ij) > lx/2.0)
				x1ij = x1ij + abs(cell[i].getx())/cell[i].getx() * lx;
			if(abs(cell[i].gety() - y1ij) > ly/2.0)
				y1ij = y1ij + abs(cell[i].gety())/cell[i].gety() * ly;
			outfile<< x0ij<<","<<y0ij<<","<<x1ij<<","<<y1ij<<endl;
			//temp_count++;
		}
		index0 = C2V_id[i][C2V_id[i].size() - 1];
		index1 = C2V_id[i][0];
		x0ij = xsolve_unique[index0];
		y0ij = ysolve_unique[index0];
		x1ij = xsolve_unique[index1];
		y1ij = ysolve_unique[index1];
		if(abs(cell[i].getx() - x0ij) > lx/2.0)
			x0ij += abs(cell[i].getx())/cell[i].getx() * lx;
		if(abs(cell[i].gety() - y0ij) > ly/2.0)
			y0ij += abs(cell[i].gety())/cell[i].gety() * ly;
		if(abs(cell[i].getx() - x1ij) > lx/2.0)
			x1ij += abs(cell[i].getx())/cell[i].getx() * lx;
		if(abs(cell[i].gety() - y1ij) > ly/2.0)
			y1ij += abs(cell[i].gety())/cell[i].gety() * ly;
		outfile<< x0ij<<","<<y0ij<<","<<x1ij<<","<<y1ij<<endl;
		//temp_count++;
	}
	outfile.close();
	
	outfile.open(filename2, ios::out);
	for(auto jter= cell.begin(); jter != cell.end(); jter++)
	{
		outfile<<jter->getx()<<","<< jter->gety()<<endl; 
	}
	outfile.close();
}

double vertex_model::force(const vector<double> &pxnew, const vector<double> &pynew\
						,vector<double> &pfxnew, vector<double> &pfynew)
{
	int cid;
	vector<int> :: iterator it;
	double force_precision(0.0);

	double vn(pxnew.size());
	
	getArea0(pxnew, pynew);
	getArea(pxnew, pynew);
	

	

	for(int i = 0; i< vn; i++)
	{
		pfxnew[i] = 0;
		pfynew[i] = 0;
		for(int j = 0; j < V2C_id[i].size();j++)
		{
			cid = V2C_id[i][j];
			pfxnew[i] -= K_alpha*(cell[cid].getArea() - cell[cid].getArea0()) * Areadx[i][j]\
						+Gama_alpha * (cell[cid].getPerimeter() - cell[cid].getPerimeter0())* Perdx[i][j];
			pfynew[i] -= K_alpha*(cell[cid].getArea() - cell[cid].getArea0()) * Aready[i][j]\
						+Gama_alpha * (cell[cid].getPerimeter() - cell[cid].getPerimeter0())* Perdy[i][j];
						
			if(abs(K_alpha) > 1000 | abs(cell[cid].getArea()) > 1000 | abs(Areadx[i][j]) > 1000 | abs(Perdx[i][j]) > 1000\
			| abs(cell[cid].getPerimeter0())> 1000| abs(cell[cid].getPerimeter()) > 1000| abs(cell[cid].getArea0()) > 1000)
				cout <<" i " <<i << " j "<<j << " area0 "<<cell[cid].getArea0() <<endl;
		}
	}
	
/*	for( int i = 0; i< pfxnew.size(); i++)
	{
		if(abs(pfxnew[i]) > 10000)
		{
			cout << i << "  "<< pfxnew[i] << "  xsolve_unique "<< xsolve_unique[i] << endl;
		}
	}
	
	for(int i = 0; i< Areadx.size(); i++)
	{
		for(int j = 0; j< Areadx[i].size(); j++)
		{
			if(abs(Areadx[i][j]) > 10000)
			{
				cout << " i "<< i<<" j " << j << " Areadx "<< Areadx[i][j]<<endl;
			}
		}
	}*/
	
	for(int i= 0; i< pfxnew.size(); i++)
	{
		pfxnew[i] += lijdx[i];
		pfynew[i] += lijdy[i];
		force_precision += sqrt(pow(pfxnew[i], 2) + pow(pfynew[i], 2));
	}
	
	return force_precision/pfxnew.size();

}

void vertex_model::checkT1( int &sig, int &Nstep)
{
	double x0ij, y0ij, x1ij, y1ij, vvr;
	int index0, index1;
	int vn(L2V_id.size()), count(0);
	

	
	for(int i = 0; i<L2V_id.size(); i++)
	{
		index0 = L2V_id[i][0];
		index1 = L2V_id[i][1];
		x0ij = xsolve_unique[index0];
		y0ij = ysolve_unique[index0];
		x1ij = xsolve_unique[index1];
		y1ij = ysolve_unique[index1];
		if(abs(x0ij - x1ij) > lx/2.0)
			x0ij = x0ij + abs(x1ij)/x1ij;
		if(abs(y0ij - y1ij) > ly/2.0)
			y0ij = y0ij + abs(y1ij)/y1ij;
		vvr = sqrt(pow(x0ij - x1ij, 2) + pow(y0ij - y1ij, 2));
		if(vvr < 0.001)
		{
			sig = 0;
			Nstep = 0;
			count++;
			xsolve_unique[index0] = - (y0ij - y1ij)/2.0+(x0ij + x1ij)/2.0;
			ysolve_unique[index0] = (x0ij - x1ij)/2.0+(y0ij + y1ij)/2.0;
			xsolve_unique[index1] = -(y1ij - y0ij)/2.0+(x0ij + x1ij)/2.0;
			ysolve_unique[index1] = (x1ij - x0ij)/2.0+(y0ij + y1ij)/2.0;
			vector<int>::iterator it0, it1, it2, it0_record, it1_record;
			int cid1, cid1_record, cid0_record, cid1_switch(0), cid0_switch(0);
			int sub_index0, sub_index1;

			for(int m = 0; m < V2C_id[index1].size(); m++)
			{
				cid1 = V2C_id[index1][m];
				it0 = find(C2V_id[cid1].begin(), C2V_id[cid1].end(), index0);
				if(it0 == C2V_id[cid1].end())
				{
					if(find(C2V_id[cid1].begin(), C2V_id[cid1].end(), index1) != C2V_id[cid1].end())
					{	
						it2 = find(C2V_id[cid1].begin(), C2V_id[cid1].end(), index1);
						cid1_record = cid1;
					}	
				}
			}
			
			for(int k = 0; k<V2C_id[index0].size(); k++)
			{
				int cid = V2C_id[index0][k];
				it0 = find(C2V_id[cid].begin(), C2V_id[cid].end(), index0);
				it1 = find(C2V_id[cid].begin(), C2V_id[cid].end(), index1);
				if(it1 == C2V_id[cid].end())
				{
					//C2V_id[cid].insert(it0+1, index1);
					it0_record = it0;
					cid0_record = cid;
				}
				else
				{
					if(distance(C2V_id[cid].begin(), it0) - distance(C2V_id[cid].begin(), it1)==1\
					|distance(C2V_id[cid].begin(), it0) - distance(C2V_id[cid].begin(), it1)==1-C2V_id[cid].size())
					{
						C2V_id[cid].erase(it1);
						cid0_switch = cid;
					}
					else if(distance(C2V_id[cid].begin(), it0) - distance(C2V_id[cid].begin(), it1)==-1\
					|distance(C2V_id[cid].begin(), it0) - distance(C2V_id[cid].begin(), it1)==C2V_id[cid].size()-1)
					{
						C2V_id[cid].erase(it0);
						cid1_switch = cid;
					}
				}
			}


			C2V_id[cid0_record].insert(it0_record+1, index1);

			C2V_id[cid1_record].insert(it2+1, index0);

			//swap(V2C_id[index0][cid1_record], V2C_id[index0][cid0_switch]);
			//swap(V2C_id[index1][cid0_record], V2C_id[index1][cid1_switch]);

		}
		
	}

	GetV2C_id();
	GetV2V_id();
	GetL2V_id();
}

int vertex_model::Honda(int nn, string str0)
{
	vector<double> &xnew = xsolve_unique, &ynew = ysolve_unique;
	vector<double> xold(xnew), yold(ynew);
	vector<double> fxold(xsolve_unique.size()), fxnew(xsolve_unique.size());
	vector<double> fyold(xsolve_unique.size()), fynew(xsolve_unique.size());
	double gamma_n(1.0e-6), temp(0.0);
	double mean_force;
	stringstream ss1, ss2;
	int sig(0), Nstep(0);
	string str1, str2;
	//str0("E:\\work\\codes\\data\\votex_model2d\\");

	int vn(xsolve_unique.size());


	
	mean_force= force(xnew, ynew, fxnew, fynew);
	cout << " iniforce "<< mean_force<<endl;
//	checkT1();

	for(int step = 0; step < 1000; step++)
	{
		if(step %10000 ==0 | mean_force < 1.0e-8)
		{
			if(step > 0)
			{
				SetCellLocs();
				ss1.str("");
				ss2.str("");
				ss1 << "edges_"<< ++nn;
				ss2 << "cells_"<< nn;
				str1 = str0 + ss1.str() + ".csv";
				str2 = str0 + ss2.str() + ".csv";
				dump(str1, str2);
				cout << "nn is "<< nn;
				cout<<"   force and gamma_n "<< mean_force<<"  "<< gamma_n<<endl;
			}

			if( mean_force < 1.0e-6)
			{
				cout << " finished "<< mean_force<<endl;
				GetV2C_id();
				GetL2V_id();
				GetV2V_id();
				return nn;
			}
		
		}
		
		
		if(step % 10 ==0& step>11)
		{
		checkT1( sig, Nstep );
		}

//		if(Nstep>20 & sig == 0)
//		{
//			sig = 1;
//			Nstep++;		
//		}
		for(int j = 0; j< vn; j++)
		{
			xold[j] = xnew[j];
			yold[j] = ynew[j];
			fxold[j] = fxnew[j];
			fyold[j] = fynew[j];
			xnew[j] += gamma_n * fxnew[j];
			ynew[j] += gamma_n * fynew[j];
		}


		mean_force = force(xnew, ynew, fxnew, fynew);


		gamma_n = 0;
		temp = 0;

		for(int i = 0; i<vn; i++)
		{
			gamma_n -= (xnew[i] - xold[i]) * (fxnew[i] - fxold[i])\
					 + (ynew[i] - yold[i]) * (fynew[i] - fyold[i]);
		}

	
		for(int i = 0; i<vn; i++)
		{
			temp += pow((fxnew[i] - fxold[i]), 2) + pow((fynew[i] - fyold[i]), 2);
		}
		gamma_n /= temp;
		//cout << "  gamma_n is "<<gamma_n <<endl;
	}

}

void vertex_model::getArea0(const vector<double> &x, const vector<double> &y)
{
	int nb_id0, nb_id1, indexc, indexv, index0, index1, index;
	double x0ij, y0ij, x1ij, y1ij, area, perimeter, dcv0, dcv1, tempLijdx, tempLijdy;
	vector<double> tempAdx, tempAdy, tempPdx, tempPdy;
	vector<int> :: iterator it;
	
	for(int i = 0; i< C2V_id.size(); i++)
	{
		area = 0;
		perimeter = 0;
	
		for(int j = 0; j<C2V_id[i].size(); j++)
		{
			nb_id0 = C2V_id[i][j%C2V_id[i].size()];
			nb_id1 = C2V_id[i][(j+1)%C2V_id[i].size()];

			x0ij = x[nb_id0];
			x1ij = x[nb_id1];
			y0ij = y[nb_id0];
			y1ij = y[nb_id1];
			
			if(abs(cell[i].getx() - x0ij)>lx/2)
				x0ij = x0ij + abs(cell[i].getx())/cell[i].getx() * lx;
			if(abs(cell[i].gety() - y0ij)>ly/2)
				y0ij = y0ij + abs(cell[i].gety())/cell[i].gety() * ly;
			if(abs(cell[i].getx() - x1ij)>lx/2)
				x1ij = x1ij + abs(cell[i].getx())/cell[i].getx() * lx;
			if(abs(cell[i].gety() - y1ij)>ly/2)
				y1ij = y1ij + abs(cell[i].gety())/cell[i].gety() * ly;
			area += (x0ij * y1ij - x1ij * y0ij)/2.0;
			perimeter += sqrt(pow((x0ij - x1ij), 2) + pow((y0ij - y1ij), 2));
		}
		cell[i].setArea0(0);
		cell[i].setPerimeter0(0);
	}
	double sumarea(0);
}



void vertex_model::getArea(const vector<double> &x, const vector<double> &y)
{
	int nb_id0, nb_id1, indexc, indexv, index0, index1, index;
	double x0ij, y0ij, x1ij, y1ij, xij, yij, area, perimeter, dcv0, dcv1, tempLijdx, tempLijdy;
	double centx, centy;
	vector<double> tempAdx, tempAdy, tempPdx, tempPdy, tempCentx, tempCenty;
	vector<int> :: iterator it;
	
	Areadx.clear();
	Aready.clear();
	Perdx.clear();
	Perdy.clear();
	lijdx.clear();
	lijdy.clear();
	
	for(int i = 0; i< C2V_id.size(); i++)
	{
		area = 0;
		perimeter = 0;
		centx = 0;
		centy = 0;
	
		for(int j = 0; j<C2V_id[i].size(); j++)
		{
			nb_id0 = C2V_id[i][j%C2V_id[i].size()];
			nb_id1 = C2V_id[i][(j+1)%C2V_id[i].size()];

			x0ij = x[nb_id0];
			x1ij = x[nb_id1];
			y0ij = y[nb_id0];
			y1ij = y[nb_id1];
			
			if(abs(cell[i].getx() - x0ij)>lx/2)
				x0ij = x0ij + abs(cell[i].getx())/cell[i].getx() * lx;
			if(abs(cell[i].gety() - y0ij)>ly/2)
				y0ij = y0ij + abs(cell[i].gety())/cell[i].gety() * ly;
			if(abs(cell[i].getx() - x1ij)>lx/2)
				x1ij = x1ij + abs(cell[i].getx())/cell[i].getx() * lx;
			if(abs(cell[i].gety() - y1ij)>ly/2)
				y1ij = y1ij + abs(cell[i].gety())/cell[i].gety() * ly;
			area += (x0ij * y1ij - x1ij * y0ij)/2.0;
			centx += (x0ij + x1ij)*(x0ij * y1ij - x1ij * y0ij);
			centy += (y0ij + y1ij)*(x0ij * y1ij - x1ij * y0ij);
			perimeter += sqrt(pow((x0ij - x1ij), 2) + pow((y0ij - y1ij), 2));
		}
		centx /= 6.0*area;
		centy /= 6.0*area;
		cell[i].set(centx, centy);
		cell[i].setArea(area);
//		cell[i].setArea0(rea);//should not be initialized dynamically
		cell[i].setPerimeter(perimeter);
//		cell[i].setPerimeter0(2.0*sqrt(PI* area));
	}

	
	for(int j = 0; j< V2C_id.size(); j++)
	{
		tempAdx.clear();
		tempAdy.clear();
		tempPdx.clear();
		tempPdy.clear();
		for(int k = 0; k < V2C_id[j].size(); k++)
		{
			indexc = V2C_id[j][k];
			it = find(C2V_id[indexc].begin(), C2V_id[indexc].end(), j);
			indexv = distance(C2V_id[indexc].begin(), it);
			index0 = C2V_id[indexc][(indexv+1+C2V_id[indexc].size())%C2V_id[indexc].size()];
			index  = C2V_id[indexc][(indexv + C2V_id[indexc].size())%C2V_id[indexc].size()];
			index1 = C2V_id[indexc][(indexv-1+C2V_id[indexc].size())%C2V_id[indexc].size()];
			
			x0ij = x[index0];
			y0ij = y[index0];
			x1ij = x[index1];
			y1ij = y[index1];
			xij = x[index];
			yij = y[index];
			
			
			if(abs(cell[indexc].getx() - x0ij)>lx/2)
				x0ij = x0ij + abs(cell[indexc].getx())/cell[indexc].getx() * lx;
			if(abs(cell[indexc].gety() - y0ij)>ly/2)
				y0ij = y0ij + abs(cell[indexc].gety())/cell[indexc].gety() * ly;
			if(abs(cell[indexc].getx() - x1ij)>lx/2)
				x1ij = x1ij + abs(cell[indexc].getx())/cell[indexc].getx() * lx;
			if(abs(cell[indexc].gety() - y1ij)>ly/2)
				y1ij = y1ij + abs(cell[indexc].gety())/cell[indexc].gety() * ly;
			if(abs(cell[indexc].getx() - xij)>lx/2)
				xij = xij + abs(cell[indexc].getx())/cell[indexc].getx() * lx;
			if(abs(cell[indexc].gety() - yij)>ly/2)
				yij = yij + abs(cell[indexc].gety())/cell[indexc].gety() * ly;
			
			tempAdx.push_back((y0ij - y1ij)/2.0);
			tempAdy.push_back((x1ij - x0ij)/2.0);
			dcv0 = sqrt(pow(xij - x0ij, 2) +\
						pow(yij - y0ij, 2));
			dcv1 = sqrt(pow(xij - x1ij, 2)+\
						pow(yij - y1ij, 2));
/*		if(abs(dcv0) >1000| abs(dcv1)>1000)
		{
			cout << "dcv0 "<<dcv0<< "  dcv1 "<<dcv1<<endl;
		}*/

			tempPdx.push_back((xij - x0ij)/dcv0 + (xij - x1ij)/dcv1);
			tempPdy.push_back((yij - y0ij)/dcv0 + (yij - y1ij)/dcv1);
		}
		Areadx.push_back(tempAdx);
		Aready.push_back(tempAdy);
		Perdx.push_back(tempPdx);
		Perdy.push_back(tempPdy);
	}
	
//======================================================================================================	
	
/*	for(int i = 0; i< Areadx.size(); i++)
	{
		for(int j =0; j < Areadx[i].size(); j++ )
		{
			if( abs(Areadx[i][j]) > 10000)
			{
			cout << " i "<< i << " j "<< j <<" Areadx "<< Areadx[i][j]<<endl;
			}
		}
	}
	
	for(int i = 0; i< Perdx.size(); i++)
	{
		for(int j =0; j < Perdx[i].size(); j++ )
		{
			if( abs(Perdx[i][j]) > 10000)
			{
			cout << " i "<< i << " j "<< j <<" Perdx "<< Perdx[i][j]<<endl;
			}
		}
	}*/
//======================================================================================================
	for(int i = 0; i< V2V_id.size(); i++)
	{
		tempLijdx = 0;
		tempLijdy = 0;
		for(int j = 0; j < V2V_id[i].size(); j++)
		{
			index = V2V_id[i][j];
			x0ij = xsolve_unique[index];
			y0ij = ysolve_unique[index];
			if(abs(x0ij - xsolve_unique[i])>lx/2.0)
				x0ij = x0ij + abs(xsolve_unique[i])/xsolve_unique[i];
			if(abs(y0ij - ysolve_unique[i])>ly/2.0)
				y0ij = y0ij + abs(ysolve_unique[i])/ysolve_unique[i];
			
			dcv0 = sqrt(pow(x0ij - xsolve_unique[i], 2)+\
						pow(y0ij - ysolve_unique[i], 2));
			tempLijdx += Lamda_ij * (x0ij - xsolve_unique[i])/dcv0;
			tempLijdy += Lamda_ij * (y0ij - ysolve_unique[i])/dcv0;
		}
		lijdx.push_back(tempLijdx);
		lijdy.push_back(tempLijdy);
	}


}
//=========================================================
//==========================================================
void vertex_model::initial()
{
    double temp0, temp1;
    cells temp_cell(0,0);
    srand(5);
    for(int i = 0; i< cell_number; i++)
    {
        temp0 = (rand()%10000)/10000.0 - 0.5;
        temp1 = (rand()%10000)/10000.0 - 0.5;
        temp_cell.set(temp0, temp1);
        cell.push_back(temp_cell);
    }
    //get();
    voronoi();
	getArea0(xsolve_unique, ysolve_unique);
	//getArea(xsolve_unique, ysolve_unique);

    //cout << "cell number is "<< cell_number<<endl;

}

void vertex_model::voronoi()
{
    int i = 0, countn = 0, index = 0, length = 0, ii = 0, jj = 0, kk = 0;
    double x0, y0, x1, y1, c0, c1;
    vector <double> xsolve, ysolve;
    vector <int> id0, id1; 
    vector <int> tempid0, tempid1;
    vector <double> tempsolvex, tempsolvey;
	fstream outfile;
    for(vector<cells>::iterator iter = cell.begin(); iter != cell.end(); iter++)
    {
        id0.clear();
        id1.clear();
        tempid0.clear(); // recod the id of the lines that intersecte at solution (xsolve, ysolve)
        tempid1.clear();
        tempsolvex.clear();
        tempsolvey.clear();
        xsolve.clear(); // record all solutions of all midperpendicular lines
        ysolve.clear();
        for(auto jter = cell.begin(); jter != cell.end(); jter++) // the loop among all midperpendicular lines about
        {                                                         // cell iter
            if(jter != iter)
            {
                for(auto kter = jter + 1; kter != cell.end(); kter++)
                {
                    if(kter != iter)
                    {
                        x0 = jter->getx() - iter->getx();
                        y0 = jter->gety() - iter->gety();
                        x1 = kter->getx() - iter->getx();
                        y1 = kter->gety() - iter->gety();
                        x0 = x0 - round(x0/lx)*lx;
                        y0 = y0 - round(y0/ly)*ly;
                        x1 = x1 - round(x1/lx)*lx;
                        y1 = y1 - round(y1/ly)*ly;
                        c0 = (y0 + pow(x0, 2)/y0)/2.0;
                        c1 = (y1 + pow(x1, 2)/y1)/2.0;
                       // cout << x0<< "  "<<y0<<"   "<< distance(cell.begin(), jter)<<endl;

                        xsolve.push_back((c0 - c1)/(x0/y0 - x1/y1)); // xsolve and ysolve record the intersecton
                        ysolve.push_back(-x0/y0*(c0 - c1)/(x0/y0 - x1/y1) + c0);//among all the midperpendicular lines of iter
                        id0.push_back(distance(cell.begin(), jter));// id0 and id1 record the cells' id that their
                        id1.push_back(distance(cell.begin(), kter));//lines belong to. 
                        //cout << x0*y1 - x1 * y0 << "   "<< *(id0.end() - 1) << "  "<< *(id1.end() - 1)<<endl;
                        if(x0*y1 - x1*y0 < 0)
                        {
                            swap(*(id0.end() - 1), *(id1.end() - 1));
							
                        }
                    }
                }
            }
        }
		


    // xsolve and ysove record the coordinates of cross points of all midperpendicular lines
    // In this loop, the right cross points will be filtered out. then put those points into neighbour list.
        length = xsolve.size(); 
        for(int i = 0; i <length; i++) //checking wether the solution is the node of cell or not
        {
            countn = 0;// The nodes of cell should be on the same side with cell center, countn record
                        // how many times the solution are on the same side with center to all midperpendicular lines.
            for(auto jter = cell.begin(); jter != cell.end(); jter++)
            {
                if(jter!=iter)
                {
                    x0 = jter->getx() - iter->getx();
                    y0 = jter->gety() - iter->gety();
                    x0 = x0 - round(x0/lx)*lx;
                    y0 = y0 - round(y0/ly)*ly;
                    
                    //cout << x0<< "   "<<y0<< "   "<<distance(cell.begin(), jter)<<endl;
                    c0 = (y0 + pow(x0, 2)/y0)/2.0;
                    double fxy = ysolve[i] + x0/y0*xsolve[i] - c0;
                    if(c0 * fxy < 1.0e-10)
                        countn++; 
                }
            }
            if(countn == cell.size() - 1) // the the times on the same side with center equatl to num of 
            {                            //midperpendicular lines which is cell.size() - 1, then it is the 
                tempid0.push_back(id0[i]); // one of cell node.
                tempid1.push_back(id1[i]);
             //tempid0 and tempid1 record the id of cell
                tempsolvex.push_back(xsolve[i]);
                tempsolvey.push_back(ysolve[i]);
            } 
        }
        for(int i = 0; i<tempsolvex.size(); i++)
        {
            tempsolvex[i] += iter->getx();
            tempsolvey[i] += iter->gety();
//            tempsolvex[i]  = tempsolvex[i] - round(tempsolvex[i]);
//            tempsolvey[i]  = tempsolvey[i] - round(tempsolvey[i]);
            // periodic boundary condition
           // cout << tempid0[i]<< "   "<< tempid1[i]<< "   "<< tempsolvex[i]<<"  "<< tempsolvey[i]<<endl;
        }
//swap the index of neighboures and nodes, make sure they are in counterclocwise order.
        for(int i = 1; i< tempid1.size() - 1; i++)
        {
            for(int j = i; j< tempid1.size(); j++)
            {
                if(tempid1[i - 1] == tempid0[j])
                {
                    swap(tempid0[i], tempid0[j]);
                    swap(tempid1[i], tempid1[j]);
                    swap(tempsolvex[i], tempsolvex[j]);
                    swap(tempsolvey[i], tempsolvey[j]);
                }
            }
        }
         
        for(int i = 0; i< tempid1.size(); i++)
        {
         //   cout << tempid0[i]<<","<< tempid1[i]<<","<< tempsolvex[i]<<","<< tempsolvey[i]<<endl; 
        }
        C2V_nodex.push_back(tempsolvex);
        C2V_nodey.push_back(tempsolvey); 
    }
	
	for(int i = 0; i< C2V_nodex.size(); i++)
	{
		for(int j = 0; j< C2V_nodex[i].size(); j++)
		{
			C2V_nodex[i][j] = C2V_nodex[i][j] - round(C2V_nodex[i][j]);
			C2V_nodey[i][j] = C2V_nodey[i][j] - round(C2V_nodey[i][j]);
		}
	}

    double temp;
    countn  = 0;
    for(auto iti : C2V_nodex)
    {
        for(auto itj: iti)
        {
            temp = itj - round(itj);
            xsolve_unique.push_back(temp);
        }
    }
    for(auto iti : C2V_nodey)
    {
        for(auto itj: iti)
        {
            temp = itj - round(itj);
            ysolve_unique.push_back(temp);
        }
    }
    length = xsolve_unique.size();
    for(int i = 0; i< length - 1; i++)
    {
        for(int j = i+1; j< length; j++)
        {
            if( abs(xsolve_unique[i] - xsolve_unique[j]) < 1.0e-10)
            {
                xsolve_unique.erase(xsolve_unique.begin() + j);
                ysolve_unique.erase(ysolve_unique.begin() + j);
                length--;
                j--;
            }
        }
    }
    vector<int> temp_id;
    for(auto iti: C2V_nodex)
    {
        temp_id.clear();
        for(auto itj: iti)
        {
            for(auto it = xsolve_unique.begin(); it != xsolve_unique.end(); it++)
            {
                if(abs(*it - itj) < 1.0e-10)
                {
                    temp_id.push_back(distance(xsolve_unique.begin(), it));
                }
            }
        }
        C2V_id.push_back(temp_id);
    }
//=======================get V2C_id=======================================
	GetV2C_id();
//=======================get V2V_id=======================================
	GetV2V_id();
//=======================get L2V_id=======================================
	GetL2V_id();
//========================================================================
}    


// output the cells' coordinates on screen
void vertex_model::get()
{
    cout << "cells' coordinates are: "<<endl;
    cout<< left<<setw(12)<<"x";
    cout << right<<setw(12)<<"y"<<endl;
    for(vector<cells>::iterator iter = cell.begin(); iter != cell.end(); iter++)
    {
        iter->get();
    }
}

//====================get L2V_id===================================
void vertex_model::GetL2V_id()
{
	vector<vector<int>> :: iterator it;
	vector<int> temp_id;
	L2V_id.clear();
	for(int i = 0; i< C2V_id.size(); i++)
	{
		for(int j = 0; j<C2V_id[i].size(); j++)
		{
			temp_id.clear();
			if(C2V_id[i][j] < C2V_id[i][(j+1)%C2V_id[i].size()])
			{
				temp_id.push_back(C2V_id[i][j]);
				temp_id.push_back(C2V_id[i][(j+1)%C2V_id[i].size()]);
			}
			else
			{
				temp_id.push_back(C2V_id[i][(j+1)%C2V_id[i].size()]);
				temp_id.push_back(C2V_id[i][j]);
			}
			it = find(L2V_id.begin(), L2V_id.end(), temp_id);
			if(it == L2V_id.end())
			{
				L2V_id.push_back(temp_id);
			}
		}
	}


}

//====================get V2C_id=====================================

void vertex_model::GetV2C_id()
{
    vector<int> :: iterator it;	
vector<int> temp_id, test{24, 9, 33};
	V2C_id.clear();
	for(int k = 0; k < xsolve_unique.size();k++)
	{
		temp_id.clear();
		for(int i = 0; i < C2V_id.size(); i++)
		{
			for(int j = 0; j < C2V_id[i].size(); j++)
			{
				it = find(temp_id.begin(), temp_id.end(), i);
				if(k==C2V_id[i][j] & it == temp_id.end())
				{
					temp_id.push_back(i);
				}

			}
		}
		V2C_id.push_back(temp_id);
	}
	for(int k = 0; k<xsolve_unique.size(); k++)
	{
		int c0 = V2C_id[k][0];
		int c1 = V2C_id[k][1];
		int c2 = V2C_id[k][2];
		double x0ij, y0ij, xc0ij, yc0ij, xc1ij, yc1ij, xc2ij, yc2ij;
		x0ij = xsolve_unique[k];
		y0ij = ysolve_unique[k];
		xc0ij = cell[c0].getx();
		yc0ij = cell[c0].gety();
		xc1ij = cell[c1].getx();
		yc1ij = cell[c1].gety();
		xc2ij = cell[c2].getx();
		yc2ij = cell[c2].gety();
		if(abs(xc0ij - x0ij)>lx/2)
			xc0ij = xc0ij + abs(x0ij)/x0ij * lx;
		if(abs(yc0ij - y0ij)>ly/2)
			yc0ij = yc0ij + abs(y0ij)/y0ij * ly;
		if(abs(xc1ij - x0ij)>lx/2)
			xc1ij = xc1ij + abs(x0ij)/x0ij * lx;
		if(abs(yc1ij - y0ij) > ly/2)
			yc1ij = yc1ij + abs(y0ij)/y0ij * ly;
		if(abs(xc2ij - x0ij) > lx/2)
			xc2ij = xc2ij + abs(x0ij)/x0ij * lx;
		if(abs(yc2ij - y0ij) > ly/2)
			yc2ij = yc2ij + abs(y0ij)/y0ij * ly;
		double xkc0 = xc0ij - xc2ij;
		double ykc0 = yc0ij - yc2ij;
		double xkc1 = xc1ij - xc2ij;
		double ykc1 = yc1ij - yc2ij;
		
		// xkc0 = xkc0 - round(xkc0/lx)*lx;
		// ykc0 = ykc0 - round(ykc0/ly)*ly;
		// xkc1 = xkc1 - round(xkc1/lx)*lx;
		// ykc1 = ykc1 - round(ykc1/ly)*ly;
		
		/*if(k == 67 )
		{
			cout << xc0ij<< " "<< yc0ij << " "<< xc1ij << " "<< yc1ij<<endl ;
			cout << xkc0 * ykc1 - xkc1 * ykc0<<endl;
			for( auto it1: V2C_id[k])
			{
				cout << it1 << "  ";
			}
			cout << endl;
		}*/

		if(xkc0 * ykc1 - xkc1 * ykc0 < 0)
		{
			swap(V2C_id[k][0],V2C_id[k][1]);

		}
		
	}


}

//=======get V2V_id==============================================
void vertex_model::GetV2V_id()
{
	vector<int> temp_id;
	vector<int> :: iterator it;
	int ii(0);
	V2V_id.clear();
    for(int i = 0; i< xsolve_unique.size(); i++)
    {
        temp_id.clear();
        for(auto iti = C2V_id.begin(); iti != C2V_id.end(); iti++)
        {
            for(auto itj = (*iti).begin(); itj != (*iti).end(); itj++)
            {
                if((*itj) == i)
                {
                    ii = distance((*iti).begin(), itj);
                    ii = (ii - 1 + (*iti).size()) % (*iti).size();
                    it = find(temp_id.begin(), temp_id.end(),*((*iti).begin() + ii));
                    if(it == temp_id.end())
                        temp_id.push_back( *((*iti).begin() +ii) );

                    ii = (ii + 2 + (*iti).size()) % (*iti).size();
                    it = find(temp_id.begin(), temp_id.end(), *((*iti).begin() + ii));
                    if(it == temp_id.end())
                        temp_id.push_back( *((*iti).begin() +ii) );
                }
            }
        }
        V2V_id.push_back(temp_id);
    }
}
