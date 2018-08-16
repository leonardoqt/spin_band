#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

int main()
{
	string tmp, label_states, label_kpoints, label_u, label_d, label_electron, label_END;
	stringstream ss; 
	int t1,t2,t3,nband, nkpoint, fnband, flag, fillbd, nele, nhsp;
	double ** band_u, ** band_d, *kx, *ky, *kz, *x, add, maxe,mine,fermi;
	double dkx1,dky1,dkz1,dkx2,dky2,dkz2,dk1,dk2,kdot, hsp[30];

	label_states = "number of Kohn-Sham states";
	label_kpoints = "number of k points";
	label_u = "SPIN UP";
	label_d = "SPIN DOWN";
	label_electron = "number of electrons";
	label_END = "JOB DONE";
	getline(cin,tmp);
	while(tmp.find(label_electron) == string::npos)
		getline(cin,tmp);
	ss << (tmp);
	getline(ss,tmp,'=');
	ss >> nele;
	ss.str("");
	ss.clear();
	
	fillbd = (nele+1)/2-1;
	
	while(tmp.find(label_states) == string::npos)
		getline(cin,tmp);
	ss << (tmp);
	getline(ss,tmp,'=');
	ss >> nband;
	ss.str("");
	ss.clear();
	
	while(tmp.find(label_kpoints) == string::npos)
		getline(cin,tmp);
	ss << (tmp);
	getline(ss,tmp,'=');
	ss >> nkpoint;
	ss.str("");
	ss.clear();

	fnband = nband;
	
	band_u = new double * [fnband];
	band_d = new double * [fnband];

	kx = new double [nkpoint];
	ky = new double [nkpoint];
	kz = new double [nkpoint];
	x  = new double [nkpoint];
	
	for (t1 = 0; t1 < fnband; t1++)
	{
		band_u[t1] = new double [nkpoint];
		band_d[t1] = new double [nkpoint];
	}
//get k points
	getline(cin,tmp);
	for (t1 = 0; t1 < nkpoint; t1++)
	{
		getline(cin,tmp,'(');
		getline(cin,tmp,'(');
		getline(cin,tmp,')');
		ss << (tmp);
		ss >> kx[t1]>>ky[t1]>>kz[t1];
		ss.str("");
		ss.clear();
		getline(cin,tmp);

		if(tmp.find(label_END) != string::npos)
		{
			cout<<"can not find all k points"<<endl;
			return 0;
		}
	}
//-------------

	while(tmp.find(label_u) == string::npos)
		getline(cin,tmp);
	getline(cin,tmp);
	for (t1 = 0; t1 < nkpoint; t1++)
	{
		getline(cin,tmp);
		getline(cin,tmp);
		getline(cin,tmp);
		for (t2 = 0; t2 < fnband; t2++)
		{
			cin>>band_u[t2][t1];
		}
		getline(cin,tmp);

	}
	
	while(tmp.find(label_d) == string::npos)
		getline(cin,tmp);
	getline(cin,tmp);
	for (t1 = 0; t1 < nkpoint; t1++)
	{
		getline(cin,tmp);
		getline(cin,tmp);
		getline(cin,tmp);
		for (t2 = 0; t2 < fnband; t2++)
		{
			cin>>band_d[t2][t1];
		}
		getline(cin,tmp);
	}

	x[0] = 0;
	for (t1 = 1; t1 < nkpoint; t1++)
	{
		add = sqrt((kx[t1]-kx[t1-1])*(kx[t1]-kx[t1-1]) + (ky[t1]-ky[t1-1])*(ky[t1]-ky[t1-1]) + (kz[t1]-kz[t1-1])*(kz[t1]-kz[t1-1]));
		if (add < 0.2)
			x[t1] = x[t1-1] + add;
		else
			x[t1] = x[t1-1];
	}
//find high symmetry point
	nhsp=0;
	dkx1 = kx[1], dky1 = ky[1], dkz1 = kz[1];
	dk1 = dkx1*dkx1+dky1*dky1+dkz1*dkz1;
	for (t1 = 2; t1 < nkpoint; t1++)
	{
		dkx2 = kx[t1]-kx[t1-1], dky2 = ky[t1]-ky[t1-1], dkz2 = kz[t1]-kz[t1-1];
		dk2 = dkx2*dkx2+dky2*dky2+dkz2*dkz2;
		kdot = dkx1*dkx2+dky1*dky2+dkz1*dkz2;
		if (kdot/sqrt(dk1*dk2) < 0.95 )
		{
			hsp[nhsp] = x[t1-1];
			nhsp++;
		}
		dkx1=dkx2,dky1=dky2,dkz1=dkz2,dk1=dk2;
	}

//find fermi level?
	for(maxe=-10000,mine=10000,t1=0; t1 < nkpoint; t1++)
	{
		if(band_u[fillbd][t1] > maxe) maxe = band_u[fillbd][t1];
		if(band_d[fillbd][t1] > maxe) maxe = band_d[fillbd][t1];
		if(band_u[fillbd+1][t1] < mine) mine = band_u[fillbd+1][t1];
		if(band_d[fillbd+1][t1] < mine) mine = band_d[fillbd+1][t1];
	}
	fermi = (maxe+mine)/2;
	cout<<"# Number of electrons = "<<nele<<endl;
	cout<<"# Number of bands = "<<fnband<<endl;
	cout<<"# VBM = "<<maxe<<'\t'<<"CBM = "<<mine<<endl;
	cout<<"# Fermi Level (approx.) = "<<fermi<<endl;
	cout<<"# Band gap = "<<mine-maxe<<endl<<endl;
	cout<<"# k\t\tup\t\tdown"<<endl;


	for (t1 = 0; t1 < fnband; t1++)
	{
		flag = t1%2;
		for (t2 = 0; t2 < nkpoint; t2++)
		{
			cout.precision(7);
			cout<<scientific<<x[(1-flag)*t2 + flag*(nkpoint-1-t2)]<<'\t'<<band_u[t1][(1-flag)*t2 + flag*(nkpoint-1-t2)]-fermi<<'\t'<<band_d[t1][(1-flag)*t2 + flag*(nkpoint-1-t2)]-fermi<<endl;
		}
	}
	for (t1 = 0; t1 < nhsp ; t1++)
		cout<<endl<<hsp[t1]<<"\t-100"<<endl<<hsp[t1]<<"\t100"<<endl;

	return 0;
}
