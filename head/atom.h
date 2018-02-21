#ifndef atom_h
#define atom_h
#include "parameter.h"
#include <list>
#include <vector>
#include "ndarrays.h"
#include <iostream>
#include <fstream>
#include <string>
class atom{
	public:
		atom()=default;
		atom(double x1,double x2,double r):x(x1),y(x2),radius(r){};
		void setx(double x1);
		void sety(double x2);
		void setr(double ra);
		double getx();
		double gety();
		double getr();
		void setstress_tensor(std::vector<double> a);
		std::vector<double> getstress();
		void printneighbor();
		void printstress();
		void printinfo();
		friend double distance(atom&,atom&);
		friend double potential(atom&,atom&);
		friend void updatelist(ndarrays<atom>&,int);
    friend void updatetensor(ndarrays<atom>&,int);
		friend std::vector<double> str_tensor(atom&,atom&);
		friend std::ostream& operator<<(std::ostream& os,atom& output);
		friend std::fstream& operator<<(std::fstream& fs,atom& output);
	private:
		double x;
		double y;
		double radius;
		std::list<int> neighborx;
		std::list<int> neighbory;
    std::vector<double> stresstensor;
};
int count(ndarrays<atom>& all,atom&,double r,int size);
void print_radial_dis(double,double,ndarrays<atom>&,int,std::string);
#endif
