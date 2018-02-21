#include "atom.h"
#include "parameter.h"
#include <cmath>
#include <string>
#include <vector>
#include "ndarrays.h"
void atom::setx(double x1){
	x=x1;
}
void atom::sety(double x2){
	y=x2;
}
void atom::setr(double ra){
	radius=ra;
}
void atom::setstress_tensor(std::vector<double> a){
        stresstensor=a;
}
std::vector<double> atom::getstress(){
	return stresstensor;
}
double atom::getx(){
	return x;
}
double atom::gety(){
	return y;
}
double atom::getr(){
	return radius;
}
void atom::printneighbor(){
		std::list<int>::iterator b=neighbory.begin();
		std::cout<<"the neighbors are:"<<std::endl;
		for(std::list<int>::iterator a=neighborx.begin();a!=neighborx.end();a++){
			std::cout<<"("<<*a<<","<<*b<<")"<<std::endl;
			b++;
		}
	}
void atom::printstress(){
	std::cout<<"the stress tensor is: (sigma_xx,sigma_xy;sigma_yx,sigma_yy):"<<std::endl;
	std::cout<<stresstensor[0]<<" "<<stresstensor[1]<<std::endl;
	std::cout<<stresstensor[2]<<" "<<stresstensor[3]<<std::endl;
}
void atom::printinfo(){
	std::cout<<"for atom position ("<<x<<","<<y<<")"<<std::endl;
	printneighbor();
	printstress();
}
double distance(atom& one,atom& two){
	double r=(one.x-two.x)*(one.x-two.x)+(one.y-two.y)*(one.y-two.y);
	return sqrt(r);
}
double potential(atom& one,atom& two){
	double r=(one.x-two.x)*(one.x-two.x)+(one.y-two.y)*(one.y-two.y);
	r=sqrt(r);
	if(r<r0){
		return 4*eps*(pow(sigma/r,12)-pow(sigma/r,6));
	}
	else if(r>r_cut){
		return 0;
	}
	else{
		return A*pow(r-r_cut,3)+B*pow(r-r_cut,2);
	}
}
//the force exerted on one by two
std::vector<double> str_tensor(atom& one,atom& two){
	std::vector<double> a(4,0);
	double r=distance(one,two);
    double deri=0;
    if(r<r0){
        deri=24*eps/r*(-2*pow(sigma/r,12)+pow(sigma/r,6));
    }
    else if(r<r_cut){
        deri=3*A*(r-r_cut)*(r-r_cut)+2*B*(r-r_cut);
    }
    else{
        deri=0.0;
    }
    a[0]=deri*(one.x-two.x)*(one.x-two.x)/r;
    a[1]=deri*(one.x-two.x)*(one.y-two.y)/r;
    a[2]=deri*(one.x-two.x)*(one.y-two.y)/r;
    a[3]=deri*(one.y-two.y)*(one.y-two.y)/r;
    return a;
}
void updatelist(ndarrays<atom>& input,int size){
	for(size_t i=0;i<size;i++)
		for(size_t j=0;j<size;j++){
			input(i,j).neighborx.clear();
			input(i,j).neighbory.clear();
		}
	for(size_t i=0;i<size;i++)
		for(size_t j=0;j<size;j++)
			for(size_t k=0;k<size;k++)
				for(size_t l=0;l<size;l++)
				{
					if(distance(input(i,j),input(k,l))<r_cut && distance(input(i,j),input(k,l))>0.0){
						input(i,j).neighborx.push_back(k);
						input(i,j).neighbory.push_back(l);
					}
				}
}
void updatetensor(ndarrays<atom>& atomall,int size){
    std::vector<double> temp(4,0);
    std::vector<double> all(4,0);
    std::list<int>::iterator b;
    for(size_t i=0;i<size;i++)
        for(size_t j=0;j<size;j++){
            for(size_t k=0;k<4;k++){
                all[k]=0.0;
            }
            b=atomall(i,j).neighbory.begin();
            for(std::list<int>::iterator a=atomall(i,j).neighborx.begin();a!=atomall(i,j).neighborx.end();a++){
                temp=str_tensor(atomall(i,j),atomall(*a,*b));
                for(size_t k=0;k<4;k++){
                    all[k]=all[k]+temp[k];
                }
                b++;
            }
            atomall(i,j).setstress_tensor(all);
        }
}
int count(ndarrays<atom>& all,atom& input,double r,int size){
    int c=0;
    for(size_t i=0;i<size;i++)
        for(size_t j=0;j<size;j++){
            if(distance(all(i,j),input)<r){
                c++;
            }
        }
    return c-1;
}
std::ostream& operator<<(std::ostream& os,atom& output){
		os<<"for atom position ("<<output.x<<","<<output.y<<")"<<std::endl;
		os<<"stress tensor(xx,xy,yx,yy): "<<output.stresstensor[0]<<" "<<output.stresstensor[1]<<" "<<output.stresstensor[2]<<" "<<output.stresstensor[3]<<std::endl;
		return os;
}
std::fstream& operator<<(std::fstream& os,atom& output){
		os<<"for atom position "<<output.x<<" "<<output.y<<std::endl;
		os<<"stress tensor(xx,xy,yx,yy): "<<output.stresstensor[0]<<" "<<output.stresstensor[1]<<" "<<output.stresstensor[2]<<" "<<output.stresstensor[3]<<std::endl;
		return os;
}
void print_radial_dis(double r_start,double r_stop,ndarrays<atom>& atomall,int size,std::string name){
	 // double r_start=0.0000001;
  //  double r_stop=15;
    int N=10000;
    std::vector<double> ra_dis(N,0.0);
    double r_delta=(r_stop-r_start)/N;
    double r_inter=0.0;
    int count_old=0;
    int count_new=0;
    int count_delta=0;
    std::fstream radis;
    radis.open(name,std::fstream::out);
    for(size_t i=0;i<N;i++){
        r_inter=i*r_delta+r_start;
        count_new=count(atomall,atomall(10,10),r_inter,size);
        count_delta=count_new-count_old;
        ra_dis[i]=count_delta/2/Pi/r_inter/r_delta;
        count_old=count_new;
        radis<<r_inter<<" "<<ra_dis[i]<<std::endl;
    }
    radis.close();
}
