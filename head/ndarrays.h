#include <stdarg.h>
#include <vector>
#include <new>
#include <iostream>
#ifndef __NDARRAYS_H_INCLUDED__
#define __NDARRAYS_H_INCLUDED__
template<typename T> class ndarrays{
	public:
	//construct differnt dimention list//
		ndarrays(int dimen,...);
		ndarrays(const ndarrays<T>& a);
		ndarrays& operator=(const ndarrays<T>& a);
  	T& operator()(int a,...);
		~ndarrays(){
		delete [] p;
  	}
	private:
		int dimention;
		std::vector<int> size;
		T *p;
		int linearlength;
	};
template<typename T>
ndarrays<T>::ndarrays(int dimen,...){
	/*construct a ndimentional array, the first argument is it's dimention, the comming arguments are the lengths of different dimentions*/
	va_list args;
	std::vector<int> a(0);
	va_start(args,dimen);
	dimention=dimen;
	int temp=1;
	int temp2;
	for(int i=0;i<dimen;i++){
		temp2=va_arg(args,int);
		a.push_back(temp2);
		temp=temp*(temp2);
	}
	va_end(args);
	size=a;
	linearlength=temp;
	p=new T[temp];
}
template<typename T>
ndarrays<T>::ndarrays(const ndarrays<T>& a){
	size=a.size;
	dimention=a.dimention;
	int temp=1;
	for(auto i=size.begin();i<size.end();i++){
		temp=temp*(*i);
	}
	p=new T[temp];
	for(size_t t=0;t<a.linearlength;t++){
		p[t]=a.p[t];
	}
}
template<typename T>
ndarrays<T>& ndarrays<T>::operator=(const ndarrays<T>& a){
	size=a.size;
	dimention=a.dimention;
	int temp=1;
	for(auto i=size.begin();i<size.end();i++){
		temp=temp*(*i);
	}
	p=new T[temp];
	for(size_t t=0;t<a.linearlength;t++){
		p[t]=a.p[t];
	}
	return *this;
}
template<typename T>
T& ndarrays<T>::operator()(int a,...){
	int temp;
	va_list arg;
	int first;
	va_start(arg,a);
	first=a;
	if(dimention ==1){
		return p[first];
	}
	std::vector<int> coord(0);
	coord.push_back(first);
	for(size_t t=1;t<dimention;t++){
		first=va_arg(arg,int);
		coord.push_back(first);
	}
	va_end(arg);
	for(size_t t=0;t<dimention-1;t++){
		if(t==0) temp=coord[t]*size[t+1]+coord[t+1];
		else temp=temp*size[t+1]+coord[t+1];
	}
	return p[temp];
}
#endif
