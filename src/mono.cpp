//[[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::vec gaussian1(arma::vec x)
{  
  vec ret = exp(-0.5*pow(x, 2) )/sqrt(2*datum::pi);
  
  return ret;
  
}




// [[Rcpp::export]]
double estimate1(arma::vec TT1,  arma::vec YYYY, arma::vec NN, double h1,  arma::vec weight, double t)
{    mat AA=zeros(3,3);
  mat BB=zeros(3,1);
  int n=NN.n_elem;
  int shu=0;
  vec A1=zeros(n);
  vec B1=zeros(n);
  
  for (int i = 1; i <= n; i++)
  {  vec TT=TT1.subvec(shu,NN(i-1)+shu-1);
    vec YY=YYYY.subvec(shu,NN(i-1)+shu-1);
    vec KK=gaussian1((TT-t)/h1);
    
  A1(i-1)=sum(KK%YY);
  B1(i-1)=sum(KK);
    
    
    shu=shu+NN(i-1); 
   
  
    
  }
  
  double CCC=sum(weight%A1)/sum(weight%B1);
  
  
  return CCC;
}







