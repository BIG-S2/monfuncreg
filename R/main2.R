gaussian3=function(x)
{  
  ret = exp(-0.5*x^2 )/sqrt(2*pi)
  ret
  
}

estimate3=function(TT1,   YYYY, NN,  h1,  h2,   weight,  N,  u)
{   
  b=0
  for (i in 1:N)
  {   a=estimate1(TT1,  YYYY, NN, h1,  weight, i/N);
  b=gaussian3((a-u)/h2)/N/h2+b
  }
  
  b
}





integrand=function(TT11,  YY11, N2, h1, h2,  weight11, N1){
  
  ff=function(u){
    
    estimate3(TT11,  YY11, N2, h1, h2,  weight11, N1, u)
    
  }
  
  return(ff)
  
  
}











#' This function estimate of monotone nonparametric function at given points for functional data
#' @param x, scaled x to [0,1],x=(x_{11},...,x_{1m_1},...,x_{1m_n})
#' @param y, y=(y_{11},...,y_{1m_1},...,y_{1m_n})
#' @param NN=(m_1,...,m_n)
#' @param hr is the bandwidth parameter for unconstraint regression
#' @param hd is the bandwidth parameter for constraint regression for inversion function
#' @param N is the N for constraint regression for inversion function
#' @param weight="OBs" or "SUBJ"
#' @param t1 are the points, for which the monotone function values will be estimated
#' @export
#' @examples
#' shuju=read.csv("Origdata1.csv",head=F,quote="")
#' YY=list()
#' TT=list()
#' jishu=0
#' for(i in 1:length(shuju[,1])){
#'  if((shuju[i,1]==1)&(is.na(shuju[i,13])==FALSE)){
#'   jishu1=1
#'   jishu=jishu+1
#'    YY[[jishu]]=numeric(0)
#'    TT[[jishu]]=numeric(0)
#'    YY[[jishu]][jishu1]=shuju[i,13]
#'   TT[[jishu]][jishu1]=shuju[i,5]
#'  }
#'  if(i>1){ if((shuju[i,1]>shuju[i-1,1])&(is.na(shuju[i,13])==FALSE)){
#'    jishu1=jishu1+1
#'   YY[[jishu]][jishu1]=shuju[i,13]
#'    TT[[jishu]][jishu1]=shuju[i,5]
#'  }
#'  }
#' }
#' NN=numeric(0)
#' for(i in 1:562){
#'  NN[i]=length(YY[[i]])
#' }
#' TT2=numeric(sum(NN))
#' YY2=numeric(sum(NN))
#' shu1=0
#' for(i in 1:562){
#'  TT2[(shu1+1):(shu1+NN[i])]=TT[[i]]
#'  YY2[(shu1+1):(shu1+NN[i])]=YY[[i]]
#'  shu1=shu1+NN[i]
#' }
#' TT1=(TT2-min(TT2))/(max(TT2)-min(TT2))
#' YY1=-(YY2-mean(YY2))/sd(YY2)
#' n=length(NN)
#' quant1=quantile(TT1,probs=seq(0,1,0.25))
#' hr=1.06*(n*mean(NN))^(-1/5)*min(sd(TT1),(quant1[4]-quant1[2])/1.34)
#' quant1=quantile(YY1,probs=seq(0,1,0.25))
#' hd=1.06*(n*mean(NN))^(-0.3)*min(sd(YY1),(quant1[4]-quant1[2])/1.34)
#' t1=seq(0.01,0.99,by=0.001)
#' N=1000
#' weight="OBS"
#' jieguo1=monfuncreg(TT1, YY1, NN,  hr, hd, N, weight,t1)
#' weight="SUBJ"
#' jieguo2=monfuncreg(TT1, YY1, NN,  hr, hd, N, weight,t1)
#' matplot(jieguo1$variable,cbind(-jieguo1$estimate,-jieguo2$estimate),type="l",lty=1,lwd=2,xlab="Scaled Age", ylab="Standardized Volume of Grey Matter", main="(B)") 


monfuncreg=function(x, y, NN,  hr, hd, N, weight,t1){
  
   n=length(NN)
   weight2=numeric(n)
  
  
  if(weight=="OBS"){
    
    weight2=rep(1/sum(NN),n) 
    
  }
  
  
  if(weight=="SUBJ"){
    
    weight2=1/n/NN 
    
  }
  
   qq=t1
   weight1=rep(1/sum(NN),n)
   jieguoNW1=numeric(0)
  
   for(i in 1:length(qq)){
    t=qq[i]
    a1=estimate1(x,y,NN,hr,weight1,t) 
    jieguoNW1[i]=a1
    
   }
  
  quant1=quantile(y,probs=seq(0,1,0.25))  
  qqq1=seq(min(jieguoNW1)-hd, max(jieguoNW1)+hd,by=0.002)
  qqq=seq(min(jieguoNW1)-hd, max(jieguoNW1)+hd,length.out = max(length(t1),length(qqq1))+100)

  ggg2=numeric(0)
  ggg2[1]=integrate(integrand(x,y,NN,hr,hd,weight2,N), lower=-Inf, upper= qqq[1],stop.on.error=TRUE)$value
  
   for(i in 2:length(qqq)){
    
    ggg2[i]=integrate(integrand(x,y,NN,hr,hd,weight2,N), lower=qqq[i-1], upper= qqq[i],stop.on.error=TRUE)$value+ggg2[i-1]
    
   }
  
   jieguo1=numeric(0)
   
   for(i in 1:length(qq)){
    t=qq[i]
    a3=qqq[which.min(abs(ggg2-t))] 
    jieguo1[i]=a3
    
   }
  
   mon1=list()
   mon1$variable=t1
   mon1$estimate=jieguo1
  
   return(mon1)
  
  
}


#  return 1: t1
#  return 2: the function values of variable


