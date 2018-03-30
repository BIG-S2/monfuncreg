\name{monfuncreg}
\alias{monfuncreg}
\alias{monfuncreg}
\docType{package}
\title{
  Monotone Nonparametric Regression for
Functional/Longitudinal Data
}
\description{
  monfuncreg provides a increasing monotone estimator of the mean regression
function in functional/longitudinal models.
}

\usage{
monfuncreg(x, y, NN, N, hr, hd, weight="OBS", t)      
 }

\arguments{
\item{x}{vector containing the x-values (design points) of a sample, scaled x to [0,1] \eqn{x=(s_{11},...,s_{1m_1},...,s_{nm_n})}.}
\item{y}{vector containing the y-values (response) of a sample, \eqn{y=(y_{11},...,y_{1m_1},...,y_{nm_n})}.}
\item{NN}{the number of observations for each subject, \eqn{NN=(m_1,...,m_n)}.}
\item{N}{the number of evaluation points of the unconstrained
nonparametric regression estimator.}
\item{hr}{bandwith of kernel \eqn{K_r} of the regression estimation step.}
\item{hd}{bandwith of kernel \eqn{K_d} of the density estimation step.}
\item{weight}{"OBS" or "SUBJ".}
\item{t}{ vector of points where the monotone estimation is computed, which is on [0,1].}
}



\value{
  \code{monfuncreg} returns a list of values
  \item{mon1$variable}{the points, for which the monotone function values will be estimated}
  \item{mon1$estimate}{the monotone estimate at mon1$variable}
}


\references{
Monotone Nonparametric Regression for
Functional/Longitudinal Data,  by Chen,  Gao, Fu and Zhu, 2018
}

\details{
Functional/Longitudinal data analysis  has wide application in the biomedical, psychometric and environmental
sciences  (Fitzmaurice et al., 2004; Yao et al., 2005; Wu and Zhang, 2006;
Wang et al., 2016; Zhu et al., 2018). In this type of analysis, subjects are repeatedly
measured over time, and measurements from the same subject are usually highly correlated.
Let \eqn{m_i} be the number of repeated measurements for subject \eqn{i}  and \eqn{n} be the total number of subjects.
 The observations from each subject
are assumed to be noisy discrete realizations of an underlying process \eqn{X(\cdot)}and given by
\deqn{
y_{ij}=X_i(s_{ij})+\sigma(s_{ij})\varepsilon_{ij}~~\mbox{for}~~  j=1,\ldots,m_i;~i=1,\ldots,n,
}
where \eqn{y_{ij}}is the response variable of interest for subject \eqn{i} measured at time \eqn{s_{ij}}, the \eqn{X_i(\cdot)}'s are independent realizations of the underlying process \eqn{\{X(\cdot)\}},  and the \eqn{\varepsilon_{ij}}'s are random errors with zero mean and variance of 1.
By using a mixed effects approach, we decompose \eqn{X_i(s_{ij})} into an unknown population mean \eqn{m(\cdot)=E\{X_i(\cdot)\}} and a subject-specific   trajectory \eqn{\eta_i(\cdot)} with zero mean and covariance function \eqn{\gamma(s,t)=\mbox{cov}\{\eta_i(s),\eta_i(t)\}}. Then,  we can rewrite the model  as
\deqn{
y_{ij}=m(s_{ij})+\eta_{i}(s_{ij})+\sigma(s_{ij})\varepsilon_{ij} ~~\mbox{for}~~  j=1,\ldots,m_i;~i=1,\ldots,n.
}

To obtain an estimator of \eqn{m^{-1}(t)}, we  need    an unconstrained estimator of \eqn{m(t)}, denoted as \eqn{\widehat{m}(s)},  as follows:
\deqn{
\widehat{m}(s)=\frac{\sum_{i=1}^n\omega_i\sum_{j=1}^{m_i}K_{r, h_r}( s_{ij}-s)Y_{ij}}{\sum_{i=1}^n \omega_i
\sum_{j=1}^{m_i}K_{r, h_r}( s_{ij}-s)},
}
where the \eqn{\omega_i}'s are weights satisfying \eqn{\sum_{i=1}^nm_i\omega_i=1}.
 We consider two commonly used weighting schemes,  OBS for equal weight per observation and SUBJ for equal weight per subject (Yao et al., 2005; Li and Hsing, 2010; Kim and Zhao, 2012; Zhang and Wang,
2016). Specifically,
 we set \eqn{\omega_i=1/(\sum_{i=1}^nm_i)}  for OBS, whereas  we set \eqn{\omega_i=1/(nm_i)} for SUBJ.
 Moreover,  \eqn{\widehat{m}(s)} is
a local constant estimator of \eqn{ m(\cdot)}.
 By plugging \eqn{\widehat{m}(s)} , we obtain
\deqn{
 \widehat{m}_I^{-1}(t)={N}^{-1}\int_{-\infty}^t\sum_{i=1}^NK_{d, h_d}( {\widehat{m}({i}/{N})-u})du.
}
  Our constrained estimator of \eqn{m(s)}, denoted as \eqn{\widehat{m}_I(s)},  is then calculated by using a numerical inversion (see also, Dette, et al., 2006).
}







\examples{
# This example is analyzing the real grey matter volume data 
# obtained from the ADNI study using our model, method and package. 
library(monfuncreg)
shuju=read.csv("Origdata1.csv",head=F,quote="")
YY=list()
TT=list()
jishu=0
for(i in 1:length(shuju[,1])){
  if((shuju[i,1]==1)&(is.na(shuju[i,13])==FALSE)){    
    jishu1=1
    jishu=jishu+1    
    YY[[jishu]]=numeric(0)
    TT[[jishu]]=numeric(0)    
    YY[[jishu]][jishu1]=shuju[i,13]
    TT[[jishu]][jishu1]=shuju[i,5]   
  }  
  if(i>1){ 
    if((shuju[i,1]>shuju[i-1,1])&(is.na(shuju[i,13])==FALSE)){    
    jishu1=jishu1+1   
    YY[[jishu]][jishu1]=shuju[i,13]
    TT[[jishu]][jishu1]=shuju[i,5]   
    }    
 }  
}

NN=numeric(0)
for(i in 1:562){ 
  NN[i]=length(YY[[i]]) 
}

TT2=numeric(sum(NN))
YY2=numeric(sum(NN))

shu1=0
for(i in 1:562){  
  TT2[(shu1+1):(shu1+NN[i])]=TT[[i]]
  YY2[(shu1+1):(shu1+NN[i])]=YY[[i]]
  shu1=shu1+NN[i] 
}

TT1=(TT2-min(TT2))/(max(TT2)-min(TT2))
YY1=-(YY2-mean(YY2))/sd(YY2)
n=length(NN)

quant1=quantile(TT1,probs=seq(0,1,0.25))
hr=1.06*(n*mean(NN))^(-1/5)*min(sd(TT1),(quant1[4]-quant1[2])/1.34)
quant1=quantile(YY1,probs=seq(0,1,0.25))
hd=1.06*(n*mean(NN))^(-0.3)*min(sd(YY1),(quant1[4]-quant1[2])/1.34)
t1=seq(0.01,0.99,by=0.001)
N=1000

weight="OBS"
jieguo1=monfuncreg(TT1, YY1, NN,  hr, hd, N, weight,t1)
weight="SUBJ"
jieguo2=monfuncreg(TT1, YY1, NN,  hr, hd, N, weight,t1)

matplot(jieguo1$variable,cbind(-jieguo1$estimate,-jieguo2$estimate),type="l",
lty=1,lwd=2,xlab="Scaled Age", ylab="Standardized Volume of Grey Matter")

}
