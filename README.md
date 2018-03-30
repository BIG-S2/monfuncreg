# monfuncreg
Monotone Nonparametric Regression for Functional/Longitudinal Data

monfuncreg provides a increasing monotone estimator of the mean regression function in functional/longitudinal models.

# Table of contents
1. [The Package method](#The Package method)
    1. [Model](#Model)
    2. [Estimating method](#method)
    3. [Data analysis](#Data analysis)
    4. [Code for the data analysis](#Code)


# The package method <a name="The packageod"></a>
## Model <a name="Model"></a>
 In functional/longitudinal data analysis, subjects are repeatedly
measured over time, and measurements from the same subject are usually highly correlated.
Let $n_i$ be the number of repeated measurements for subject $i$  and $n$ be the total number of subjects.
 The observations from each subject
are assumed to be noisy discrete realizations of an underlying process $\{X(\cdot)\}$ and given by
\begin{eqnarray}\label{I1}
y_{ij}=X_i(s_{ij})+\sigma(s_{ij})\varepsilon_{ij}  j=1,\ldots,n_i;~i=1,\ldots,n,
\end{eqnarray}
where $y_{ij}$ is the response variable of interest for subject $i$ measured at time $s_{ij}$, the $X_i(\cdot)$'s are independent realizations of the underlying process $\{X(\cdot)\}$,  and the $\varepsilon_{ij}$'s are random errors with zero mean and variance of 1.
By using a mixed effects approach, we decompose $X_i(s_{ij})$ into an unknown population mean $m(\cdot)=E\{X_i(\cdot)\}$ and a subject-specific   trajectory $\eta_i(\cdot)$ with zero mean and covariance function $\gamma(s,t)=\mbox{cov}\{\eta_i(s),\eta_i(t)\}$. Then,  we can rewrite (\ref{I1})  as
\begin{eqnarray}\label{I2}
y_{ij}=m(s_{ij})+\eta_{i}(s_{ij})+\sigma(s_{ij})\varepsilon_{ij} ~~\mbox{for}~~  j=1,\ldots,n_i;~i=1,\ldots,n.
\end{eqnarray}

## Estimating method  <a name="method"></a>
We define an unconstrained estimator of $m(t)$, denoted as $\widehat{m}(s)$,  as follows:
\begin{eqnarray}
\widehat{m}(s)=\frac{\sum_{i=1}^n\omega_i\sum_{j=1}^{n_i}K_{r, h_r}( s_{ij}-s)Y_{ij}}{\sum_{i=1}^n \omega_i
\sum_{j=1}^{n_i}K_{r, h_r}( s_{ij}-s)},\label{initial}\nonumber
\end{eqnarray}
where the $\omega_i$'s are weights satisfying $\sum_{i=1}^nn_i\omega_i=1$.
 We consider two commonly used weighting schemes,  OBS for equal weight per observation and SUBJ for equal weight per subject  (Yao et al., 2005; Li and Hsing, 2010; Kim and Zhao, 2012; Zhang and Wang,
2016). Specifically,
 we set $\omega_i=1/(\sum_{i=1}^nn_i)$  for OBS, whereas  we set $\omega_i=1/(nn_i)$ for SUBJ.
 Moreover,  $\widehat{m}(s)$ is
a local constant estimator of $m(\cdot)$  (Kim and Zhao, 2012; Zhang and Wang,
2016).
 By plugging $\widehat{m}(s)$, we obtain
 $$
 \widehat{m}_I^{-1}(t)={N}^{-1}\int_{-\infty}^t\sum_{i=1}^NK_{d, h_d}( {\widehat{m}({i}/{N})-u})du. \nonumber
$$
  Our constrained estimator of $m(s)$, denoted as $\widehat{m}_I(s)$,  is then calculated by using a numerical inversion.


# Data analysis <a name="Data analysis"></a>
Our work is motivated by an analysis of structural brain Magnetic Resonance Imaging
(MRI) data extracted from the Alzheimer's Disease Neuroimaging Initiative (ADNI). The
ADNI is an ongoing public-private partnership to test whether genetic, structural and
functional neuroimaging and clinical data can be combined to measure the progression of
mild cognitive impairment and early Alzheimer's disease. Subjects in the ADNI have been
recruited from over 50 sites across the United States and Canada. Our problem of interest is
to study the effect of aging on the progression of Alzheimer's disease. Aging can generally be
referred to as a progressive deterioration of physiological function, leading to impairments
in cognitive function and the ability to execute and learn new movements. In addition to
devastating cognitive impairment, disorders of degenerative dementia such as Alzheimer's
disease are characterized by accelerating cerebral atrophy. Recent methodological advances
in MRI allow for the characterization of the structural changes that accompany aging of the
healthy brain, such as changes in the volume of grey matter (GM). MRI is often used to
differentiate normal aging from the neurodegeneration that is seen with early Alzheimer's
disease. Imaging studies have shown that cerebellar GM volumes are less in elderly people
than in younger adults (Henkenius et al., 2003; Hoogendam et al., 2012), which suggests a
monotonic relationship between GM volume and age.

##Code for the data analysis <a name="Code"></a>
 # This example is analyzing the real grey matter volume data  obtained from the ADNI study using our model, method and package. 

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

matplot(jieguo1\$variable, cbind(-jieguo1\$estimate, -jieguo2\$estimate), type="l", lty=1,lwd=2,xlab="Scaled Age", ylab="Standardized Volume of Grey Matter")

