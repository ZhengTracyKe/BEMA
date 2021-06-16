library(MASS)
library(RMTstat)
library(rmatio)
library(pracma)

SQM <- function(alpha){
  S <- t(data) %*% data  /n
  l <- eigen(S)$values
  gamma=min(p,n)/max(p,n)
  k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  predictor=qmp(k/min(p,n),max(n,p),min(n,p))*max(p,n)/n
  
  sigma2hat=lm(rev(l[k])~predictor-1)$coef[[1]]
  plot(l)
  k0=min(p,n):1
  lines(qmp(k0/min(p,n),max(n,p),min(n,p))*max(p,n)/n*sigma2hat,col='red')
  cutoff=(sigma2hat*((1+sqrt(gamma))^2+qtw(0.9)*max(p,n)^(-2/3)*gamma^(-1/6)*(1+sqrt(gamma))^(4/3)))*max(p,n)/n
  K=sum(l>cutoff)
  return(K)
}

# demo
n=500
p=100
K_true=5
B=qr.Q(qr(mvrnorm(p,rep(0,K_true),diag(K_true))))
Sigma=diag(p)+5*B%*%t(B)
data=mvrnorm(n,rep(0,p),Sigma)
SQM(0.1)
