library(MASS)
library(RMTstat)
library(rmatio)
library(pracma)

loss <- function(proposal,l,alpha){
  L=matrix(data=0,nrow=10,ncol=p)
  for (i in 1:10){
    #wi=rWishart(n, p, diag(rgamma(p,proposal,1)))/n
    x1=mvrnorm(n,rep(0,p),diag(rgamma(p,proposal,1)))
    if (p<=n){
      l1=eigen(t(x1)%*%x1)$value
      L[i,]=l1/n
    }
    if (p>n){
      l1=eigen(x1%*%t(x1))$value
      L[i,1:n]=l1/n
    }  
    
  }
  l1=colMeans(L)
  
  k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  s1=lm(l[k]~l1[k]-1)$coef[[1]]
  l1=s1*l1
  return(sum((l1-l)[k]^2))
}

SQM_Gamma <- function(alpha){
  S <- t(data) %*% data  /n
  l <- eigen(S)$values
  
  o=optimize(loss,interval=c(0.1, 50),l=l,alpha=alpha)
  a=o$minimum
  L=matrix(data=0,nrow=500,ncol=p)
  for (i in 1:500){
    x1=mvrnorm(n,rep(0,p),diag(rgamma(p,a,1)))
    if (p<=n){
      l1=eigen(t(x1)%*%x1)$value
      L[i,]=l1/n
    }
    if (p>n){
      l1=eigen(x1%*%t(x1))$value
      L[i,1:n]=l1/n
    }  
  }
  
  l1=c()
  l2=c()
  l3=c()
  l4=c()
  for (i in 1:p){
    l1[i]=mean(L[,i])
    l2[i]=quantile(L[,i],0.9)
    l3[i]=quantile(L[,i],0.1)
    l4[i]=quantile(L[,i],0.5)
  }
  k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  s1=lm(l[k]~l1[k]-1)$coef[[1]]
  plot(l)
  lines(l1*s1,col='red')
  lines(l2*s1,col='orange')
  lines(l3*s1,col='blue')
  lines(l4*s1,col='black')
  abline(max(l2*s1),0,col='orange')
  return(sum(l>max(l2*s1)))
}

# demo
n=500
p=100
K_true=5
B=qr.Q(qr(mvrnorm(p,rep(0,K_true),diag(K_true))))
Sigma=diag(p)+5*B%*%t(B)
data=mvrnorm(n,rep(0,p),Sigma)
SQM_Gamma(0.1)