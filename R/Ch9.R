################
#Chapter 9: Model Calibration
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2026). Experimental Design for Data Science and Engineering. CRC Press.
######################################

#######################################
#Section 9.1: Nonlinear Optimal Designs
#######################################

#################################
#Section 9.1.1: Analytical Models
#################################

setwd("~/R commands/Book/Ch9")
so2=read.csv(file = "so2.csv", header = TRUE)
a=nls(y1~y0+((beta0+beta1*x1)/(beta0+beta1*x1+x2)-y0)*(1-exp(-(beta0+beta1*x1+x2)*dt/m)),start=list(beta0=0,beta1=1,m=1000),data=so2)
summary(a)

h=function(u,eta){
  x=u[1:3]
  val=u[4]+((eta[1]+eta[2]*x[1])/(eta[1]+eta[2]*x[1]+x[2])-u[4])*(1-exp(-(eta[1]+eta[2]*x[1]+x[2])*x[3]/eta[3]))
  return(val)
}
sse=function(eta){
  y=so2[,5]
  yhat=apply(so2[,1:4],1,h,eta=eta)
  val=sum((y-yhat)^2)
  return(val)
}
optim(c(0,1,1000),sse)

library(numDeriv)
n=dim(so2)[1]
S=matrix(0,nrow=n,ncol=3)
for(i in 1:n) {
  hg=function(x) h(u=as.numeric(so2[i,1:4]),eta=x)
  S[i,]=grad(hg,x=c(0,1,1000))
}

gradient=function(x,eta){
  term1=((eta[1]+eta[2]*x[1])/(eta[1]+eta[2]*x[1]+x[2])-x[4])
  term2=exp(-(eta[1]+eta[2]*x[1]+x[2])*x[3]/eta[3])
  v1=term1*term2*x[3]/eta[3]+x[2]/(eta[1]+eta[2]*x[1]+x[2])^2*(1-term2)
  v2=term1*term2*x[1]*x[3]/eta[3]+x[1]*x[2]/(eta[1]+eta[2]*x[1]+x[2])^2*(1-term2)
  v3=-term1*term2*(eta[1]+eta[2]*x[1]+x[2])*x[3]/eta[3]^2
  return(c(v1,v2,v3))
}
S=t(apply(so2[,1:4],1,gradient,eta=c(0,1,1000)))
library(Matrix)
rankMatrix(S)

#An example (two-factor exponential model)

library(ICAOD)
f=function(x){
  eta0+exp(-eta1*x1)+exp(-eta2*x2)
}
set.seed(1)
m=4
opdes=locally(formula=~eta0+exp(-eta1*x1)+exp(-eta2*x2), predvars=c("x1","x2"),parvars=c("eta0","eta1","eta2"),
              lx=c(0,0), ux=c(1,1),inipars=c(0,2,2),iter=100,k=m, family=gaussian(), ICA.control = list("trace"=FALSE))
opdes$design
D=matrix(as.numeric(opdes$design[2:(2*m+1)]),nrow=m,ncol=2)
w=as.numeric(opdes$design[(2*m+2):(2*m+1+m)])
round(D,2)
round(w,3)
opdes$design[2*m+2+m]

#Psuedo Bayesian D-optimal design
library(ICAOD)
pr=uniform(lower=c(-10,1,1), upper = c(10,3,3))
m=4
opdes=bayes(formula=~eta0+exp(-eta1*x1)+exp(-eta2*x2), predvars=c("x1","x2"),parvars=c("eta0","eta1","eta2"),
            lx=c(0,0), ux=c(1,1),prior = pr,iter=100,k=m, family=gaussian(), ICA.control = list("trace"=FALSE))
opdes$design
D=matrix(as.numeric(opdes$design[2:(2*m+1)]),nrow=m,ncol=2)
w=as.numeric(opdes$design[(2*m+2):(2*m+1+m)])
round(D,2)
round(w,3)
opdes$design[2*m+2+m]

###############################################
#Section 9.1.2 Computationally Expensive Models
###############################################

#Figure 9.1

rm(list=ls())
library(rkriging)
library(FNN)
library(ICAOD)
library(numDeriv)
set.seed(1)
setwd("~/R commands/Book/Ch9/functional")
MP=read.csv(file="results.csv", h=TRUE)
q=5
D=MP[-51,1:q]
colnames(D) <- NULL
# transform it back to 0 and 1 scale
# Diffusivity, (1e-12,1e-9), log
D[,1] <- (log(D[,1]) - log(1e-12)) / (log(1e-9) - log(1e-12))
# Surface Concentration, (4e-3,5e-3), original
D[,2] <- (D[,2] - 4e-3) / (5e-3 - 4e-3)
# Accessible Polymer Concentration: (5e-3,6e-3), original
D[,3] <- (D[,3] - 5e-3) / (6e-3 - 5e-3)
# Hindering, (500,2500), original
D[,4] <- (D[,4] - 500) / (2500 - 500)
# Reaction rate, (1e-3,1e1), log
D[,5] <- (log(D[,5]) - log(1e-3)) / (log(1e1) - log(1e-3))


N=dim(D)[1]
q=dim(D)[2]
pdet <- c(seq(0, 5000, length.out = 2001), seq(5000, 120000, length.out = 4601)[-1])
M=length(pdet)
filename=paste(paste("run",sep="",1:50),sep="",".txt")
Y=matrix(0,nrow=N,ncol=M)
for(i in 1:N) Y[i,]=read.table(file=filename[i])$V1

tsamp=seq(0,max(pdet),length=1001)
Ysamp=Y[,1:1001]
for(i in 1:N) {
  Ysamp[i,]=knn.reg(train=cbind(pdet),test=cbind(tsamp),y=Y[i,],k=1)$pred
}

M=length(tsamp)
tsamp01=tsamp/max(tsamp)

r.x=function(u, theta.x)
{
  A=t(t(D1)-u)
  basis.x=function(h) exp(-sum((h/theta.x)^2))
  vec=apply(A,1,basis.x)
  return(vec)
}

r.t=function(v, theta.t)
{
  basis.t=function(h) exp(-abs(h/theta.t))
  vec=basis.t(v-tsamp01)
  return(vec)
}
Rtinv.mat=function(rho)
{
  A=diag(M)
  for(i in 2:(M-1))
  {
    A[i,i-1]=A[i,i+1]=-rho
    A[i,i]=1+rho^2
  }
  A[1,2]=A[M,M-1]=-rho
  return(A/(1-rho^2))
}
ML=function(para)
{
  theta.x=para[1:q]
  theta.t=para[q+1]
  Rx=apply(D1,1,r.x,theta=theta.x)
  Rxinv=solve(Rx+10^(-6)*diag(N))
  rho=exp(-1/((M-1)*theta.t))
  Rtinv=Rtinv.mat(rho)
  a=c(t((Rxinv%*%rep(1,N))%*%(t(rep(1,M))%*%Rtinv)))
  b=c(t(Rxinv%*%Ysamp%*%Rtinv))
  mu=sum(b)/sum(a)
  sigma2=1/(M*N)*sum((c(t(Ysamp))-mu)*(b-mu*a))
  val=M*N*log(sigma2)+M*determinant(Rx,logarithm = TRUE)$mod[1]+M*N*log(1-rho^2)
  return(val)
}

y=apply(Ysamp,1,mean)
D1=as.matrix(D[,1:q])
a=Fit.Kriging(D1,y, kernel.parameters = list(type="Gaussian"))
ini.x=Get.Kriging.Parameters(a)$lengthscale*sqrt(2)

y=apply(Ysamp,2,mean)
a=Fit.Kriging(tsamp/max(tsamp),y, kernel.parameters = list(type="Matern12"))
ini.t=Get.Kriging.Parameters(a)$lengthscale

ini=c(ini.x,ini.t)
a.opt=optim(ini,ML, lower=ini/10, upper=ini*10, method="L-BFGS-B")
theta=a.opt$par
theta

theta.x=theta[1:q]
theta.t=theta[q+1]
Rx=apply(D1,1,r.x,theta=theta.x)
A=solve(Rx+10^(-6)*diag(N))
Rxinv=A
for(i in 1:3)  Rxinv=A%*%(diag(N)+10^(-6)*Rxinv)
rho=exp(-1/((M-1)*theta.t))
Rtinv=Rtinv.mat(rho)
a=c(t((Rxinv%*%rep(1,N))%*%(t(rep(1,M))%*%Rtinv)))
b=c(t(Rxinv%*%Ysamp%*%Rtinv))
mu=sum(b)/sum(a)
sigma2=1/(M*N)*sum((c(t(Ysamp))-mu)*(b-mu*a))
coef=c(b-mu*a)
COEF=matrix(coef,nrow=N,ncol=M,byrow = TRUE)


library(numDeriv)
hhat=function(ti,eta) drop(mu+t(r.x(eta, theta.x)%*%COEF%*%r.t(ti,theta.t)))

fim=function(x,w,param){
  S=matrix(0,nrow=m,ncol = q)
  X=matrix(x,nrow=m)
  for(i in 1:m) {
    heta=function(eta) hhat(X[i,],eta) 
    S[i,]=grad(heta,x=param)
  }
  M=t(S)%*%diag(w)%*%S
  return(M)
}
m=q
opdes=locally(fimfunc = fim,
              lx=c(0), ux=c(1),inipars=rep(.5,q),iter=100,k=m, family=gaussian(), ICA.control = list("trace"=FALSE))
opdes$design
D=matrix(as.numeric(opdes$design[2:(m+1)]),nrow=m,ncol=1)
w=as.numeric(opdes$design[(m+2):(2*m+1)])
round(D,2)
round(w,3)
opdes$design[2*m+2]
det(fim(x=c(D),w=w,param=rep(.5,q)))
DVPI=D
pred=apply(cbind(tsamp01),1,hhat, eta=rep(.5,q))
yD=knn.reg(train=cbind(tsamp01),test=cbind(D),y=pred,k=1)$pred

par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
matplot(pdet,t(Y),"l", xlab="time in seconds", ylab="mass uptake", main="functional output", cex.lab=2,cex.main=3,cex.axis=1)
plot(tsamp,pred,"l", ylim = range(Y), xlab="time in seconds", ylab="mass uptake", main="optimal design", cex.lab=2,cex.main=3,cex.axis=1,lwd=2)
points(cbind(D*120000,yD),pch=1,col="blue",cex=2,lwd=4)
par(mfrow=c(1,1))

###########################
#Section 9.2 Robust designs
###########################

######################################
#Section 9.2.1 Parameter Uncertainties
######################################

#Figure 9.2 (two-factor exponential model)

library(support)
#from https://github.com/BillHuang01/PQMC/blob/master/scripts/lib.R
sp.sample <- function(x, n, prob, tol = 1e-6, iter.max = 10){
  # support points resampling
  if (is.null(dim(x))) stop("x must be a matrix!")
  N <- nrow(x)
  x.dist <- as.matrix(dist(x))
  x.measure <- c(x.dist %*% prob)
  # initial sample
  idx <- c(which.min(x.measure))
  if (prob[idx[1]] < 1/N) idx[1] <- which.max(prob)
  for (i in 2:n){
    measure <- x.measure - apply(matrix(x.dist[,idx],ncol=(i-1)), 1, sum) / i
    idx <- c(idx, which.min(measure))
  }
  edist <- 2 * sum(x.measure[idx]) / n - sum(x.dist[idx,idx]) / n^2
  # iterative update for improvement
  iter <- 0
  while (TRUE){
    iter <- iter + 1
    for (i in 1:n){
      measure <- x.measure - apply(x.dist[,idx[-i]], 1, sum) / n
      idx[i] <- which.min(measure)
    }
    edist.new <- 2 * sum(x.measure[idx]) / n - sum(x.dist[idx,idx]) / n^2
    if ((edist-edist.new) < tol || iter > iter.max){
      break
    } else {
      edist <- edist.new
    }
  }
  return (idx)
}
L=matrix(c(0,0,1/2,0,0,1/2,1/2,1/2),ncol=2,byrow=TRUE)
set.seed(1)
N=20
p=2
samp=sp(n=N, p=2, dist.str=c("uniform","normal"))$sp
samp[,1]=1+2*samp[,1]
samp[,2]=2+1/2*samp[,2]
D=NULL
for(i in 1:20){
  D=rbind(D,matrix(c(0,0,1/samp[i,1],0,0,1/samp[i,2],1/samp[i,1],1/samp[i,2]),ncol=2,byrow=TRUE))
}
ind=sp.sample(D,8,prob=rep(.25/20,80))
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
plot(L,,pch=16,col=4,cex=4,main="Locally D-optimal Design",xlim=c(0,1),ylim=c(0,1),cex.main=2.5, xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=2,cex.axis=2)
text(L+.025,"2",col=2,cex=2)
plot(D,xlim=c(0,1),ylim=c(0,1),main="Robust Design",cex.main=2.5, xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=2,cex.axis=2,lwd=2)
points(D[ind,],pch=16,col=4,cex=4)
text(.025,.025,"2",col=2,cex=2)
par(mfrow=c(1,1))

######################################
#Section 9.2.2 Model Uncertainties
######################################

#Figure 9.3

setwd("~/R commands/Book/Ch9/functional")
best=read.table(file="run51.txt")$V1
pdet <- c(seq(0, 5000, length.out = 2001), seq(5000, 120000, length.out = 4601)[-1])
sorp=read.csv(file="exp_sorption_130c_8.7torr_483nm.csv")
desorp=read.csv(file="exp_desorption_130c_8.7torr_483nm.csv")
correct=35358.89275-34644.0896-1
yexp=c(sorp[,2],desorp[,2]+correct)
texp=c(sorp[,1]^2,desorp[,1]^2+62973)
plot(sqrt(texp),yexp, main="model prediction vs real data", cex.main=2, xlab=expression(sqrt(t)),ylab="mass uptake",lwd=3,cex.lab=2,cex.axis=1.5, type="l", lty=2)
lines(sqrt(pdet),best,col=2,lwd=3)
legend("bottomright", legend=c("measured", "predicted"), col=c(1,2),lty=c(2,1), bty="n",cex=2, lwd=c(3,3))

#Figure 9.4 (two-factor exponential model)

set.seed(1)
N=20
p=2
samp=sp(n=N, p=2, dist.str=c("uniform","normal"))$sp
samp[,1]=1+2*samp[,1]
samp[,2]=2+1/2*samp[,2]
D=NULL
for(i in 1:20){
  D=rbind(D,matrix(c(0,0,1/samp[i,1],0,0,1/samp[i,2],1/samp[i,1],1/samp[i,2]),ncol=2,byrow=TRUE))
}
ind=sp.sample(D,4,prob=rep(.25/20,80))
library(MaxPro)
cand=MaxProLHD(100,2)$Design
a=MaxProAugment(ExistDesign = D[ind,], CandDesign = cand, nNew = 4)
Dnew=a$Design[5:8,]
plot(D,xlim=c(0,1),ylim=c(0,1),main="Robust Design",cex.main=2.5, xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=2,cex.axis=2,lwd=2)
points(D[ind,],pch=16,col=4,cex=4)
points(Dnew,col=2,pch=4,cex=4,lwd=4)

