################
#Chapter 10: Data Subsampling
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2026). Experimental Design for Data Science and Engineering. CRC Press.
######################################

##############################################
#Section 10.1 Support Points-based Subsampling
##############################################

#Figure 10.1

rm(list=ls())
library(support)
library(SPlit)
library(MASS)
set.seed(123)
N=1000
p=2
rho=.5
Sigma=diag(p)
Sigma=rho^abs(row(Sigma)-col(Sigma))
D=mvrnorm(n=N, mu=rep(0,p), Sigma=Sigma)

n=100
S=sp(n,p,dist.samp = D)$sp

par(mfrow=c(1,2))
plot(D, pch=4,cex=.75,lwd=1,main="Support Points", cex.main=2, xlab="x",ylab="y",cex.lab=2,cex.axis=2)
points(S,pch=1,col=2,lwd=3, cex=1.5)
ind=subsample(D,S)

plot(D, pch=4,cex=.75,lwd=1,main="Subsample", cex.main=2, xlab="x",ylab="y",cex.lab=2,cex.axis=2)
points(D[ind,],pch=1,col=2,lwd=3,cex=1.5)
par(mfrow=c(1,1))

#Figure 10.2 (wind energy data)
library(rkriging)
setwd("~/R commands/Book/Ch10")
wind=read.csv("wind10.csv", header=TRUE)
set.seed(123)
p=2
D=as.matrix(wind[,-1])
n=392
system.time({S=sp(n,p,dist.samp = D)$sp
ind=subsample(D,S)
})
S=D[ind,]
system.time({a=Fit.Kriging(X=cbind(S[,1]),y=S[,2],kernel.parameters = list(type="Gaussian"),interpolation = FALSE)})

ev=seq(min(D[,1]),max(D[,1]),length=1000)
pred=Predict.Kriging(a,cbind(ev))$mean
plot(D,pch=4,cex=.5, main="wind power data", xlab="speed", ylab="power", cex.main=2,cex.lab=2, cex.axis=2)
points(S,col=2,lwd=2)
lines(ev,pred,col=4,lwd=3)
legend("bottomright",legend=c("data","subsample","GP"),bty="n",pch=c(4,1,NA),lty=c(NA,NA,1),col=c(1,2,4), cex=2, lwd=c(1,2,3))

############################
#Section 10.2 Data Splitting
############################

#Figure 10.3

library(SPlit)
rm(list=ls())
set.seed(102)
N=100
x1=scale(rnorm(N))
x2=scale(x1^2+rnorm(N))
D=matrix(cbind(x1,x2),ncol=2)
Ntest=20
y=rep(2,N)
y[rbinom(N,1,pnorm(-(6*x1+x2+5)))==1]=1
y[rbinom(N,1,pnorm(-(-6*x1+x2+5)))==1]=3

Dy=data.frame(x1=D[,1],x2=D[,2],y=as.factor(y))
ind.opt=SPlit(Dy,splitRatio = .2)
ind=sample(1:N,Ntest)

par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
plot(x1,x2,pch=y+1,col=y+1,cex.lab=2, main="Random",xlab=expression(x[1]), ylab=expression(x[2]),cex.main=3, cex.axis=2, cex=2,lwd=2)
points(D[ind,],col=1,pch=1,lwd=3,cex=4)
Dtest=D[ind.opt,1:2]
plot(x1,x2,pch=y+1,col=y+1,cex.lab=2, main="SPlit",xlab=expression(x[1]), ylab=expression(x[2]),cex.main=3, cex.axis=2, cex=2,lwd=2)
points(Dtest,col=1,pch=1,lwd=3,cex=4)
par(mfrow=c(1,1))

###########################
#Section 10.3 Data Twinning
###########################

#Figure 10.4 (only the final iteration is shown here)

library(twinning)
set.seed(6)
N = 50
x = scale(rnorm(N))
y = scale(x^2 + rnorm(N)) 
D = cbind(x, y)
ind.twin=twin(D,r=5)
plot(D,pch=16,cex=.75, xlab="x", ylab="y")
points(D[ind.twin,], col=2,cex=2,lwd=2)




#Figure 10.5

rm(list=ls())
set.seed(4)
library(twinning)
library(SPlit)
library(BalancedSampling)

N=1000
X=matrix(rnorm(2*N),ncol=2)
n=100

edsplit=edtwin=edlpm1=edlpm2=numeric(30)
sbsplit=sbtwin=sblpm1=sblpm2=numeric(30)

for(i in 1:30){
  
  ind=SPlit(X,splitRatio = n/N)
  edsplit[i]=energy(X,X[ind,])
  sbsplit[i]=sb(rep(n/N,N), X, ind)
  
  ind=twin(X, r=N/n)
  edtwin[i]=energy(X,X[ind,])
  sbtwin[i]=sb(rep(n/N,N), X, ind)
  
  ind=lpm1(rep(n/N,N),X)
  edlpm1[i]=energy(X,X[ind,])
  sblpm1[i]=sb(rep(n/N,N), X, ind)
  
  ind=lpm2(rep(n/N,N),X)
  edlpm2[i]=energy(X,X[ind,])
  sblpm2[i]=sb(rep(n/N,N), X, ind)
}
par(mfrow=c(1,2))
ed=data.frame(SPlit=edsplit,twin=edtwin,LPM1=edlpm1,LPM2=edlpm2)
boxplot(ed,main="Energy Distance", xlab="sampling methods", ylab="energy distance",cex.main=3, cex.axis=1.5, cex=2, cex.lab=2, col=c("yellow",2:4))
sbmeasure=data.frame(SPlit=sbsplit,twin=sbtwin,LPM1=sblpm1,LPM2=sblpm2)
boxplot(sbmeasure,main="Spatial Balance", xlab="sampling methods", ylab="spatial balance measure",cex.main=3, cex.axis=1.5, cex=2, cex.lab=2, col=c("yellow",2:4))
par(mfrow=c(1,1))

####################################
#Section 10.4 Supervised Compression
####################################

#Figure 10.6

rm(list=ls())
library(supercompress)
library(deldir)
library(FNN)
library(fields)
f=function(x) {
  p=length(x)
  x=pi*x
  val=-sum(sin(x)*(sin((1:p)*x^2/pi))^(20))
  return(val)
}
p=2
N=10000*p

set.seed(1)

x=NULL
for(i in 1:p) x=cbind(x,runif(N))
y=apply(x,1,f)+.01*rnorm(N)
n=200
true=apply(x,1,f)

N.plot=250
p1=seq(0,1,length=N.plot)
p2=seq(0,1,length=N.plot)
fc=matrix(apply(expand.grid(p1,p2),1,f),nrow = N.plot, ncol= N.plot)

system.time({
  a0=kmeans(x,n,iter.max = 100)})
D0=a0$centers
cluster=a0$cluster
y0=sapply(split(y,cluster),mean)
pred.nn=knn.reg(D0,test = x,y0,k=1)$pred
sqrt(mean((pred.nn-true)^2))

system.time({a=supercompress(n,x,y)})
Dopt=a$D
ybar=a$ybar
pred.nn=knn.reg(Dopt,test = x,ybar,k=1)$pred
sqrt(mean((pred.nn-true)^2))

par(mar=c(5,5,4,6))
par(mfrow=c(1,2))
imagePlot(p1,p2,fc,col=cm.colors(5),xlab=expression(x[1]),ylab=expression(x[2]),main="unsupervised",cex.main=3,cex.lab=2, cex.axis=2)
points(D0,pch=16,col=1,cex=2)
vor=deldir(D0[,1],D0[,2])
plot.deldir(vor,add=TRUE,wlines = "tess")
imagePlot(p1,p2,fc,col=cm.colors(5),xlab=expression(x[1]),ylab=expression(x[2]),main="supervised",cex.main=3,cex.lab=2, cex.axis=2)
points(Dopt,pch=16,col=1,cex=2)
vor=deldir(Dopt[,1],Dopt[,2])
plot.deldir(vor,add=TRUE,wlines = "tess")
par(mfrow=c(1,1))


