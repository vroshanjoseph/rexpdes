################
#Chapter 2: Modeling Techniques
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2026). Experimental Design for Data Science and Engineering. CRC Press.
############
#Section 2.1: Interpolation
############


########################################
#Section 2.1.1 Polynomial Interpolation
########################################
#Figure 2.1
rm(list=ls())
f=function(x) sin(10*pi*x)/(1+64*(x-.25)^2)+x^2
test=seq(0,1,length=301)
true=f(test)

#polynomial interpolation
par(mfrow=c(1,2))
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="polynomial interpolation equi-spaced points",cex.main=1.5,cex.lab=2)
n=10

#equi-spaced design
D=((1:n)-1)/(n-1)
y=f(D)
points(cbind(D,y),pch=16,cex=3,col="blue")

d=2*D-1
X=NULL
for(i in 1:n) X=cbind(X,cos((i-1)*acos(d)))
a=solve(X,y)
a
u=seq(0,1,length=301)
yhat=u
us=2*u-1
for(i in 1:301){
  yhat[i]=0
  for(j in 1:n) yhat[i]=yhat[i]+a[j]*cos((j-1)*acos(us[i]))
}
lines(u,yhat, col=3,lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(2,3),bty="n",cex=2)

#Chebyshev nodes
d=cos((2*(1:n)-1)/(2*n)*pi)
D=(d+1)/2
y=f(D)
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="polynomial interpolation with Chebyshev nodes",cex.main=1.5,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")

X=NULL
for(i in 1:n) X=cbind(X,cos((i-1)*acos(d)))
a=solve(X,y)
a
u=seq(0,1,length=301)
yhat=u
us=2*u-1
for(i in 1:301){
  yhat[i]=0
  for(j in 1:n) yhat[i]=yhat[i]+a[j]*cos((j-1)*acos(us[i]))
}
lines(u,yhat, col=3,lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(2,3),bty="n",cex=2)
par(mfrow=c(1,1))


#Figure 2.2

#two-dimensional
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
D=matrix(c(0,0,0,1,1,0,1,1),nrow=4,ncol=2,byrow=TRUE)
plot(D,pch=16,cex=4,col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.axis=2,cex.lab=2)
D=matrix(c(0,0,.2,.9,.9,.1,1,1),nrow=4,ncol=2,byrow=TRUE)
plot(D,pch=16,cex=4,col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.axis=2,cex.lab=2)
par(mfrow=c(1,1))

######################################
#Section 2.1.2 Radial Basis Functions
#####################################

f=function(x) sin(10*pi*x)/(1+64*(x-.25)^2)+x^2
test=seq(0,1,length=301)
true=f(test)

n=10
D=((1:n)-1)/(n-1)
y=f(D)
par(mar=c(5,5,4,2))
par(mfrow=c(1,3))
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main=expression(theta==0.01),cex.main=3,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")

E=as.matrix(dist(D))
theta=.01
R=exp(-(E/theta)^2)
coef=solve(R,y)
basis=function(h) exp(-(h/theta)^2)
r=function(x) basis(x-D)
fhat=function(x) sum(r(x)*coef)
pred=test
for(i in 1:301) pred[i]=fhat(test[i])
lines(test,pred,col=3, lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(2,3),bty="n",cex=3)

plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main=expression(theta==0.1),cex.main=3,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")
theta=.1
R=exp(-(E/theta)^2)
coef=solve(R,y)
basis=function(h) exp(-(h/theta)^2)
r=function(x) basis(x-D)
fhat=function(x) sum(r(x)*coef)
pred=test
for(i in 1:301) pred[i]=fhat(test[i])
lines(test,pred,col=3, lwd=3)

plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main=expression(theta==1),cex.main=3,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")
theta=1
R=exp(-(E/theta)^2)
coef=solve(R,y)
basis=function(h) exp(-(h/theta)^2)
r=function(x) basis(x-D)
fhat=function(x) sum(r(x)*coef)
pred=test
for(i in 1:301) pred[i]=fhat(test[i])
lines(test,pred,col=3, lwd=3)
par(mfrow=c(1,1))


n=10
D=((1:n)-1)/(n-1)
y=f(D)
E=as.matrix(dist(D))
MSCV=function(theta)
{
  R=exp(-(E/theta)^2)
  Rinv=solve(R)
  cv=c(Rinv%*%y)/diag(Rinv)
  val=log(mean(cv^2))
  return(val)
}

ev=seq(.01,.2,length=100)
val=ev
for(i in 1:100) val[i]=MSCV(ev[i])
plot(ev,val,"l")

a.opt=optimize(MSCV,c(.01,.2))
theta=a.opt$min
theta


######################
#Section 2.1.3 Kriging
######################

#Figure 2.4
rm(list=ls())
f=function(x) sin(10*pi*x)/(1+64*(x-.25)^2)+x^2
test=seq(0,1,length=301)
true=f(test)
n=10
D=((1:n)-1)/(n-1)
y=f(D)
par(mar=c(5,5,4,2))

plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="ordinary kriging",cex.main=3,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")

E=as.matrix(dist(D))
theta=.08
R=exp(-(E/theta)^2)
basis=function(h) exp(-(h/theta)^2)
r=function(x) basis(x-D)
one=rep(1,n)
Rinv=solve(R+10^(-10)*diag(n))
mu=drop(t(one)%*%Rinv%*%y/(t(one)%*%Rinv%*%one))
mu
coef=Rinv%*%(y-mu)

tau2=drop(1/n*t(y-mu)%*%Rinv%*%(y-mu))
f.var=function(x)
{ 
  fac=1/(t(one)%*%Rinv%*%one)
  val=tau2*(1-t(r(x))%*%Rinv%*%r(x)+fac*(1-t(r(x))%*%Rinv%*%one)^2)
  return(val)
}

fhat=function(x) mu+t(r(x))%*%coef
N=301
pred=low=up=numeric(N)
for(i in 1:N)
{
  pred[i]=fhat(test[i])
  SD=sqrt(f.var(test[i]))
  low[i]=pred[i]-2*SD
  up[i]=pred[i]+2*SD
}
polygon(c(test,rev(test)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(test,low,lwd=2, col="pink")
lines(test,up,lwd=2, col="pink")
lines(test,pred,col=3,lwd=4,lty=2)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(test,true,lty=2,lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(3,4),bty="n",cex=2)

#########################################
#Section 2.1.4 Gaussian Process Modeling
#########################################

#Figure 2.5

rm(list=ls())
f=function(x) sin(10*pi*x)/(1+64*(x-.25)^2)+x^2
N=301
test=seq(0,1,length=N)
true=f(test)
n=10
D=((1:n)-1)/(n-1)
y=f(D)

library(rkriging)
kernel=Get.Kernel(.08/sqrt(2),type="Gaussian")

a=Fit.Kriging(D,y,fit=FALSE,kernel = kernel)
tau2=Get.Kriging.Parameters(a)$nu2
pred=Predict.Kriging(a,test)

set.seed(1)
library(MASS)
S0=Evaluate.Kernel(kernel,test)
u=mvrnorm(n=5,mu=rep(0,N),Sigma = tau2*S0)
par(mfrow=c(1,3))
matplot(test,t(u),type="l",xlab="x",ylab="y",col=1:5,lty=1,lwd=3, ylim=c(min(true)-.25,max(true)+.25),main="prior distribution",cex.main=3,cex.axis=2,cex.lab=2)

matplot(test,t(u),type="l",xlab="x",ylab="y",col=1:5,lty=1,lwd=3, ylim=c(min(true)-.25,max(true)+.25),main="data",cex.main=3,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=4,col="blue")


R=Evaluate.Kernel(kernel,D)
L=chol(R+10^(-10)*diag(n))
library(pdist)
M=t(as.matrix(pdist(cbind(D),cbind(test))))
M=exp(-M^2/.08^2)
b=forwardsolve(t(L),t(M))
S=tau2*(S0-t(b)%*%b)
library(Matrix)
S=nearPD(S)$mat
me=pred$mean
set.seed(1)
u=mvrnorm(n=5,mu=me,Sigma = S)
matplot(test,t(u),type="l",xlab="x",ylab="y",col=1:5,lty=1,lwd=3, ylim=c(min(true)-.25,max(true)+.25),main="posterior distribution",cex.main=3,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=4,col="blue")
par(mfrow=c(1,1))


#Gaussian process prediction (not shown in the book)
a<-Fit.Kriging(D,y,kernel.parameters = list(type="Gaussian"))
Get.Kriging.Parameters(a)
theta=Get.Kriging.Parameters(a)$lengthscale*sqrt(2)
theta
pred=Predict.Kriging(a,test)
low=pred$mean-2*pred$sd
up=pred$mean+2*pred$sd
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="Gaussian process",cex.main=3,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")
polygon(c(test,rev(test)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(test,low,lwd=2, col="pink")
lines(test,up,lwd=2, col="pink")
lines(test,pred$mean,col=3,lwd=4,lty=2)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(test,true,lty=2,lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(3,4),bty="n",cex=2)


##########################################################################################
#Section 2.2 Regression
###########################################################################################

################################
#Section 2.2.1 Linear Regression
################################
#Figure 2.6
rm(list=ls())
f=function(x) sin(10*pi*x)/(1+64*(x-.25)^2)+x^2
test=seq(0,1,length=301)
true=f(test)

set.seed(5)
n=10
r=2

#equi-spaced design
D=((1:n)-1)/(n-1)
D=rep(D,r)
e=rnorm(n*r,sd=.1)
y=f(D)+e

par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
i=3
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="3rd degree polynomial regression (equi-spaced)",cex.main=1.5,cex.lab=2)
d=2*D-1
a=lm(y~poly(d,i))
u=seq(0,1,length=301)
us=2*u-1
pred=predict(a,newdata = data.frame(d=us), interval="confidence")
low=pred[,2]
up=pred[,3]
lines(u,pred[,1], col=3,lwd=3)
polygon(c(u,rev(u)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(u,low,lwd=2, col="pink")
lines(u,up,lwd=2, col="pink")
lines(u,pred[,1],col=3,lwd=4,lty=1)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(u,f(u),lty=2,lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(2,3),bty="n",cex=2)

i=8
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="8th degree polynomial regression (equi-spaced)",cex.main=1.5,cex.lab=2)
d=2*D-1
a=lm(y~poly(d,i))
u=seq(0,1,length=301)
us=2*u-1
pred=predict(a,newdata = data.frame(d=us), interval="confidence")
low=pred[,2]
up=pred[,3]
lines(u,pred[,1], col=3,lwd=3)
polygon(c(u,rev(u)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(u,low,lwd=2, col="pink")
lines(u,up,lwd=2, col="pink")
lines(u,pred[,1],col=3,lwd=4,lty=1)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(u,f(u),lty=2,lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(2,3),bty="n",cex=2)
par(mfrow=c(1,1))


#Figure 2.7

library(AlgDesign)
cand=data.frame(x=seq(-1,1,length=301))
a=optFederov(~ 1 + x + I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6)+I(x^7)+I(x^8)+I(x^9), nTrials=n,data=cand, approximate=FALSE)
D=rep(c(a$design[,1]),r)
D=(D+1)/2
y=f(D)+e


par(mfrow=c(1,2))
i=3
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="3rd degree polynomial regression (D-optimal)",cex.main=1.5,cex.lab=2)
d=2*D-1
a=lm(y~poly(d,i))
u=seq(0,1,length=301)
us=2*u-1
pred=predict(a,newdata = data.frame(d=us), interval="confidence")
low=pred[,2]
up=pred[,3]
lines(u,pred[,1], col=3,lwd=3)
polygon(c(u,rev(u)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(u,low,lwd=2, col="pink")
lines(u,up,lwd=2, col="pink")
lines(u,pred[,1],col=3,lwd=4,lty=1)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(u,f(u),lty=2,lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(2,3),bty="n",cex=2)

i=8
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="8th degree polynomial regression (D-optimal)",cex.main=1.5,cex.lab=2)
d=2*D-1
a=lm(y~poly(d,i))
u=seq(0,1,length=301)
us=2*u-1
pred=predict(a,newdata = data.frame(d=us), interval="confidence")
low=pred[,2]
up=pred[,3]
lines(u,pred[,1], col=3,lwd=3)
polygon(c(u,rev(u)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(u,low,lwd=2, col="pink")
lines(u,up,lwd=2, col="pink")
lines(u,pred[,1],col=3,lwd=4,lty=1)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(u,f(u),lty=2,lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(2,3),bty="n",cex=2)
par(mfrow=c(1,1))


#######################################
#Section 2.2.2 Kernel ridge regression
#######################################

#Figure 2.8

rm(list=ls())
f=function(x) sin(10*pi*x)/(1+64*(x-.25)^2)+x^2
test=seq(0,1,length=301)
true=f(test)

set.seed(5)
n=10
re=2

#equi-spaced design
D=((1:n)-1)/(n-1)
D=rep(D,re)
e=rnorm(n*re,sd=.1)
y=f(D)+e
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
E=as.matrix(dist(D))
MSCV=function(para)
{
  theta=para[1]
  lam=para[2]
  R=exp(-(E/theta)^2)
  Rinv=solve(R+lam*diag(n*re))
  cv=c(Rinv%*%y)/diag(Rinv)
  val=log(mean(cv^2))
  return(val)
}

a=optim(c(.1,.1),MSCV,lower=c(.05,.001),upper=c(1,1),method = "L-BFGS-B")
theta=a$par[1]
theta
lam=a$par[2]
lam
R=exp(-(E/theta)^2)
coef=solve(R+lam*diag(n*re),y)
basis=function(h) exp(-(h/theta)^2)
r=function(x) basis(x-D)
fhat=function(x) sum(r(x)*coef)
pred=test
for(i in 1:301) pred[i]=fhat(test[i])
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="kernel ridge regression (equi-spaced)",cex.main=1.5,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(test,pred,col=3, lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(2,3),bty="n",cex=2)

GCV=function(para)
{
  theta=para[1]
  lam=para[2]
  R=exp(-(E/theta)^2)
  Rinv=solve(R+lam*diag(n*re))
  cv=c(Rinv%*%y)/mean(diag(Rinv))
  val=log(mean(cv^2))
  return(val)
}
a=optim(c(.1,.1),GCV,lower=c(.05,.001),upper=c(1,1),method = "L-BFGS-B")
theta=a$par[1]
theta
lam=a$par[2]
lam
R=exp(-(E/theta)^2)
coef=solve(R+lam*diag(n*re),y)
basis=function(h) exp(-(h/theta)^2)
r=function(x) basis(x-D)
fhat=function(x) sum(r(x)*coef)
pred=test
for(i in 1:301) pred[i]=fhat(test[i])

#D-optimal design
library(AlgDesign)
cand=data.frame(x=seq(-1,1,length=301))
a=optFederov(~ 1 + x + I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6)+I(x^7)+I(x^8)+I(x^9), nTrials=n,data=cand, approximate=FALSE)
D=rep(c(a$design[,1]),re)
D=(D+1)/2
y=f(D)+e
E=as.matrix(dist(D))
a=optim(c(.1,.1),MSCV,lower=c(.05,.001),upper=c(1,1),method = "L-BFGS-B")
theta=a$par[1]
lam=a$par[2]
R=exp(-(E/theta)^2)
coef=solve(R+lam*diag(n*re),y)
basis=function(h) exp(-(h/theta)^2)
r=function(x) basis(x-D)
fhat=function(x) sum(r(x)*coef)
pred=test
for(i in 1:301) pred[i]=fhat(test[i])
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="kernel ridge regression (D-optimal)",cex.main=1.5,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(test,pred,col=3, lwd=3)
par(mfrow=c(1,1))

##########################################
#Section 2.2.3 Gaussian process Regression
##########################################

#Figure 2.9

library(rkriging)
rm(list=ls())
f=function(x) sin(10*pi*x)/(1+64*(x-.25)^2)+x^2
test=seq(0,1,length=301)
true=f(test)

set.seed(5)
n=10
re=2

#equi-spaced design
D=((1:n)-1)/(n-1)
D=rep(D,re)
e=rnorm(n*re,sd=.1)
y=f(D)+e
a=Fit.Kriging(D,y,interpolation=FALSE, kernel.parameters = list(type="Gaussian"))
Get.Kriging.Parameters(a)
pred=Predict.Kriging(a,test)
low=pred$mean-2*pred$sd
up=pred$mean+2*pred$sd
par(mfrow=c(1,2))
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="Gaussian process regression (equi-spaced)",cex.main=1.5,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")
polygon(c(test,rev(test)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(test,low,lwd=2, col="pink")
lines(test,up,lwd=2, col="pink")
lines(test,pred$mean,col=3,lwd=4,lty=1)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(test,true,lty=2,lwd=3)
legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),lwd=c(3,4),bty="n",cex=2)

#D-optimal design
library(AlgDesign)
cand=data.frame(x=seq(-1,1,length=301))
a=optFederov(~ 1 + x + I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6)+I(x^7)+I(x^8)+I(x^9), nTrials=n,data=cand, approximate=FALSE)
D=rep(c(a$design[,1]),re)
D=(D+1)/2
D=rep(D,re)
e=rnorm(n*re,sd=.1)
y=f(D)+e
a=Fit.Kriging(D,y,interpolation=FALSE, kernel.parameters = list(type="Gaussian"))
Get.Kriging.Parameters(a)
pred=Predict.Kriging(a,test)
low=pred$mean-2*pred$sd
up=pred$mean+2*pred$sd
plot(test,true,type="l",xlab="x",ylab="y",lwd=3,lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="Gaussian process regression (D-optimal)",cex.main=1.5,cex.axis=2,cex.lab=2)
points(cbind(D,y),pch=16,cex=3,col="blue")
polygon(c(test,rev(test)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(test,low,lwd=2, col="pink")
lines(test,up,lwd=2, col="pink")
lines(test,pred$mean,col=3,lwd=4,lty=1)
points(cbind(D,y),pch=16,cex=3,col="blue")
lines(test,true,lty=2,lwd=3)
par(mfrow=c(1,1))
