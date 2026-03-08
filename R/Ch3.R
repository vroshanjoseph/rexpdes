################
#Chapter 3: Model-based Designs
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2026). Experimental Design for Data Science and Engineering. CRC Press.
######################################

#######################################
#Section 3.1: Prediction-based Designs
######################################

rm(list=ls())
r=function(x) exp(-((x-D)/theta)^2)
s=function(x){
  a=forwardsolve(t(L),r(x))
  return(sqrt(1-sum(a^2)))
}
imse=function(D){
  E=as.matrix(dist(D))
  R=exp(-(E/theta)^2)
  L=chol(R+10^(-6)*diag(n))
  r=function(x) exp(-((x-D)/theta)^2)
  s=function(x){
    a=forwardsolve(t(L),r(x))
    return(1-sum(a^2))
  }
  val=integrate(Vectorize(s),0,1,abs.tol=10^(-10))$val
  return(val)
}
mmse=function(D){
  E=as.matrix(dist(D))
  R=exp(-(E/theta)^2)
  L=chol(R+10^(-6)*diag(n))
  r=function(x) exp(-((x-D)/theta)^2)
  s=function(x){
    a=forwardsolve(t(L),r(x))
    return(1-sum(a^2))
  }
  test=seq(0,1,length=1001)
  val=apply(cbind(test),1,s)
  return(max(val))
}

#Figure 3.1

p=2
n=7
par(mar=c(5,5,4,2))
par(mfrow=c(1,3))
ind=c(9,45,31)
for(i in 1:3){
  set.seed(ind[i])
  D=matrix(runif(n*p),nrow=n,ncol=p)
  plot(D,xlim=c(0,1),ylim=c(0,1), xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=3,cex.axis=2.5,pch=16,col="blue",cex=4)
}

#Figure 3.2

par(mfrow=c(2,2))
n=10
D=D1=((1:n)-1)/(n-1)
theta=.1
E=as.matrix(dist(D1))
R=exp(-(E/theta)^2)
L=chol(R+10^(-6)*diag(n))
N=301
test=seq(0,1,length=N)
sa=apply(cbind(test),1,s)
plot(test,sa,"l",xlab="x",ylab="RMSE",ylim=c(0,.7),main="equi-spaced design",cex.main=2.5,cex.lab=2,cex.axis=2,lwd=3)
points(cbind(D1,0),pch=16,col="blue",cex=3)
imse_value=imse(D1)
text(.2,.7,paste0("IMSE=",round(imse_value,4)),cex=2,col=2)
mmse_value=mmse(D1)
text(.8,.7,paste0("MMSE=",round(mmse_value,4)),cex=2,col=3)
round(sort(D1),3)

d=cos((2*(1:n)-1)/(2*n)*pi)
D=D2=(d+1)/2
E=as.matrix(dist(D2))
R=exp(-(E/theta)^2)
L=chol(R+10^(-6)*diag(n))
sa=apply(cbind(test),1,s)
plot(test,sa,"l",xlab="x",ylab="RMSE",ylim=c(0,.7),main="Chebyshev design",cex.main=2.5,cex.lab=2,cex.axis=2,lwd=3)
points(cbind(D2,0),pch=16,col="blue",cex=3)
imse_value=imse(D2)
text(.2,.7,paste0("IMSE=",round(imse_value,4)),cex=2,col=2)
mmse_value=mmse(D2)
text(.8,.7,paste0("MMSE=",round(mmse_value,4)),cex=2,col=3)
round(sort(D2),3)

D0=((1:n)-.5)/n
a=optim(D0,imse,lower=rep(0,n),upper=rep(1,n),method="L-BFGS-B")
D=D3=a$par
E=as.matrix(dist(D3))
R=exp(-(E/theta)^2)
L=chol(R+10^(-6)*diag(n))
sa=apply(cbind(test),1,s)
plot(test,sa,"l",xlab="x",ylab="RMSE",ylim=c(0,.7), main="IMSE-optimal design",cex.main=2.5,cex.lab=2,cex.axis=2,lwd=3)
points(cbind(D3,0),pch=16,col="blue",cex=3)
imse_value=imse(D3)
text(.2,.7,paste0("IMSE=",round(imse_value,4)),cex=2,col=2)
mmse_value=mmse(D3)
text(.8,.7,paste0("MMSE=",round(mmse_value,4)),cex=2,col=3)
round(sort(D3),3)

a=optim(D0,mmse,lower=rep(0,n),upper=rep(1,n),method="L-BFGS-B")
D=D4=a$par
E=as.matrix(dist(D4))
R=exp(-(E/theta)^2)
L=chol(R+10^(-6)*diag(n))
sa=apply(cbind(test),1,s)
plot(test,sa,"l",xlab="x",ylab="RMSE",ylim=c(0,.7), main="MMSE-optimal design",cex.main=2.5,cex.lab=2,cex.axis=2,lwd=3)
points(cbind(D4,0),pch=16,col="blue",cex=3)
imse_value=imse(D4)
text(.2,.7,paste0("IMSE=",round(imse_value,4)),cex=2,col=2)
mmse_value=mmse(D4)
text(.8,.7,paste0("MMSE=",round(mmse_value,4)),cex=2,col=3)
round(sort(D4),3)
par(mfrow=c(1,1))

#Figure 3.3

ev=seq(.01,.4,length=100)
IR=MR=matrix(0,nrow=100,ncol=4)
for(i in 1:100){
  theta=ev[i]
  IR[i,]=c(imse(D1),imse(D2),imse(D3),imse(D4))
  MR[i,]=c(mmse(D1),mmse(D2),mmse(D3),mmse(D4))
}
par(mar=c(5,5,4,2))
par(mfrow=c(2,2))
matplot(ev,IR,"l",xlab=expression(theta), ylab="IMSE",lwd=4, main="IMSE",cex.main=2.5,cex.lab=2,cex.axis=2)
legend("topright",legend=c("equi-spaced","Chebyshev","IMSE-optimal","MMSE-optimal"),lty=1:4,col=1:4,lwd=4,bty = "n",cex=2)
matplot(ev,MR,"l",xlab=expression(theta), ylab="MMSE",lwd=4, main="MMSE",cex.main=2.5,cex.lab=2,cex.axis=2)
IRR=sweep(IR,1,apply(IR,1,min),"/")
matplot(ev,1/IRR,"l",xlab=expression(theta), ylab="RE",lwd=4, main="Relative efficiency-IMSE",cex.main=2.5,cex.lab=2,cex.axis=2)
MRR=sweep(MR,1,apply(MR,1,min),"/")
matplot(ev,1/MRR,"l",xlab=expression(theta), ylab="RE",lwd=4, main="Relative efficiency-MMSE",cex.main=2.5,cex.lab=2,cex.axis=2)
par(mfrow=c(1,1))

####################################
#Section 3.2 Maximum Entropy designs
####################################

#Figure 3.4

r=function(x) exp(-((x-D)/theta)^2)
s=function(x){
  a=forwardsolve(t(L),r(x))
  return(sqrt(max(1-sum(a^2),0)))
}
ent=function(D){
  E=as.matrix(dist(D))
  R=exp(-(E/theta)^2)
  val=det(R)
  return(val)
}
negent=function(D) -ent(D)

set.seed(1)
par(mfrow=c(1,3))
n=10

D0=((1:n)-.5)/n
a=optim(D0,imse,lower=rep(0,n),upper=rep(1,n),method="L-BFGS-B")
D=D3=a$par
E=as.matrix(dist(D3))
theta=.1
R=exp(-(E/theta)^2)
L=chol(R+10^(-6)*diag(n))
sa=apply(cbind(test),1,s)
plot(test,sa,"l",xlab="x",ylab="RMSE",ylim=c(0,.7), main="IMSE-optimal design",cex.main=4,cex.lab=2,cex.axis=2,lwd=3)
points(cbind(D3,0),pch=16,col="blue",cex=4)
ent_value=ent(D3)
text(.5,.7,paste0("Entropy=",round(ent_value,4)),cex=3,col=4)
imse_value=imse(D3)
text(.5,.65,paste0("IMSE=",round(imse_value,4)),cex=3,col=2)
mmse_value=mmse(D3)
text(.5,.6,paste0("MMSE=",round(mmse_value,4)),cex=3,col=3)
round(sort(D3),3)


a=optim(D0,mmse,lower=rep(0,n),upper=rep(1,n),method="L-BFGS-B")
D=D4=a$par
E=as.matrix(dist(D4))
R=exp(-(E/theta)^2)
L=chol(R+10^(-6)*diag(n))
sa=apply(cbind(test),1,s)
plot(test,sa,"l",xlab="x",ylab="RMSE",ylim=c(0,.7), main="MMSE-optimal design",cex.main=4,cex.lab=2,cex.axis=2,lwd=3)
points(cbind(D4,0),pch=16,col="blue",cex=4)
ent_value=ent(D4)
text(.5,.7,paste0("Entropy=",round(ent_value,4)),cex=3,col=4)
imse_value=imse(D4)
text(.5,.65,paste0("IMSE=",round(imse_value,4)),cex=3,col=2)
mmse_value=mmse(D4)
text(.5,.6,paste0("MMSE=",round(mmse_value,4)),cex=3,col=3)
round(sort(D4),3)


D0=((1:n)-1)/(n-1)
a=optim(D0,negent,lower=rep(0,n),upper=rep(1,n),method="L-BFGS-B")
D=D5=a$par
E=as.matrix(dist(D5))
R=exp(-(E/theta)^2)
L=chol(R+10^(-6)*diag(n))
sa=apply(cbind(test),1,s)
plot(test,sa,"l",xlab="x",ylab="RMSE",ylim=c(0,.7), main="Maximum entropy design",cex.main=4,cex.lab=2,cex.axis=2,lwd=3)
points(cbind(D5,0),pch=16,col="blue",cex=4)
ent_value=ent(D5)
text(.5,.7,paste0("Entropy=",round(ent_value,4)),cex=3,col=4)
imse_value=imse(D5)
text(.5,.65,paste0("IMSE=",round(imse_value,4)),cex=3,col=2)
mmse_value=mmse(D5)
text(.5,.6,paste0("MMSE=",round(mmse_value,4)),cex=3,col=3)
round(sort(D5),3)
par(mfrow=c(1,1))


