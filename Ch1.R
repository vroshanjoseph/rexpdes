################
#Chapter 1: Introduction
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2025). Experimental Design for Data Science and Engineering. CRC Press.
############
#Section 1.2
############

#Figure 1.4: different objectives of experiment

par(mar=c(5,5,4,2))
par(mfrow=c(2,2))
x=seq(0,1,length=301)
y=sin(10*pi*x)/(1+64*(x-.5)^2)
plot(x,y,type="l",axes=FALSE,xlab="",ylab="",lwd=3,col=1,ylim=c(min(y)-.05,max(y)),main="approximate a rough function",cex.main=3)
n=15
D=((1:n)-.5)/n
points(cbind(D,min(y)-.05),pch=16,cex=2,lwd=3,col="blue")

y=1/(1+64*(x-.5)^2)
plot(x,y,type="l",axes=FALSE,xlab="",ylab="",lwd=3,col=1,ylim=c(min(y)-.05,max(y)),main="approximate a smooth function",cex.main=3)
D=qbeta(((1:n)-.5)/n,.5,.5)
points(cbind(D,min(y)-.05),pch=16,cex=2,lwd=3,col="blue")

n1=n*2/3
n2=n/3
y=2/3*dbeta(x,2,10)+1/3*dbeta(x,10,2)
plot(x,y,type="l",axes=FALSE,xlab="",ylab="",lwd=3,col=1,ylim=c(min(y)-.05,max(y)),main="optimization",cex.main=3)
D0=((1:n2)-.5)/n2
D=c(D0,qbeta(((1:n1)-.5)/n1,2,10))
points(cbind(D,min(y)-.05),pch=16,cex=2,lwd=3,col="blue")

library(support)
y=2/3*dbeta(x,2,10)+1/3*dbeta(x,10,2)
plot(x,y,type="l",axes=FALSE,xlab="",ylab="",lwd=3,col=1,ylim=c(min(y)-.05,max(y)),main="uncertainty propagation",cex.main=3)
x=rbeta(1000*n1,2,10)
x=c(x,rbeta(1000*n2,10,2))
a=sp(n,1,dist.samp = cbind(x))
D=a$sp
points(cbind(D,min(y)-.05),pch=16,cex=2,lwd=3,col="blue")
par(mfrow=c(1,1))

############
#Section 1.3
############

#Figure 1.5: linear regression vs GP

set.seed(7)
n=6
x=rep(c(0,1),c(n/2,n/2))
f=function(x){
  x=(x+.2)/1.4
  val=2/3*dbeta(x,2,10)+1/3*dbeta(x,10,2)
  return(val)
}
f=function(x) sin(10*pi*x)/(1+64*(x-.25)^2)+x^2

y=f(x)+.1*rnorm(n)


par(mfrow=c(1,2))
l=-.25
u=1.25
test=seq(l,u,.001)
truey=f(test)
plot(test,truey,type="l",xlim=c(l,u),ylim=c(-.5,1.5),main="Linear regression",xlab="x",ylab="y",cex.main=2,cex.lab=2)
points(x,y)
a.lm=lm(y~x)

pred=predict.lm(a.lm,newdata = data.frame(x=test),interval = "prediction")
lines(test,pred[,1],col=3,lwd=3,lty=2)
low=pred[,2]
up=pred[,3]
polygon(c(test,rev(test)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(test,low,lwd=2, col="pink")
lines(test,up,lwd=2, col="pink")
lines(test,pred[,1],col=3,lwd=3,lty=2)
lines(test,truey,xlim=c(l,u),lwd=2)
points(x,y,pch=16,col="blue",cex=2)
legend("topleft", legend = c("truth","prediction"),lty=c(1,2),col=c(1,3),lwd=c(2,3),bty = "n",cex = 2)


library(rkriging)
a.gp=Fit.Kriging(x,y,interpolation=FALSE,fit=TRUE,model="OK",kernel.parameters=list(type="Gaussian"))
pred <- Predict.Kriging(a.gp,test)
plot(test,truey,type="l",xlim=c(l,u),ylim=c(-.5,1.5),main="Gaussian process regression",xlab="x",ylab="y",cex.main=2,cex.lab=2)
lines(test,pred$mean,col=3,lwd=3,lty=2)
low=pred$mean-2*pred$sd
up=pred$mean+2*pred$sd
polygon(c(test,rev(test)),c(low,rev(up)),col="pink", border="NA",lwd=2)
lines(test,low,lwd=2, col="pink")
lines(test,up,lwd=2, col="pink")
lines(test,pred$mean,col=3,lwd=3,lty=2)
lines(test,truey,xlim=c(l,u),lwd=2)
points(x,y,pch=16,col="blue",cex=2)
legend("topleft", legend = c("truth","prediction"),lty=c(1,2),col=c(1,3),lwd=c(2,3),bty = "n",cex = 2)
par(mfrow=c(1,1))


