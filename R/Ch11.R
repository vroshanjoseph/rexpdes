################
#Chapter 11: Data Analysis
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2026). Experimental Design for Data Science and Engineering. CRC Press.
######################################

##############################################
#Section 11.1 Factor Selection and Ranking
##############################################

#Figure 11.1

rm(list=ls())
par(mfrow=c(1,2))
curve(sqrt(2)*x,from=0, to=1, lwd=3, main="Effects", cex.main=2, xlab="x",ylab="y",cex.lab=2,cex.axis=2, col=2)
curve(2/sqrt(3)*x,from=0, to=1,add=TRUE,col=3, lty=1, lwd=3)
curve(1*x,from=0, to=1,add=TRUE,col=4, lty=1, lwd=3)
text(.8,1.22, expression(x[1]),cex=2, col=2)
text(.8,1, expression(x[2]),cex=2, col=3)
text(.8,.7, expression(x[3]),cex=2, col=4)

curve(2*(1-x^2), from=-1,to=1, lwd=3, main="Total Sobol' ", cex.main=2, xlab=expression(rho),ylab=expression(V^tot),cex.lab=2,cex.axis=2,col=2)
curve(4/3*(1-x^2), from=-1,to=1, add=TRUE, col=3,lty=1, lwd=3)
abline(h=1,col=4,lty=1, lwd=3)
text(0,1.9, expression(x[1]),cex=2, col=2)
text(0,1.4, expression(x[2]),cex=2, col=3)
text(0,.9, expression(x[3]),cex=2, col=4)
par(mfrow=c(1,1))

library(MASS)
library(first)
set.seed(123)
N=1000
p=2
rho=.9
Sigma=diag(p)
Sigma=rho^abs(row(Sigma)-col(Sigma))
D=mvrnorm(n=N, mu=rep(0,p), Sigma=Sigma)
D=cbind(D,rnorm(N))
y=3*D[,1]+2*D[,2]+sqrt(3)*D[,3]
first(X=D,y=y)
first(X=D[,-2],y=y)
first.rank(X=D,y=y, noise=FALSE)

#Figure 11.2 (concrete data)--This will take time!

library(SPlit)
library(randomForest)
library(rkriging)
library(first)
setwd("~/R commands/Book/Ch11")
concrete=read.csv(file="Concrete_Data.csv", fileEncoding="UTF-8-BOM")
dim(concrete)
names(concrete)
concrete.scale=data.frame(scale(concrete))

set.seed(123)
rmse.RF=rmse.GP=matrix(0,nrow=30,ncol=2)
for(i in 1:30){
ind=SPlit(concrete, splitRatio = 0.1)
train=concrete.scale[-ind,]
test=concrete.scale[ind,]
#all the factors
a=randomForest(x=train[,1:8],y=train[,9])
pred=predict(a,test[,1:8])
rmse.RF[i,1]=sqrt(mean((test$y-pred)^2))

a=Fit.Kriging(X=train[,1:8], y=train[,9], interpolation = FALSE, kernel.parameters = list(type="Gaussian"))
pred=Predict.Kriging(a,test[,1:8])$mean
rmse.GP[i,1]=sqrt(mean((test$y-pred)^2))

#selected factors
sel=first(X=train[,-9],y=train[,9])
sel=(sel>0)
a=randomForest(x=train[,-9][,sel],y=train[,9])
pred=predict(a,test[,-9][,sel])
rmse.RF[i,2]=sqrt(mean((test$y-pred)^2))

a=Fit.Kriging(X=train[,-9][,sel], y=train[,9], interpolation = FALSE, kernel.parameters = list(type="Gaussian"))
pred=Predict.Kriging(a,test[,-9][,sel])$mean
rmse.GP[i,2]=sqrt(mean((test$y-pred)^2))
}
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
boxplot(data.frame(full=rmse.RF[,1], select=rmse.RF[,2]), ylim=c(0,.37), main="Random Forest", ylab="RMSE", cex.main=2,cex.lab=2,cex.axis=2, col=c(2,3))
boxplot(data.frame(full=rmse.GP[,1], select=rmse.GP[,2]), ylim=c(0,.37), main="Gaussian process", ylab="RMSE", cex.main=2,cex.lab=2,cex.axis=2, col=c(2,3))
par(mfrow=c(1,1))

###################################
#Section 11.2 Twin-Gaussian Process
###################################

#Figure 11.3

library(twingp)
setwd("~/R commands/Book/Ch11")
wind=read.csv("wind10.csv", header=TRUE)
D=as.matrix(wind[,-1])
setwd("~/R commands/Book/Ch11")
wind=read.csv("wind10.csv", header=TRUE)
D=as.matrix(wind[,-1])
ev=seq(min(D[,1]),max(D[,1]),length=1000)
set.seed(123)
system.time({a.twin=twingp(x=cbind(D[,1]),y=D[,2],x_test = cbind(ev))})
pred.twin=a.twin$mu
par(mar=c(5,5,4,2))
plot(D,pch=4,cex=.5, main="TwinGP", xlab="speed", ylab="power", cex.main=3,cex.lab=2,cex.axis=2)
lines(ev,pred.twin,col=4,lwd=3)

