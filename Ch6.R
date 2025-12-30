################
#Chapter 6: Screening Designs
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2026). Experimental Design for Data Science and Engineering. CRC Press.
######################################

#######################################
#Section 6.1: Sensitivity Analysis
######################################

##########################################
#6.1.1 Variance-based Sensitivity Measures
##########################################

#Figure 6.1 (left)

rm(list=ls())
#borehole function
f=function(x) {
  lower<-c(0.05,100,63070,990,63.1,700,1120,9855)
  upper<-c(0.15,50000,115600,1110,116,820,1680,12045)
  x=lower+x*(upper-lower)
  val=2*pi*x[3]*(x[4]-x[6])/(log(x[2]/x[1])*(1+2*x[7]*x[3]/(log(x[2]/x[1])*x[1]^2*x[8])+x[3]/x[5]))
  return(val)
}
set.seed(1)
p=8
library(sensitivity)
m=10000
A=matrix(runif(m*p),nrow=m)
B=matrix(runif(m*p),nrow=m)
borehole=function(X) apply(X,1,f)
a.sen=soboljansen(model = borehole, X1 = data.frame(A), X2 = data.frame(B))
print(a.sen)
sf=pmax(unlist(a.sen$S),0)
st=pmax(unlist(a.sen$T),0)
x=matrix(c(sf,st-sf),nrow=2,byrow=TRUE)
colnames(x)=paste("x",sep="",1:p)
rownames(x)=c("first order","total-first order")
#barplot(x,col=c("blue","pink"),legend.text=TRUE,density = 20,cex.names = 2,angle = c(45,135),args.legend = list("cex"=3,"bty"="n"))
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
barplot(x,col=c("lightblue","pink"),legend.text=TRUE,cex.names = 2,args.legend = list("cex"=3,"bty"="n"),ylab="total Sobol' index",cex.lab=2)

#Figure 6.1 (right)

library(SFDesign)
p=8
n=10*p
D=maxpro.optim(maxproLHD(n,p)$design)$design
y=apply(D,1,f)

a.fit=Fit.Kriging(D,y,kernel.parameters = list(type="Gaussian"))
fhat=function(x,fit) Predict.Kriging(fit,matrix(x,ncol=p))$mean



N=10*p
D.eval=uniformLHD(N,p)$design
f.eval=apply(D.eval,1,fhat, a.fit)
f0=mean(f.eval)
x.ev=seq(0,1,length=20)
val=x.ev
main.effect=function(ind)
{
  D0=D.eval
  for(i in 1:20)
  {
    D0[,ind]=rep(x.ev[i],N)
    val[i]=mean(apply(D0,1,fhat, a.fit))
  }
  return(val-f0)
}
M=matrix(0,nrow=20,ncol=p)
for(j in 1:p) M[,j]=main.effect(j)
matplot(x.ev,M,type="l",lty=1:p,col=1:p,lwd=3, xlab="x",ylab="main effects", cex.lab=2)
legend("topleft",c("x1","x2","x3","x4","x5","x6","x7","x8"),bty="n",lty=1:p,col=1:p,lwd=3, cex=1.5)
par(mfrow=c(1,1))

####################################################
#Section 6.1.2 Derivative-based Sensitivity Measures
####################################################

#Figure 6.2

library(SFDesign)
set.seed(1)
f=function(x) {
  lower<-c(0.05,100,63070,990,63.1,700,1120,9855)
  upper<-c(0.15,50000,115600,1110,116,820,1680,12045)
  x=lower+x*(upper-lower)
  val=2*pi*x[3]*(x[4]-x[6])/(log(x[2]/x[1])*(1+2*x[7]*x[3]/(log(x[2]/x[1])*x[1]^2*x[8])+x[3]/x[5]))
  return(val)
}
borehole=function(X) apply(X,1,f)
p=8
m=4
nu=matrix(0,nrow=100,ncol=p)
for(i in 1:100){
  X=matrix(runif(m*p),ncol=p)
  ad=delsa(model=borehole, X0=X,varprior = rep(1,p))
  nu[i,]=colMeans(ad$deriv^2)
  nu[i,]=nu[i,]/sum(nu[i,])
}
lab=paste("x",sep="",1:p)
colnames(nu)=lab
boxplot(nu,main="Derivative-based sensitivity measure",cex.main=2,cex.axis=1.5,col="lightblue",ylab = expression(bar(nu)),cex.lab=2)


#Figure 6.3

p=4
f=function(x) {
  val=0
  for(k in 1:p) val=val+sin((k^2+1)*pi*x[k])/k
  return(val)
}
f1=function(X) apply(X,1,f)
X=uniformLHD(10*p,p)$design
a=delsa(model=f1, X0=X,varprior = rep(1,p))
nu=colMeans(a$deriv^2)
nu=nu/sum(nu)

m=10000
A=matrix(runif(m*p),nrow=m)
B=matrix(runif(m*p),nrow=m)
a.sen=soboljansen(model = f1, X1 = data.frame(A), X2 = data.frame(B))
tot=unlist(a.sen$T)

par(mfrow=c(1,2))
curve(sin(2*pi*x),from=0,to=1,lwd=3,ylab="y",cex.lab=2)
for(k in 2:p) curve(sin((k^2+1)*pi*x)/k,from=0,to=1,add=TRUE,col=k,lty=1,lwd=3)
legend("topright",legend = paste("x",sep="",1:p),col=1:p,lty=1,lwd=3,bty="n",cex=1.5)

plot(1:p,tot,col="blue",axes="F",cex=3,pch=16,xlab="",ylab="sensitivity measure", main="Variance-based vs Derivative-based",cex.main=2,cex.lab=2)
lab=paste("x",sep="",1:p)
axis(1,at=1:p,labels=lab,cex.axis=2)
axis(2,seq(0,.8,.1),cex.axis=2)
points(1:p,nu, col="red",cex=3,pch=4,lwd=3)
legend("topright",legend=c("T",expression(nu)),col=c("blue","red"),pch=c(16,4),cex=1.5,lwd=3,lty=NA)
par(mfrow=c(1,1))

####################################
#Section 6.2 Morris screening design
####################################

#Figure 6.4

set.seed(10)
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
x=((1:4)-1)/3
A=cbind(x,sample(x))
delta=1/3
C1=C2=A
for(i in 1:4){
  C1[i,1]=A[i,1]+delta
  if(C1[i,1]<0||C1[i,1]>1) C1[i,1]=A[i,1]-delta
  C2[i,1]=C1[i,1]
  C2[i,2]=C1[i,2]+delta
  if(C2[i,2]<0||C2[i,2]>1) C2[i,2]=C1[i,2]-delta
}
D=rbind(A,C1,C2)
dim(unique(D))[1]

plot(A,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),lwd=3,cex=2,main="Morris Design with Strict OFAT",cex.main=2,cex.lab=2)
points(C1,col=2,pch=2,lwd=3,cex=2)
points(C2,col=3,pch=3,lwd=3,cex=2)
for(i in 1:4){
  arrows(A[i,1],A[i,2],C1[i,1],C1[i,2],col=4,lty=1,lwd=2)
  arrows(C1[i,1],C1[i,2],C2[i,1],C2[i,2],col=4,lty=1,lwd=2)
} 


C1=C2=A
for(i in 1:4){
  C1[i,1]=A[i,1]+delta
  if(C1[i,1]<0||C1[i,1]>1) C1[i,1]=A[i,1]-delta
  C2[i,2]=C1[i,2]+delta
  if(C2[i,2]<0||C2[i,2]>1) C2[i,2]=C1[i,2]-delta
}
D=rbind(A,C1,C2)
dim(unique(D))[1]

plot(A,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),lwd=3,cex=2,main="Morris Design with Standard OFAT",cex.main=2,cex.lab=2)
points(C1,col=2,pch=2,lwd=3,cex=2)
points(C2,col=3,pch=3,lwd=3,cex=2)
for(i in 1:4){
  arrows(A[i,1],A[i,2],C1[i,1],C1[i,2],col=4,lty=1,lwd=2)
  arrows(A[i,1],A[i,2],C2[i,1],C2[i,2],col=4,lty=1,lwd=2)
} 
par(mfrow=c(1,2))

#Figure 6.5

library(sensitivity)
set.seed(1)
f=function(x) {
  lower<-c(0.05,100,63070,990,63.1,700,1120,9855)
  upper<-c(0.15,50000,115600,1110,116,820,1680,12045)
  x=lower+x*(upper-lower)
  val=2*pi*x[3]*(x[4]-x[6])/(log(x[2]/x[1])*(1+2*x[7]*x[3]/(log(x[2]/x[1])*x[1]^2*x[8])+x[3]/x[5]))
  return(val)
}
p=8
borehole=function(X) apply(X,1,f)
a <- morris(model = borehole, factors = 8, r = 4,
            design = list(type = "oat", levels = 4, grid.jump = 2))
print(a)

mustar=apply(a$ee, 2, function(x) mean(abs(x)))
sigma <- apply(a$ee, 2, sd)
mu=abs(apply(a$ee, 2, function(x) mean(x)))

par(mfrow=c(1,2))
sel=(mu>10)
plot(mu,sigma,type="n",xlab = expression(abs(mu)),ylab=expression(sigma),cex.lab=1.5,main="Borehole function",cex.main=2,cex=2)
lab=paste("x",sep="",1:p)
text(mu[sel],sigma[sel],lab[sel],cex=2,col=2)
points(mu[!sel],sigma[!sel])
mustar=matrix(0,nrow=100,ncol=p)
for(i in 1:100){
am <- morris(model = borehole, factors = p, r = 4, design = list(type = "oat", levels = 4, grid.jump = 2))
mustar[i,]=apply(am$ee, 2, function(x) mean(abs(x)))
}
colnames(mustar)=lab
boxplot(mustar, main="Morris screening design",cex.main=2,cex.axis=1.5,col="lightblue",ylab = expression(mu^"*"),cex.lab=1.5)
par(mfrow=c(1,1))

##########################
#Section 6.3 MOFAT designs
##########################


#Figure 6.6 (Sobol designs)
par(mfrow=c(1,2))
set.seed(3)
par(mar=c(5,5,3,2))
A=matrix(runif(4*2),nrow=4,ncol=2)
B=matrix(runif(4*2),nrow=4,ncol=2)
C1=cbind(B[,1],A[,2])
C2=cbind(A[,1],B[,2])
D=rbind(A,C1,C2)

plot(A,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),lwd=3,cex=2,main="Sobol' Design",cex.main=2, cex.lab=2)
points(C1,col=2,pch=2,lwd=3,cex=2)
points(C2,col=3,pch=3,lwd=3,cex=2)
for(i in 1:4){
  arrows(D[i,1],D[i,2],D[i+4,1],D[i+4,2],col=4,lwd=2)
  arrows(D[i,1],D[i,2],D[i+8,1],D[i+8,2],col=4,lwd=2)
} 
set.seed(1)
par(mar=c(5,5,3,2))
A=matrix(runif(4*2),nrow=4,ncol=2)
B=matrix(runif(4*2),nrow=4,ncol=2)
C1=cbind(B[,1],A[,2])
C2=cbind(A[,1],B[,2])
D=rbind(A,C1,C2)

plot(A,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),lwd=3,cex=2,main="Sobol' Design",cex.main=2, cex.lab=2)
points(C1,col=2,pch=2,lwd=3,cex=2)
points(C2,col=3,pch=3,lwd=3,cex=2)
for(i in 1:4){
  arrows(D[i,1],D[i,2],D[i+4,1],D[i+4,2],col=4,lwd=2)
  arrows(D[i,1],D[i,2],D[i+8,1],D[i+8,2],col=4,lwd=2)
} 

par(mfrow=c(1,1))

#Figure 6.7 (MOFAT designs)
library(MOFAT)
p=2
l=4
s=1
set.seed(s)
D=mofat(p,l,method="uniform")
A=D[1:l,]
C1=D[(l+1):(2*l),]
C2=D[(2*l+1):(3*l),]

par(mfrow=c(1,2))
plot(A,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),lwd=3,cex=2,main="MOFAT Design",cex.main=2,cex.main=2, cex.lab=1.5)
points(C1,col=2,pch=2,lwd=4,cex=2)
points(C2,col=3,pch=3,lwd=4,cex=2)
for(i in 1:l){
  arrows(A[i,1],A[i,2],C1[i,1],C1[i,2],col=4,lty=1,lwd=2)
  arrows(A[i,1],A[i,2],C2[i,1],C2[i,2],col=4,lty=1,lwd=2)
}
min(dist(D))
dim(unique(D))[1]

set.seed(s)
D=mofat(p,l,method="projection")
A=D[1:l,]
C1=D[(l+1):(2*l),]
C2=D[(2*l+1):(3*l),]
plot(A,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),lwd=3,cex=2,main="MOFAT Design with projections",cex.main=2, cex.lab=1.5)
points(C1,col=2,pch=2,lwd=4,cex=2)
points(C2,col=3,pch=3,lwd=4,cex=2)
for(i in 1:l){
  arrows(A[i,1],A[i,2],C1[i,1],C1[i,2],col=4,lty=1,lwd=2)
  arrows(A[i,1],A[i,2],C2[i,1],C2[i,2],col=4,lty=1,lwd=2)
}
min(dist(D))
dim(unique(D))[1]
par(mfrow=c(1,1))

#Figure 6.8

#borehole function
library(MOFAT)
f=function(x) {
  lower<-c(0.05,100,63070,990,63.1,700,1120,9855)
  upper<-c(0.15,50000,115600,1110,116,820,1680,12045)
  x=lower+x*(upper-lower)
  val=2*pi*x[3]*(x[4]-x[6])/(log(x[2]/x[1])*(1+2*x[7]*x[3]/(log(x[2]/x[1])*x[1]^2*x[8])+x[3]/x[5]))
  return(val)
}

p=8
m=4
mustar=matrix(0,nrow=100,ncol=p)
for(i in 1:100){
  d=mofat(p=8,l=m)
  y=apply(d,1,f)
  mustar[i,]=measure(d,y)$mustar
  mustar[i,]=mustar[i,]/sum(mustar[i,])
}
lab=paste("x",sep="",1:p)
colnames(mustar)=lab
par(mar=c(5,5,3,2))
boxplot(mustar, main="MOFAT design",cex.main=2,cex.axis=1.5,col="lightblue",ylab = expression(mu^"*"),cex.lab=1.75)


