################
#Chapter 8: Fractional Factorial Designs
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2026). Experimental Design for Data Science and Engineering. CRC Press.
######################################

#######################################
#Section 8.1: Two-level Designs
######################################

#Figure 8.1

rm(list=ls())
A=rep(c(-1,1),c(4,4))
B=rep(rep(c(-1,1),c(2,2)),2)
C=rep(c(-1,1),4)
D=A*B*C
X=model.matrix(rep(1,8)~(A+B+C)^3)
E=model.matrix(rep(1,8)~(A+B+C+D)^4)
solve(t(X)%*%X)%*%t(X)%*%E
y=c(504, 984, 928, 808, 992, 784, 464, 976)
a=lm(y~(A+B+C)^3)
summary(a)

eff=a$coef[-1]
names(eff)=c("1","2","3","12","13","23","123")
I=length(eff)
u=qnorm(.5+.5*(1:I-.5)/I)
saeff=sort(abs(eff))
par(mar=c(5,5,4,2))
plot(u,saeff,type="n",xlab="half-normal quantiles",ylab="absolute coefficients", main="Half-normal plot",cex.main=3,cex.axis=2,cex.lab=2, xlim=c(0,max(u)+.1))
text(u,saeff,names(saeff),cex=2, col=c(rep(4,4),rep(2,3)))

#######################################
# Section 8.2 Bayesian-inspired Designs
#######################################

#Figure 8.2

library(rkriging)
d=cbind(A,B,C,D)
a=Fit.Kriging(d,y, kernel.parameters = list(type="Gaussian"))
theta=Get.Kriging.Parameters(a)$lengthscale
mu=Get.Kriging.Parameters(a)$mu
gam=(1-exp(-2^2/(2*theta)))/(1+exp(-2^2/(2*theta)))

X=model.matrix(lm(y~(A+B+C+D)^4),data=data.frame(A,B,C,D))
G=diag(c(model.matrix(lm(1~(gam[1]+gam[2]+gam[3]+gam[4])^4))))
alpha=G%*%t(X)%*%solve(X%*%G%*%t(X))%*%(y-mu)
eff=alpha[-1]
names(eff)=c("1","2","3","4","12","13","14","23","24","34","123","124","134","234","1234")
I=length(eff)
u=qnorm(.5+.5*(1:I-.5)/I)
saeff=sort(abs(eff))
par(mar=c(5,5,4,2))
plot(u,saeff,type="n",xlab="half-normal quantiles",ylab="absolute coefficients", main="Half-normal plot",cex.main=3,cex.axis=2,cex.lab=2, xlim=c(0,max(u)+.1))
size=c(rep(1.5,12),rep(2,3))
text(u,saeff,names(saeff),cex=size,col=c(rep(4,12),rep(2,3)))
gam[4]/(gam[4]+gam[1]*gam[2]*gam[3])*165

library(HiGarrote)
HiGarrote(d,y)

#Analysis of cast fatigue experiment
setwd("~/R commands/Book")
cast=read.table(file="cast.txt",h=TRUE)
summary(lm(y~.,data=cast))
summary(lm(y~F1,data=cast))

HiGarrote(cast[,1:7],y=cast[,8])
summary(lm(y~F1*G+D+D:G,data=cast))

###############################################
# Section 8.3 Designs with more than two levels
###############################################

#Figure 8.3 (this will take time!)

#L9 with two replicates
library(DoE.base)
library(SFDesign)
library(twinning)

full=full.factorial(4,3)
D=undesign(oa.design(nfactors=4, nlevels=3, randomize = FALSE))
D=data.matrix(D)
D1=rbind(D,D)
GWLP(D1)
D1=(D1-1)/2
mp1=maxpro.crit(D1,delta=1/3^2)/3^2
en1=energy(full,D1)

#L18
D=undesign(oa.design(ID=L18,randomize=FALSE))
D=data.matrix(D)
D2=D[,2:5]
GWLP(D2)
D2=(D2-1)/2
mp2=maxpro.crit(D2,delta=1/3^2)/3^2
en2=energy(full,D2)

D=undesign(oa.design(ID=L18,randomize=FALSE))
D=data.matrix(D)
D3=D[,3:6]
GWLP(D3)
mp3=maxpro.crit((D3-1)/2,delta=1/3^2)/3^2
en3=energy(full,(D3-1)/2)

A=array(data=rep(0,18*4*27),dim=c(18,4,27))
D=undesign(oa.design(ID=L18,randomize=FALSE))
D=data.matrix(D)
D=D[,3:6]
D1 <- (D[,1]-1)/2
A1 <- lapply(1:3, function(i) {
  D2 <- D[,2]
  if(i == 2) {
    D2 <- ifelse(D2 == 1, 1, ifelse(D2 == 2, 3, 2))
  } else if(i == 3) {
    D2 <- ifelse(D2 == 1, 2, ifelse(D2 == 2, 1, 3))
  }
  D2=(D2-1)/2
  return(D2)
})

A2 <- lapply(1:3, function(i) {
  D3 <- D[,3]
  if(i == 2) {
    D3 <- ifelse(D3 == 1, 1, ifelse(D3 == 2, 3, 2))
  } else if(i == 3) {
    D3 <- ifelse(D3 == 1, 2, ifelse(D3 == 2, 1, 3))
  }
  D3=(D3-1)/2
  return(D3)
})

A3 <- lapply(1:3, function(i) {
  D4 <- D[,4]
  if(i == 2) {
    D4 <- ifelse(D4 == 1, 1, ifelse(D4 == 2, 3, 2))
  } else if(i == 3) {
    D4 <- ifelse(D4 == 1, 2, ifelse(D4 == 2, 1, 3))
  }
  D4=(D4-1)/2
  return(D4)
})

a <- expand.grid(c(1,2,3), c(1,2,3), c(1,2,3))
for(i in 1:27) {
  idx <- a[i,]
  A[,,i] <- cbind(D1, A1[[idx$Var1]], A2[[idx$Var2]], A3[[idx$Var3]])
}

val.mp=val.en=numeric(27)
for(k in 1:27) {
  val.mp[k]=maxpro.crit(A[,,k],delta=1/3^2)/3^2
  val.en[k]=energy(full,A[,,k])
}
val.mp
val.en

full=full.factorial(4,3)
N=1000000
val.mp.sim=val.en.sim=numeric(N)
for(i in 1:N){
  set.seed(i)
  D=full[sample(1:81,18),]
  val.mp.sim[i]= maxpro.crit(D,delta=1/3^2)/3^2
  val.en.sim[i]=energy(full,D)
}
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
plot(density(val.mp.sim),xlim=c(.4,.52),lwd=2,main = "MaxPro crietrion for 18-run designs",xlab = "MaxPro measure",cex.lab=2,cex.axis=2, cex.main=1.75)
abline(v=unique(val.mp),lwd=2,col=4)
abline(v=unique(mp1),lwd=3,lty=3,col=2)
abline(v=unique(mp2),lwd=3,lty=2,col=3)
text(mp1,3,expression(D[1]),cex=2)
text(mp2,3,expression(D[2]),cex=2)
text(mp3,3,expression(D[3]),cex=2)

plot(density(val.en.sim),xlim=c(2.7,3.05),lwd=2,main = "Energy distance for 18-run designs",xlab = "energy distance",cex.lab=2,cex.axis=2, cex.main=1.75)
abline(v=unique(val.en),lwd=2,col=4)
abline(v=unique(en1),lwd=3,lty=3,col=2)
abline(v=unique(en2),lwd=3,lty=2,col=3)
text(en1+.015,3,expression(D[1]),cex=2)
text(en2+.015,3,expression(D[2]),cex=2)
text(en3-.015,3,expression(D[3]),cex=2)
par(mfrow=c(1,1))


#CCD
x1=rep(c(-1,1),c(4,4))
x2=rep(rep(c(-1,1),c(2,2)),2)
x3=rep(c(-1,1),4)
x4=x1*x2
D=cbind(x1,x2,x3,x4)
axial=matrix(rep(0,8*4),nrow=8)
for(i in 1:4) {
  axial[2*i-1,i]=-1
  axial[2*i,i]=1
}
D=rbind(D,axial,rep(0,4))
D=(D+1)/2
maxpro.crit(D,delta=1/3^2)/3^2
energy(full,D)

############################
#Section 8.4 Mixture designs
############################

#Figure 8.4

rm(list=ls())
library(mixexp)
D1=SLD(3,2)
D2=SCD(3)
#DesignPoints(D1)
#DesignPoints(D2)

library("Ternary")
par(mar=c(.1,.1,.1,.1))
par(mfrow=c(1,2))
TernaryPlot(axis.labels = seq(0, 1, by = 0.1),
            atip = expression(x[1]), btip = expression(x[3]), ctip = expression(x[2]),
            alab = "Fraction x1", blab = "Fraction x3", clab = "Fraction x2",
            lab.cex = 1, tip.cex = 2)
TernaryPoints(D1,
              cex = 3,  # Point size
              col = 4,  # Point colour
              pch = 16  # Plotting symbol (16 = filled circle)
)

TernaryPlot(axis.labels = seq(0, 1, by = 0.1),
            atip = expression(x[1]), btip = expression(x[3]), ctip = expression(x[2]),
            alab = "Fraction x1", blab = "Fraction x3", clab = "Fraction x2",
            lab.cex = 1, tip.cex = 2)
TernaryPoints(D2,
              cex = 3,  # Point size
              col = 4,  # Point colour
              pch = 16  # Plotting symbol (16 = filled circle)
)
par(mfrow=c(1,1))

#Figure 8.5

library(spacefillr)
library(support)
p=3
N=10000
u=generate_sobol_set(N,p-1)
x=matrix(0,nrow=N,ncol=p)
x[,1]=(1-u[,1]^(1/2))
x[,2]=u[,1]^(1/2)*(1-u[,2])
x[,3]=u[,1]^(1/2)*u[,2]

set.seed(1)
n=11
D1=sp(n,p,dist.samp = x)$sp

par(mar=c(.1,.1,.1,.1))
par(mfrow=c(1,2))
TernaryPlot(axis.labels = seq(0, 1, by = 0.1),
            atip = expression(x[1]), btip = expression(x[3]), ctip = expression(x[2]),
            alab = "Fraction x1", blab = "Fraction x3", clab = "Fraction x2",
            lab.cex = 1, tip.cex = 2)
TernaryPoints(D1,
              cex = 3,  # Point size
              col = 4,  # Point colour
              pch = 16  # Plotting symbol (16 = filled circle)
)

CAND=x
C1=(x[,1]+x[,2]<.7)
C2=(x[,1]+x[,2]>.3)
CAND=CAND[C1&C2,]

D2=sp(n,p,dist.samp = CAND)$sp

par(mar=c(.1,.1,.1,.1))
TernaryPlot(axis.labels = seq(0, 1, by = 0.1),
            atip = expression(x[1]), btip = expression(x[3]), ctip = expression(x[2]),
            alab = "Fraction x1", blab = "Fraction x3", clab = "Fraction x2",
            lab.cex = 1, tip.cex = 2)
TernaryPoints(D2,
              cex = 3,  # Point size
              col = 4,  # Point colour
              pch = 16  # Plotting symbol (16 = filled circle)
)
par(mfrow=c(1,1))
