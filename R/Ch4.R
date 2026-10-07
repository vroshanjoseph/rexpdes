################
#Chapter 4: Space-filling Designs
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2026). Experimental Design for Data Science and Engineering. CRC Press.
######################################

#######################################
#Section 4.1: Clustering-based Designs
######################################
rm(list=ls())
library(FNN)

library(LHD)

#Figure 4.1

library(SFDesign)
library(mined)
library(deldir)
p=2
n=7
N=n*p*10000
X=Lattice(10000,2)

a2=kmeans(X,centers = n,iter.max = 1000,nstart = 100)
D2=a2$centers
cluster.error(D2)
D=clustering.design(n,p)$design
cluster.error(D)
D=clustering.design(n,p, D.ini = D2)$design
cluster.error(D)

par(mar=c(5,5,4,2))
plot(D,xlim=c(0,1),ylim=c(0,1), main="clustering-based design", xlab=expression(x[1]),ylab=expression(x[2]), cex.main=3,cex.lab=2,cex.axis=2,pch=16,col="blue",cex=3)
vor=deldir(D[,1],D[,2])
plot.deldir(vor,add=TRUE,wlines = "tess",xlim=c(0,1),ylim=c(0,1),showpoints = FALSE, lwd=3)

############################
#Section 4.2 Minimax designs
############################

#Figure 4.2

library(MASS)
library(plotrix)
p=2
n=7
DmM=matrix(c(.5,.5,.5,0,.5,1,1/3-1/12*sqrt(7),.25,1/3-1/12*sqrt(7),.75,2/3+1/12*sqrt(7),.25,2/3+1/12*sqrt(7),.75),nrow=7,ncol=2,byrow=T)
eqscplot(DmM[,1],DmM[,2],xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",cex=3,cex.main=3,cex.lab=2,cex.axis=2, main="minimax design")
for(i in 1:7) draw.circle(DmM[i,1],DmM[i,2],radius=1/6*(sqrt(7)-1),lwd=2)
mMdist(DmM)$dist

#Generate minimax designs

library(minimaxdesign)
set.seed(2)

for(i in 1:10){
  set.seed(i)
D=minimax(N=n,p=2,q=10)
print(c(i,mMdist(D)$dist))
}

set.seed(10) #best seed
D=minimax(N=n,p=2,q=10)
points(D,pch=2,col=3)
mMdist(D)$dist

library(LatticeDesign)
a=InterleavedMinimaxD(p,n)
mMdist(a$Design)$dist

#############################
#Section 4.3 Maximin designs
#############################



#Figure 4.3
p=2
n=7
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
D1=matrix(c(0,.1,.1,0,0,1,1,1,1,.5,1,0,.5,.5),nrow=n,ncol=p,byrow = TRUE)
plot(D1,xlim=c(0,1),ylim=c(0,1), xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=2,cex.axis=2,pch=16,col="blue",cex=3, main=paste0("index=",1),cex.main=3)
segments(0,.1,.1,0,col=2,lwd=3)
D2=matrix(c(0,.1,.1,0,0,1,.9,1,1,.9,1,0,.5,.5),nrow=n,ncol=p,byrow = TRUE)
plot(D2,xlim=c(0,1),ylim=c(0,1), xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=2,cex.axis=2,pch=16,col="blue",cex=3, main=paste0("index=",2),cex.main=3)
segments(0,.1,.1,0,col=2,lwd=3)
segments(.9,1,1,.9,col=2,lwd=3)
par(mfrow=c(1,1))

#Figure 4.4

library(SFDesign)
set.seed(1)
D=maximinLHD(n,p)$design
D = maximin.optim(D, sa=TRUE, find.best.ini = TRUE)$design
r=min(dist(D))/2
r
eqscplot(D[,1],D[,2],xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",cex=3,cex.main=3,cex.lab=2,cex.axis=2, main="maximin design")
for(i in 1:n) draw.circle(D[i,1],D[i,2],radius=r,lwd=2)


#Figure 4.5

rm(list=ls())
p=2
n=7

set.seed(1)
library(SFDesign)
D0=matrix(c(.5,.5),nrow=1,ncol=p)
CAND=full.factorial(p=2,level=5)
D=maximin.augment(n=7,p=2,D.ini=D0, candidate = CAND)
min(dist(D))/2
D2=maximin.augment(n=13,p=2,D.ini=D,candidate = CAND)
min(dist(D2))/2
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
plot(D,xlim=c(0,1),ylim=c(0,1),cex=3,lwd=3, xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=2,cex.axis=2,pch=16,col="blue")
plot(D,xlim=c(0,1),ylim=c(0,1),cex=3,lwd=3, xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=2,cex.axis=2,pch=16,col="blue")
points(D2[8:13,], cex=3,lwd=5, xlab=expression(x[1]),ylab=expression(x[2]), cex.lab=2,cex.axis=2,pch=4,col="red")
par(mfrow=c(1,1))

#optimal 13-point maximin design
D1=maximin.optim(D2, sa=TRUE, find.best.ini = TRUE)$design
min(dist(D1))/2

####################################
#Section 4.4 Latin hypercube designs
####################################

#Figure 4.6

par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
D=matrix(c(0,0,0,1,1,0,1,1),nrow=4,ncol=2,byrow=TRUE)
plot(D,pch=16,xaxt="n",yaxt="n",cex=4,bty="n",col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2, main="maximin design",cex.main=3)
axis(1,at=seq(0,1,by=.25),cex.axis=2)
axis(2,at=seq(0,1,by=.25),cex.axis=2)
points(cbind(D[,1],-.15),col=2,pch=16,cex=4)
arrows(0,0,0,-.125,lwd=2,lty=2)
arrows(1,0,1,-.125,lwd=2,lty=2)
arrows(0,1,0,.875,lwd=2,lty=2)
arrows(1,1,1,.875,lwd=2,lty=2)
points(cbind(-.15,D[,2]),col=2,pch=16,cex=4)
arrows(0,0,-.125,0,lwd=2,lty=2)
arrows(0,1,-.125,1,lwd=2,lty=2)
arrows(1,0,.875,0,lwd=2,lty=2)
arrows(1,1,.875,1,lwd=2,lty=2)

D=matrix(c(0.25,0.25,.25,.75,.75,.25,.75,.75),nrow=4,ncol=2,byrow=TRUE)
plot(D,pch=16,xaxt="n",yaxt="n",cex=4,bty="n",col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2, main="minimax design",cex.main=3)
axis(1,at=seq(0,1,by=.25),cex.axis=2)
axis(2,at=seq(0,1,by=.25),cex.axis=2)
points(cbind(D[,1],-.15),col=2,pch=16,cex=4)
arrows(.25,.25,.25,.125,lwd=2,lty=2)
arrows(.75,.25,.75,.125,lwd=2,lty=2)
arrows(.25,.75,.25,.625,lwd=2,lty=2)
arrows(.75,.75,.75,.625,lwd=2,lty=2)
points(cbind(-.15,D[,2]),col=2,pch=16,cex=4)
arrows(.25,.25,.125,.25,lwd=2,lty=2)
arrows(.25,.75,.125,.75,lwd=2,lty=2)
arrows(.75,.25,.625,.25,lwd=2,lty=2)
arrows(.75,.75,.625,.75,lwd=2,lty=2)
par(mfrow=c(1,1))


#Figure 4.7

p=2
n=4
library(SFDesign)
D=maximinLHD(n,p)$design

par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
plot(D,pch=16,xaxt="n",yaxt="n",cex=4,bty="n",col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2)
axis(1,at=seq(0,1,by=.25),cex.axis=2)
axis(2,at=seq(0,1,by=.25),cex.axis=2)
for(i in 1:(n+1))
{
  segments((i-1)/n,0,(i-1)/n,1)
  segments(0,(i-1)/n,1,(i-1)/n)
}
points(cbind(D[,1],-.15),col=2,pch=16,cex=4)
points(cbind(-.15,D[,2]),col=2,pch=16,cex=4)

D=(cbind(1:n,1:n)-.5)/n
plot(D,pch=16,xaxt="n",yaxt="n",cex=4,bty="n",col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2)
axis(1,at=seq(0,1,by=.25),cex.axis=2)
axis(2,at=seq(0,1,by=.25),cex.axis=2)
for(i in 1:(n+1))
{
  segments((i-1)/n,0,(i-1)/n,1)
  segments(0,(i-1)/n,1,(i-1)/n)
}
points(cbind(D[,1],-.15),col=2,pch=16,cex=4)
points(cbind(-.15,D[,2]),col=2,pch=16,cex=4)
par(mfrow=c(1,1))


#Figure 4.8

p=2
n=7
D=maximinLHD(n,p)$design
plot(D,pch=16,cex=4,bty="n",col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2,cex.axis=2, main="MmLHD", cex.main=3)
for(i in 1:(n+1))
{
  segments((i-1)/n,0,(i-1)/n,1)
  segments(0,(i-1)/n,1,(i-1)/n)
}

#Figure 4.9

p=2
n=20
set.seed(1)
D=maximinLHD(n,p)$design
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
plot(D,pch=16,cex=4,bty="n",col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2,cex.axis=2, main="MmLHD", cex.main=3)
for(i in 1:(n+1))
{
  segments((i-1)/n,0,(i-1)/n,1)
  segments(0,(i-1)/n,1,(i-1)/n)
}
D1=(1+cos(pi*D))/2
plot(D1,pch=16,cex=4,bty="n",col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2,cex.axis=2, main="Chebyshev-MmLHD", cex.main=3)
d1=sort(D1[,1])
l1=c(0,(d1[-1]+d1[-n])/2,1)
d2=sort(D1[,2])
l2=c(0,(d2[-1]+d2[-n])/2,1)
for(i in 2:(n+2))
{
  segments(l1[i-1],0,l1[i-1],1)
  segments(0,l2[i-1],1,l2[i-1])
}
par(mfrow=c(1,1))

###########################
#Section 4.5 MaxPro designs
###########################

#Figure 4.10

library(MaxPro)
p=2
n=7
a=MaxPro(MaxProLHD(n,p)$Design)
D1=a$Design

library(SFDesign)
p=2
n=20
set.seed(10)
D = maxpro.optim(maxproLHD(n, p)$design, iteration = 100)$design
maxpro.crit(D)

par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
plot(D1,pch=16,cex=4,bty="n",col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2,cex.axis=2, main="MaxPro(7,2)", cex.main=3)
n=7
d1=sort(D1[,1])
l1=c(0,(d1[-1]+d1[-n])/2,1)
d2=sort(D1[,2])
l2=c(0,(d2[-1]+d2[-n])/2,1)
for(i in 2:(n+2))
{
  segments(l1[i-1],0,l1[i-1],1)
  segments(0,l2[i-1],1,l2[i-1])
}

n=20
plot(D,pch=16,cex=4,bty="n",col="blue",xlim=c(-.1,1.1),ylim=c(-.1,1.1),xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2,cex.axis=2, main="MaxPro(20,2)", cex.main=3)
d1=sort(D[,1])
l1=c(0,(d1[-1]+d1[-n])/2,1)
d2=sort(D[,2])
l2=c(0,(d2[-1]+d2[-n])/2,1)
for(i in 2:(n+2))
{
  segments(l1[i-1],0,l1[i-1],1)
  segments(0,l2[i-1],1,l2[i-1])
}
par(mfrow=c(1,1))


#Figure 4.11

p=5
n=50
set.seed(2)
library(SFDesign)
D=maxpro.optim(maxproLHD(n,p)$design)$design
maxpro.crit(D)

#pairs(D)

library(car)
par(mar=c(5,5,4,2))
scatterplotMatrix(D,diagonal=list(method ="histogram"),smooth=FALSE, regLine = FALSE,pch=16)

#Figure 4.12

D0=maximinLHD(n,p)$design
library(gtools)
min_mindist=function(M)
{
  n=dim(M)[1]
  p=dim(M)[2]
  a=rep(p,(p-1))
  for(j in 1:(p-1))
  {
    comb=combinations(p,j)
    k=dim(comb)[1]
    temp=mat.or.vec(k,1)
    for(i in 1:k){
      temp[i]=min(dist(M[,comb[i,]]))
    }
    a[j]=min(temp)
  }
  minfull=min(dist(M))
  return(c(a,minfull))
}
d1=min_mindist(D)
d2=min_mindist(D0)
par(mar=c(5,5,4,2))
matplot(1:p,cbind(d1,d2),type="pp",pch=c(1,4),col=c(4,2),ylim=c(0,.55),bty="n",axes=F,xlab="projection dimension",ylab="minimum of minimum distances",cex=2, cex.lab=2,lwd=2)
axis(1,1:p,cex.axis=2)
axis(2,seq(0,.52,.1),cex.axis=2)
lines(d1,lwd=3,col=4)
lines(d2,lwd=3,lty=2,col=2)
legend("topleft",c("MaxPro","MmLHD"),pch=c(1,4),lty=c(1,2),bty='n',lwd=3,col=c(4,2),cex=2)

#Figure 4.13

rm(list=ls())
library(SFDesign)
library(fitdistrplus)
p=6
n=64
set.seed(1)
D0=maximin.optim(maximinLHD(n,p)$design, find.best.ini = TRUE)$design
D1=maxpro.optim(maxproLHD(n,p)$design)$design
d0=c(dist(D0))
d1=c(dist(D1))
a=fitdist(d1/sqrt(p),"beta")
a$estimate
par(mar=c(5,5,4,2))
hist(d0,prob=TRUE,xlim=c(0,sqrt(p)), col=rgb(1,0,0,1/4),xlab="pairwise distance", main="distance distribution",cex.lab=2,cex.axis=2, cex.main=3)
hist(d1,prob=TRUE, col=rgb(0,0,1,1/4),breaks=30,add=TRUE)
curve(dbeta(x/sqrt(p),a$estimate[1],a$estimate[2])/sqrt(p),from=0,to=sqrt(p),add=TRUE,col=3,lwd=3)
mycol=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4),3)
legend("topleft", legend = c("Maximin","MaxPro", "Beta(10.1,12.7)"),col=mycol, lty = c(1,1,1),bty="n",lwd=c(15,15,3),cex=2)

#Figure 4.14

branin=function(x){
    x1 <- x[1] * 15 - 5
    x2 <- x[2] * 15
    val1=(x2 - 5.1/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2
    val2=10 * (1 - 1/(8 * pi)) * cos(x1) + 10
    return(val1+val2)
}
library(MaxPro)
p=2
n=10
set.seed(123)
D=MaxPro(MaxProLHD(n,p)$Design)$Design
CAND=CandPoints(N=10000,2)
A=MaxProAugment(ExistDesign = D,CandDesign = CAND,nNew = 10)$Design[(n+1):20,]
y=apply(D,1,branin)
library(rkriging)
a<-Fit.Kriging(D,y,kernel.parameters = list(type="Gaussian"))
true=apply(A,1,branin)
pred=Predict.Kriging(a,A)$mean
rmse0=sqrt(mean((pred-true)^2))

rmse=numeric(100)
for(i in 1:100){
  u=matrix(runif(p*10),nrow=10,ncol=p)
  pred=Predict.Kriging(a,u)$mean
  true=apply(u,1,branin)
  rmse[i]=sqrt(mean((pred-true)^2))
}

N.plot=250
p1=seq(0,1,length=N.plot)
p2=seq(0,1,length=N.plot)
fc=matrix(0,N.plot,N.plot)
for(i in 1:N.plot)
{
  for(j in 1:N.plot)
    fc[i,j]=branin(c(p1[i],p2[j]))
}
library(fields)
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
imagePlot(p1,p2,fc, xlab=expression(x[1]), ylab=expression(x[2]),col=cm.colors(12, rev = TRUE),main="Validation points",cex.main=3,cex.lab=2,cex.axis=2)
points(D,pch=16,col="blue",lwd=4,cex=2.5)
points(A,pch=2,col="darkred",lwd=4,cex=2.5)

boxplot(rmse,main="RMSE",cex.main=3,cex.axis=2)
abline(h=rmse0,col="darkred",lwd=3)
text(.75,rmse0+2,"sequential MaxPro",cex=2)
par(mfrow=c(1,1))

##################################
#Section 4.5.1 Constrained regions
##################################

rm(list=ls())
constraint <- function(x){
  c1 <- (x[1] - sqrt(50 * (x[2] - 0.52)^2 + 2) + 1)
  c2 <- (sqrt(120 * (x[2] - 0.48)^2 + 1) - 0.75 - x[1])
  c3 <- (0.65^2 - x[1]^2 - x[2]^2)
  return (c(c1,c2,c3))
}
# define constraint contour
x1 <- x2 <- matrix(NA, nrow = 3, ncol = 1001)
x2.seq <- seq(0,1,length.out = 1001)
x2[1,] <- x2.seq
x1[1,] <- sqrt(50 * (x2.seq - 0.52)^2 + 2) - 1
x2[2,] <- x2.seq
x1[2,] <- sqrt(120 * (x2.seq - 0.48)^2 + 1) - 0.75
x2[3,] <- x2.seq
x1[3,] <- sqrt(0.65^2 - x2.seq^2)
contour <- list(x1 = x1, x2 = x2)

#Figure 4.15 (plot feasible region)

plot(NULL, type = 'n', xlim = c(0,1), ylim = c(0,1),xlab = expression(x[1]), ylab = expression(x[2]),main="feasible region",cex.lab=2,cex.main=3,cex.axis=2)
for (i in 1:3) lines(contour$x1[i,], contour$x2[i,],lwd=2)
polygon(x = c(contour$x1[1,383:536], contour$x1[2,536:529], 
              contour$x1[3,529:412], contour$x1[2,412:383]),
        y = c(contour$x2[1,383:536], contour$x2[2,536:529], 
              contour$x2[3,529:412], contour$x2[2,412:383]),
        col = "red")


#Figure 4.16

library(lhs)
set.seed(1)
N=1e5
lhs.all <- randomLHS(N, 2)
lhs.gval <- t(apply(lhs.all, 1, constraint))
lhs.out.idx <- apply(lhs.gval, 1, function(x) return(any(x>0)))
lhs.feasible <- lhs.all[!lhs.out.idx,]
nrow(lhs.feasible) # number of feasible samples

layout(matrix(c(1:3), nrow = 1, byrow = T))
# all Latin Hypercube Samples
plot(lhs.all, col = "green", pch = 18, cex = 0.5, 
     xlim = c(0,1), ylim = c(0,1),xlab = expression(x[1]), ylab = expression(x[2]),main="initial points",cex.lab=2,cex.main=3,cex.axis=2)
for (i in 1:3) lines(contour$x1[i,], contour$x2[i,])
# only feasaible Latin Hypercube Samples
plot(lhs.feasible, col = "red", pch = 16, cex = 0.75, 
     xlim = c(0.3,0.8), ylim = c(0.35,0.55),xlab = expression(x[1]), ylab = expression(x[2]),main="candidate points",cex.lab=2,cex.main=3,cex.axis=2)
for (i in 1:3) lines(contour$x1[i,], contour$x2[i,])

library(MaxPro)
n=20
a=MaxProAugment(ExistDesign = lhs.feasible[1,], CandDesign = lhs.feasible[-1,],nNew = n-1)
D=a$Design
plot(D, col = "blue", pch = 16, cex = 2, 
     xlim = c(0.3,0.8), ylim = c(0.35,0.55),xlab = expression(x[1]), ylab = expression(x[2]),main="MaxPro design",cex.lab=2,cex.main=3,cex.axis=2)
for (i in 1:3) lines(contour$x1[i,], contour$x2[i,])
layout(matrix(c(1:1), nrow = 1, byrow = T))

##################################
#Section 4.5.2 Qualitative factors
##################################

#Figure 4.17

p=2
n=9
#Bad LHD
set.seed(1)
D0=cbind(rep(c(1,2,3),rep(3,3)),1:9,sample(1:9))
D0[,2:3]=(D0[,2:3]-.5)/9
#SLHD
library(SLHD)
a=maximinSLHD(3,3,2)
D1=a$StandDesign
#MaxPro
library(MaxPro)
set.seed(3)
ini=cbind(MaxPro(MaxProLHD(n,p)$Design)$Design,rep(1:3,c(3,3,3)))
a=MaxProQQ(InitialDesign = ini,p_nom = 1)
a$measure
D=a$Design

par(mar=c(5,5,4,2))
par(mfrow=c(1,3))
plot(D0[,2:3],bty="n",pch=D0[,1],col=D0[,1]+1,lwd=5,xlab=expression(x[1]),ylab=expression(x[2]),main="random LHD",cex=4,cex.lab=2,cex.main=3,cex.axis=2,xlim=c(0,1),ylim=c(0,1))
for(i in 1:n)
{
  segments(i/n,0,i/n,1,lty=3,lwd=2)
  segments(0,i/n,1,i/n,lty=3,lwd=2)
}
for(i in 0:(n/3)){
  segments(3*i/n,0,3*i/n,1,lwd=3)
  segments(0,3*i/n,1,3*i/n,lwd=3)
}

plot(D1[,2:3],bty="n",pch=D1[,1],col=D1[,1]+1,lwd=5,xlab=expression(x[1]),ylab=expression(x[2]),main="SLHD",cex=4,cex.lab=2,cex.main=3,cex.axis=2,xlim=c(0,1),ylim=c(0,1))
for(i in 1:n)
{
  segments(i/n,0,i/n,1,lty=3,lwd=2)
  segments(0,i/n,1,i/n,lty=3,lwd=2)
}
for(i in 0:(n/3)){
  segments(3*i/n,0,3*i/n,1,lwd=3)
  segments(0,3*i/n,1,3*i/n,lwd=3)
}

plot(D[,1:2],bty="n",pch=D[,3],col=D[,3]+1,lwd=5,xlab=expression(x[1]),ylab=expression(x[2]),main="MaxPro",cex=4,cex.lab=2,cex.main=3,cex.axis=2)
for(i in 0:n)
{
  segments(i/n,0,i/n,1,lty=3,lwd=2)
  segments(0,i/n,1,i/n,lty=3,lwd=2)
}
for(i in 0:(n/3)){
  segments(3*i/n,0,3*i/n,1,lwd=3)
  segments(0,3*i/n,1,3*i/n,lwd=3)
}
par(mfrow=c(1,1))

################################
#Section 4.5.3 Branching designs
################################

#Figure 4.19

library(SFDesign)
library(MaxPro)
for(k in 1:20){
  set.seed(k)
  #create an initial design for shared+nested factors
  shared_nested.ini=maxpro.optim(maxproLHD(30,6+2)$design)$design
  #augment the initial design with a column for the branching factor
  branch=rep(c(0,1),c(15,15))
  ini=cbind(shared_nested.ini,branch)
  #find the MaxProQQ design
  D0=MaxProQQ(ini,p_nom=1)$Design
  na_nested=matrix(NA,nrow=15,ncol=2)
  val=matrix(Inf,nrow=8,ncol=8)
  D1=D0[1:15,-9]
  D2=D0[16:30,-9]
  #find the optimal pair of columns for the nested factors
  for(i in 1:7){
    for(j in (i+1):8){
      #make two columns of D2 as NA
      B=D2
      B[,c(i,j)]=na_nested
      #rescale and optimize the corresponsing columns of D1
      A=D1
      A[,c(i,j)]=(apply(D1[,c(i,j)],2,rank)-.5)/15
      A[,c(i,j)]=maxpro.optim(A[,c(i,j)])$design
      #combine  them to obtain the branching design
      D=cbind(rbind(A,B), branch)
      val[i,j]=MaxProMeasure(D[,-c(i,j)])+MaxProMeasure(D[1:15,-9])+MaxProMeasure(D[16:30,-c(9,i,j)])
    }
  }
  print(c(k,min(val)))
} 
#find the best initial design (best seed)
opt.seed=12
#construct the optimal branching design
set.seed(opt.seed)
shared_nested.ini=maxpro.optim(maxproLHD(30,6+2)$design)$design
branch=rep(c(0,1),c(15,15))
ini=cbind(shared_nested.ini,branch)
D0=MaxProQQ(ini,p_nom=1)$Design
na_nested=matrix(NA,nrow=15,ncol=2)
val=matrix(Inf,nrow=8,ncol=8)
D1=D0[1:15,-9]
D2=D0[16:30,-9]
for(i in 1:7){
  for(j in (i+1):8){
    B=D2
    B[,c(i,j)]=na_nested
    A=D1
    A[,c(i,j)]=(apply(D1[,c(i,j)],2,rank)-.5)/15
    A[,c(i,j)]=maxpro.optim(A[,c(i,j)])$design
    D=cbind(rbind(A,B), branch)
    val[i,j]=MaxProMeasure(D[1:15,-9])+MaxProMeasure(D[16:30,-c(9,i,j)])
  }
}
ind=c(which(val == min(val), arr.ind = TRUE))
B=D2
B[,ind]=na_nested
A=D1
A[,ind]=(apply(D1[,ind],2,rank)-.5)/15
A[,ind]=maxpro.optim(A[,ind])$design
D=cbind(rbind(A,B), branch)
print(c(MaxProMeasure(D[,-ind]),MaxProMeasure(D[1:15,-9]),MaxProMeasure(D[16:30,-c(9,ind)])))

#re-order the columns
temp=D
temp[,1]=D[,9]
temp[,2:3]=D[,ind]
temp[,4:9]=D[,-c(9,ind)]
D=temp
colnames(D)=c("branch", "nested1","nested2", paste("shared",sep="", 1:6))
pairs(D, pch=16, col="blue",cex.labels = 1.1)


#############################
#Section 4.5.4 Nested designs
#############################

#Figure 4.21

rm(list=ls())
library(MaxPro)
library(SFDesign)
n1=100
n2=10
p=2

set.seed(1)
D1=maxpro.optim(maxproLHD(n1,p)$design)$design
maxpro.crit(D1)

D2=maxpro.remove(D1,n.remove = n1-n2,delta = 1/n1^2)
maxpro.crit(D2)

par(mar=c(5,5,4,2))
plot(D1,xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",lwd=4, cex=2,main="MaxPro-based nested design",xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=2,cex.main=2.5)
points(D2,pch=2,col=2,lwd=4,cex=2.5)

###################################
#Section 4.6 Minimum energy designs
###################################

#Figure 4.22

rm(list=ls())
library(mined)
library(zipcodeR)

dat=na.omit(zip_code_db[,c(9,8,14)])
ind=(dat[,2]<52 & dat[,2]>25)
dat=dat[ind,]
CAND=as.matrix(dat[,1:2])
y=dat[,3]
sizecode=.5+5*sqrt(y)/max(sqrt(y))
D=SelectMinED(candidates = CAND,candlf = log(y),n=200,s=1)$points

library(maximin)
p=2
D0=matrix(CAND[sample(1:dim(CAND)[1],1),],nrow=1,ncol=p)
D2=CAND[maximin.cand(199,Xcand = CAND,Xorig = D0)$inds,]
D2=rbind(D0,D2)

sizecode=.25+log(y+1)/max(log(y+1))
par(mfrow=c(1,2))
plot(CAND,cex=sizecode,pch=16,col="pink",xlab="",ylab="",main="maximin design",cex.main=3)
points(D2,pch=16,col="blue",cex=2)
plot(CAND,cex=sizecode,pch=16,col="pink",xlab="",ylab="",main="minimum energy design",cex.main=3)
points(D,pch=16,col="blue",cex=2)
par(mfrow=c(1,1))


###################################
#Section 4.6.1 Bayesian computation
###################################

#Figure 4.24

rm(list=ls())
library(mined)
library(MaxPro)
logf <- function(para)
{
  l1 <- -40
  u1 <- 40
  l2 <- -25
  u2 <- 10
  x1 <- l1 + (u1 - l1) * para[1]
  x2 <- l2 + (u2 - l2) * para[2]
  val <- -.5 * (x1 ^2 / 100 + (x2+ .03 * x1^2 -3)^2)
  return(val)
}
N.plot=300
p1 <- seq(0, 1, length.out = N.plot)
p2 <- seq(0, 1, length.out = N.plot)
fc <- matrix(0.0, N.plot, N.plot)
for(i in 1:N.plot)
{
  for(j in 1:N.plot)
  {
    fc[i, j] = exp(logf(c(p1[i], p2[j])))
  }
}
set.seed(8)
p=2
n=20
ini=MaxPro(MaxProLHD(n,p)$Design)$Design
res=mined(ini,logf,K_iter = 5)
D=res$points
cand=res$cand

library(fields)
par(mar=c(5,5,4,5))
par(mfrow=c(1,2))
imagePlot(p1,p2,fc, xlab=expression(theta[1]), ylab=expression(theta[2]),col=cm.colors(5),main="initial MaxPro design",cex.main=2.5,cex.lab=2,cex.axis=2)
points(ini,pch=16,col="darkred",cex=2)
imagePlot(p1,p2,fc, xlab=expression(theta[1]), ylab=expression(theta[2]),col=cm.colors(5),main="MED points with annealing",cex.main=2.5,cex.lab=2,cex.axis=2)
ind=match(ini[,1],cand[,1])
points(ini,pch=16,col="darkred",cex=2)
points(cand[-ind,],pch=16,col="blue",cex=2)
par(mfrow=c(1,1))

#Figure 4.25

nlogh=function(para) -logf(para)
h=function(para) exp(-nlogh(para))
target=function(para) .5*logf(para)
set.seed(8)
p=2
n=20
ini=MaxPro(MaxProLHD(n,p)$Design)$Design
res=mined(ini,target,K_iter = 5)
D=res$points
cand=res$cand
nu=cand
hev=apply(nu,1,h)

G.m=function(para)
{
  A=nu%*%diag(sqrt(1/para))
  A=as.matrix(dist(A, diag=T, upper=T))
  return(exp(-.5*A^2))
}
exp1=function(e) exp(-.5*sum(e^2/sigma2))
g=function(gamma)
{
  A=t(t(nu)-gamma)
  vec=apply(A,1,exp1)
  return(vec)
}
wmscv.sigma=function(para)
{
  G=G.m(para)
  Ginv=solve(G+.0001*diag(dim(G)[1]))
  di=diag(Ginv)
  cv=c(Ginv%*%sqrt(hev))/di
  val=mean(di*cv^2)
  return(log(val))
}
ini=(median(dist(nu)))^2/2
sigma2=optim(rep(ini,p), wmscv.sigma,lower=rep(ini/100,p), upper=rep(100*ini,p), method="L-BFGS-B")$par
sigma2
Sigma=diag(sigma2)
s=sqrt(diag(Sigma))
B=nu%*%diag(1/s)
B=as.matrix(dist(B, diag=T, upper=T))
G=exp(-.5*B^2)

Ginv=solve(G+.0001*diag(dim(G)[1]))
coef=Ginv%*%sqrt(hev)
G2=exp(-.25*B^2)
denom=pi*prod(s)*as.numeric(t(coef)%*%G2%*%coef)
Sigma.inv=diag(1/sigma2)
m=dim(nu)[1]
M=array(0,dim = c(m,m,p))
for(k in 1:p) M[,,k]=outer(nu[,k],nu[,k],"+")/2
d=outer(c(coef),c(coef),"*")*G2
denk=function(theta,k){
  num=sum(d*dnorm(theta,M[,,k],s[k]/sqrt(2)))
  return(num/sum(d))
}

den=matrix(0,nrow=N.plot,ncol=p)
v=cbind(p1,p2)
for(k in 1:p) den[,k]=apply(cbind(v[,k]),1,denk,k=k)

#exact density by numerical integration
library(cubature)
denom=adaptIntegrate(h, c(0,0),c(1,1))$int
h1=function(theta1) apply(cbind(theta1,theta2),1,h)
val=p2
for(i in 1:N.plot)
{
  theta2=p2[i]
  val[i]=integrate(h1,0,1)$val
}
exact=matrix(0,nrow=N.plot,ncol=p)
exact[,2]=val/denom

h2=function(theta2) apply(cbind(theta1,theta2),1,h)
val=p1
for(i in 1:N.plot)
{
  theta1=p1[i]
  val[i]=integrate(h2,0,1)$val
}
exact[,1]=val/denom

#MCMC
library(adaptMCMC)
image(p1,p2,fc,col=cm.colors(5),xlab=expression(theta[1]),ylab=expression(theta[2]))
s=(2.4/sqrt(2))^2*diag(p)#specify covariance matrix of the jumping distribution
out=MCMC(logf,n=10000,init=rep(.5,p),scale=s,adapt=T,acc.rate=.05)
theta=out$samples
points(theta)

par(mfrow=c(1,2))
fac=c(1.2,1)
k=1
u=max(c(den[,k],exact[,k]))*fac[k]
plot(density(theta[,k]),type="l",col=2,lwd=3,ylim=c(0,u),lty=1, main = "", xlab=expression(theta[1]),cex.lab=2,cex.axis=2)
lines(v[,k],exact[,k],lwd=4,col=3,lty=1)
lines(v[,k],den[,k],col="blue",lwd=3)
legend("topleft",col=c(3,"blue",2),lwd=c(4,3,3),legend=c("True", "DoIt","MCMC"), bty="n",cex=2)
k=2
u=max(c(den[,k],exact[,k]))*fac[k]
plot(density(theta[,k]),type="l",col=2,lwd=3,ylim=c(0,u),lty=1, main = "", xlab=expression(theta[2]),cex.lab=2,cex.axis=2)
lines(v[,k],exact[,k],lwd=4,col=3,lty=1)
lines(v[,k],den[,k],col="blue",lwd=3)
legend("topleft",col=c(3,"blue",2),lwd=c(4,3,3),legend=c("True", "DoIt","MCMC"), bty="n",cex=2)
par(mfrow=c(1,1))

#######################################
#Section 4.6.2 Candidate set generation
#######################################
#https://github.com/BillHuang01/CoMinED

#Figure 4.26

rm(list=ls())
constraint <- function(x){
  c1 <- (x[1] - sqrt(50 * (x[2] - 0.52)^2 + 2) + 1)
  c2 <- (sqrt(120 * (x[2] - 0.48)^2 + 1) - 0.75 - x[1])
  c3 <- (0.65^2 - x[1]^2 - x[2]^2)
  return (c(c1,c2,c3))
}
# define constraint contour
x1 <- x2 <- matrix(NA, nrow = 3, ncol = 1001)
x2.seq <- seq(0,1,length.out = 1001)
x2[1,] <- x2.seq
x1[1,] <- sqrt(50 * (x2.seq - 0.52)^2 + 2) - 1
x2[2,] <- x2.seq
x1[2,] <- sqrt(120 * (x2.seq - 0.48)^2 + 1) - 0.75
x2[3,] <- x2.seq
x1[3,] <- sqrt(0.65^2 - x2.seq^2)
contour <- list(x1 = x1, x2 = x2)
# plot feasible region
plot(NULL, type = 'n', xlim = c(0,1), ylim = c(0,1),xlab = expression(x[1]), ylab = expression(x[2]),main="feasible region",cex.lab=2,cex.main=3,cex.axis=2)
for (i in 1:3) lines(contour$x1[i,], contour$x2[i,],lwd=2)
polygon(x = c(contour$x1[1,383:536], contour$x1[2,536:529], 
              contour$x1[3,529:412], contour$x1[2,412:383]),
        y = c(contour$x2[1,383:536], contour$x2[2,536:529], 
              contour$x2[3,529:412], contour$x2[2,412:383]),
        col = "pink")
# log distance function
logdis <- function(x, s=2){
  if (s > 0) length(x)/s*log(sum(abs(x)^s)) else sum(log(abs(x)))
}

# mined by one-point-at-a-time greedy algorithm
mined.seq <- function(n, cand, cand.lf, s = 2, return.obj = F){
  N <- nrow(cand)
  if (n > N) stop(message("Not enough candidate points!"))
  p <- ncol(cand)
  idx <- which.max(cand.lf)
  xi <- c(idx)
  val <- NULL
  for (i in 2:n){
    cand.dist <- cand - rep(1,N) %*% t(cand[idx,])
    val <- cbind(val, 0.5*cand.lf[idx] + 0.5*cand.lf + apply(cand.dist, 1, logdis, s=s))
    idx <- which.max(apply(val, 1, min))
    xi <- c(xi, idx)
  }
  if (return.obj){
    val <- val[xi,]
    val <- cbind(val, c(val[n,],-Inf))
    obj <- min(val[upper.tri(val)])
    return (list(idx = xi, obj = obj))
  }
  return (xi)
}
comined <- function(n, p, tau, constraint, n.aug, auto.scale = F, 
                    s = 2, visualization = F, visualization.params = list()){
  # initialize sample
  samp <- Lattice(n*n.aug, p)
  min.dist <- min(dist(samp))
  min.ddist <- Inf # dimensional minimum
  for (i in 1:p){
    ddist <- min(dist(samp[,i]))
    if (ddist < min.ddist) min.ddist <- ddist
  }
  samp.gval <- matrix(t(apply(samp,1,constraint)), nrow=nrow(samp))
  # scale parameter if auto-scale is on
  scale <- rep(1, ncol(samp.gval))
  if (auto.scale) scale <- apply(samp.gval,2,mad,center=0)
  # find first set of minimum energy design points
  samp.lf <- apply(samp.gval,1,function(x) sum(pnorm(-tau[2]*x/scale,log.p=T)))
  med.op <- mined.seq(n, samp, samp.lf, s = s, return.obj = TRUE)
  samp.med <- samp[med.op$idx,]
  if (visualization){
    vparams <- list(unit.scale = FALSE, contour = NULL)
    vparams[(nm<-names(visualization.params))] <- visualization.params
    xlim <- ylim <- NULL
    if (vparams$unit.scale) xlim <- ylim <- c(0,1)
    plot(samp, col = "green", pch = 18, cex = 0.75, 
         xlim = xlim, ylim = ylim, xlab = "", ylab = "")
    points(samp.med, col = "red", pch = 16, cex = 1)
    if (!is.null(vparams$contour)){
      contours <- vparams$contour
      for (l in 1:nrow(contour$x1)){
        lines(contour$x1[l,], contour$x2[l,], lty = 2)
      }
    }
  }
  
  for (k in 3:length(tau)){
    min.dist <- min.dist / 2
    min.ddist <- min.ddist / 2
    no.decimal <- attr(regexpr("(?<=\\.)0+", format(min.ddist, scientific = F), 
                               perl = TRUE), "match.length") + 1
    # sample augmentation
    samp.med.dist <- as.matrix(dist(samp.med))
    samp.aug <- NULL
    for (i in 1:n){
      nn.idx <- order(samp.med.dist[i,])[2:(n.aug+1)]
      samp.aug <- rbind(samp.aug, 0.5*(samp.med[nn.idx,]+rep(1,n.aug)%*%t(samp.med[i,])))
      samp.aug <- rbind(samp.aug, 0.5*(3*samp.med[nn.idx,]-rep(1,n.aug)%*%t(samp.med[i,])))
    }
    # remove repeated samples
    samp.aug.rep <- round(samp.aug, digits = no.decimal)
    samp.aug.dp <- duplicated(samp.aug.rep)
    samp.aug <- samp.aug[!samp.aug.dp,]
    # remove samples outside of boundary
    samp.aug.out <- apply(samp.aug, 1, function(x) (any(x<0)|any(x>1)))
    samp.aug <- samp.aug[!samp.aug.out,]
    # remove nearby point due to numerical issue
    no.aug <- nrow(samp.aug)
    samp.rep <- rbind(samp.aug, samp)
    samp.rep <- round(samp.rep, digits = no.decimal)
    samp.rep.dp <- duplicated(samp.rep, fromLast = TRUE)
    samp.aug <- samp.aug[!(samp.rep.dp[1:no.aug]),]
    if (nrow(samp.aug) > 0){
      # compute constraint value for augmented samples
      samp.aug.gval <- matrix(t(apply(samp.aug,1,constraint)), 
                              nrow=nrow(samp.aug))
      # combine augmented samples with original samples
      samp <- rbind(samp, samp.aug)
      samp.gval <- rbind(samp.gval, samp.aug.gval)
    }
    # find minimum energy design points
    if (auto.scale) scale <- apply(samp.gval,2,mad,center=0)
    samp.lf <- apply(samp.gval,1,function(x) sum(pnorm(-tau[k]*x/scale,log.p=T)))
    if (s == 0){
      hp <- floor(p/2)
      nl <- sqrt((min.dist^2-min.ddist^2*hp)/(p-hp))
      min.logdis <- logdis(c(rep(min.ddist,hp),rep(nl,(p-hp))),s=s)
    } else {
      min.logdis <- logdis(rep(min.dist/sqrt(p),p),s=s)
    }
    samp.lf.cv <- sort(samp.lf, decreasing = T)[n] + 2.5 * min.logdis
    samp.cand.idx <- (samp.lf > samp.lf.cv)
    samp.cand <- samp[samp.cand.idx,]
    samp.cand.lf <- samp.lf[samp.cand.idx]
    med.op <- mined.seq(n, samp.cand, samp.cand.lf, s = s, return.obj = TRUE)
    samp.med <- samp.cand[med.op$idx,]
    if (visualization){
      plot(samp.cand, col = "green", pch = 18, cex = 0.75, 
           xlim = xlim, ylim = ylim, xlab = "", ylab = "")
      points(samp.med, col = "red", pch = 16, cex = 1)
      if (!is.null(vparams$contour)){
        contour <- vparams$contour
        for (l in 1:nrow(contour$x1)){
          lines(contour$x1[l,], contour$x2[l,], lty = 2)
        }
      }
    }
    
  }
  
  samp.med.lf <- samp.cand.lf[med.op$idx]
  feasible.idx <- !(apply(samp.gval, 1, function(x) any(x>0)))
  
  return (list(med = samp.med, 
               med.lf = samp.med.lf,
               cand = samp, 
               cand.lf = samp.lf,
               cand.min.dist = min.dist,
               cand.min.ddist = min.ddist,
               feasible.idx = feasible.idx))
  
}

# Run CoMinED (the intermediate plots it draws are not included in the book)
set.seed(20210329)
layout(matrix(c(1:8), nrow = 2, byrow = T))
tau <- c(0,exp(c(1:7)),1e6)
comined.output <- comined(n = 53, p = 2, tau = tau, constraint = constraint,
                          n.aug = 5, auto.scale = F, s = 2, visualization = T,
                          visualization.params = list(unit.scale=TRUE,contour=contour))

par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
cand=comined.output$cand
plot(cand, col = "blue", pch = 16, cex = 0.75, 
     xlim = c(0,1), ylim = c(0,1),xlab = expression(x[1]), ylab = expression(x[2]),main="full candidate set",cex.lab=2,cex.main=3,cex.axis=2)
for (i in 1:3) lines(contour$x1[i,], contour$x2[i,],col=2,lwd=2)

cand.feasible=comined.output$cand[comined.output$feasible.idx,]
plot(cand.feasible, col = "green", pch = 16, cex = 0.75, 
     xlim = c(0.3,0.8), ylim = c(0.35,0.55), xlab = expression(x[1]), ylab = expression(x[2]),main="feasible candidate set",cex.lab=2,cex.main=3,cex.axis=2)
for (i in 1:3) lines(contour$x1[i,], contour$x2[i,],col=2,lwd=2)
par(mfrow=c(1,1))


############################################
#Section 4.6.3 Continuous fidelity parameter
############################################

#Figure 4.27

rm(list=ls())
set.seed(1)
library(MaxPro)
library(SFDesign)
library(rkriging)
f=function(u){
  x=u[1:2]
  M=1/u[3]
  term1=(1-exp(-1/(2*x[2])))*(2300*x[1]^3+1900*x[1]^2+2092*x[1]+60)/(100*x[1]^3+500*x[1]^2+4*x[1]+20)
  term2=16*exp(-1.4*x[1])*cos(3.5*pi*x[2])
  val=term1+term2*u[3]
  return(val)
}
p=2
Mt=2
Mb=16
lam=Mb/Mt
nb=50
nu=(nb*lam-1)/(lam*(nb-1))
n=1+round(log(lam)/log(nu))
u=((1:n)-1)/(n-1)
M=Mt*(Mb/Mt)^u
h=1/M
h=sort(h/max(h))

D=maxpro.optim(maxproLHD(n,p+1)$design)$design
D=D[order(D[,3]),]
D[,3]=h
MaxProMeasure(D)
a=MaxProQQ(D)
a$measure
D=a$Design
y=apply(D,1,f)

library(car)
scatterplotMatrix(D,diagonal=list(method ="histogram",breaks=15),smooth=FALSE, regLine = FALSE,pch=16,var.labels = c(expression(x[1]),expression(x[2]),"t"),cex.axis = 2)


#Generate testing points
test=MaxProAugment(ExistDesign = D[,1:p],CandDesign = CandPoints(100*n,p),nNew = 100)$Design[-(1:n),]
true=apply(cbind(test,0),1,f)

#Fit the model in (4.58)

OKfit=function(D,y, ini,nug=10^(-6)){
  n=length(y)
  if(!is.matrix(D)) D=matrix(D,ncol=1)
  p=dim(D)[2]-1
  
  r=function(x,theta){
    theta0=theta[1:p]
    theta1=theta[(p+1):(2*p)]
    gam=theta[2*p+1]
    lam=theta[2*p+2]
    if(p>1) {
      A=sweep(D[,1:p],2,x[1:p],"-")
      g0=exp(-apply((A%*%diag(1/theta0))^2,1,sum))
      g1=exp(-apply((A%*%diag(1/theta1))^2,1,sum))
      b=exp(-((D[,p+1]-x[p+1])/gam)^2)-exp(-((D[,p+1])/gam)^2)-exp(-((x[p+1])/gam)^2)+1
      g=g0+lam*g1*b
    }else {
      g0=exp(-((D[,1]-x[1])/theta0)^2)
      g1=exp(-((D[,1]-x[1])/theta1)^2)
      b=pmin(D[,2],x[2])
      g=g0+lam*g1*b
    }
    return(g)
  }
  
  ML=function(theta)
  {
    R=apply(D,1,r,theta)+nug*diag(n)
    L=chol(R)
    a=forwardsolve(t(L),one)
    b=forwardsolve(t(L),y)
    mu=sum(a*b)/sum(a^2)
    nu2=1/n*sum((b-mu*a)^2)
    nu2=max(nu2,10^(-15))
    val=log(nu2)+2*mean(log(diag(L)))
    return(val)
  }
  one=rep(1,n)
  a.opt=optim(ini,ML,lower=ini/10, upper = 10*ini, method = "L-BFGS-B",control=list("maxit"=100))
  theta=a.opt$par
  R=apply(D,1,r,theta)+nug*diag(n)
  L=chol(R)
  a=forwardsolve(t(L),one)
  b=forwardsolve(t(L),y)
  mu=sum(a*b)/sum(a^2)
  nu2=1/n*sum((b-mu*a)^2)
  return(list(mu=mu,nu2=nu2, theta=theta, L=L,a=a,b=b, D=D))
}

OKpredict=function(obj,new){
  mu=obj$mu
  nu2=obj$nu2
  theta=obj$theta
  L=obj$L
  a=obj$a
  b=obj$b
  D=obj$D
  r=function(x,theta){
    theta0=theta[1:p]
    theta1=theta[(p+1):(2*p)]
    gam=theta[2*p+1]
    lam=theta[2*p+2]
    if(p>1) {
      A=sweep(D[,1:p],2,x[1:p],"-")
      g0=exp(-apply((A%*%diag(1/theta0))^2,1,sum))
      g1=exp(-apply((A%*%diag(1/theta1))^2,1,sum))
      b=exp(-((D[,p+1]-x[p+1])/gam)^2)-exp(-((D[,p+1])/gam)^2)-exp(-((x[p+1])/gam)^2)+1
      g=g0+lam*g1*b
    }else {
      g0=exp(-((D[,1]-x[1])/theta0)^2)
      g1=exp(-((D[,1]-x[1])/theta1)^2)
      b=pmin(D[,2],x[2])
      g=g0+lam*g1*b
    }
    return(g)
  }
  ok=function(x,para) {
    d=forwardsolve(t(L),r(x,para))
    t1=mu+sum(d*(b-mu*a))
    t2=sqrt(nu2)*sqrt(max(1-sum(d^2)+(1-sum(d*a))^2/sum(a^2),0))
    return(c(t1,t2))
  }
  pred=apply(cbind(new),1,ok, para=theta)
  return(list(mean=pred[1,],sd=pred[2,]))
}
a=OKfit(D,y, ini=rep(1,2*p+2))
pred=OKpredict(a,cbind(test,0))$mean
sqrt(mean((pred-true)^2))

#Single-fidelity simulation
D1=maxpro.optim(maxproLHD(nb,p)$design)$design
y1=apply(cbind(D1,1/Mb),1,f)
a1<-Fit.Kriging(D1,y1,kernel.parameters = list(type="Gaussian"))
pred1=Predict.Kriging(a1,test)$mean
sqrt(mean((pred1-true)^2))
