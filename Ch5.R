################
#Chapter 5: Representative Points
#Copyright: V. Roshan Joseph
#Please cite: V. R. Joseph (2026). Experimental Design for Data Science and Engineering. CRC Press.
######################################

#######################################
#Section 5.1: Uniform Design
######################################

#################
#5.1.1 Uniformity
#################

#Figure 5.1

rm(list=ls())
set.seed(1)
n=7
x=sort(runif(n))
d=matrix(0,nrow=n,ncol=2)
d[,1]=x
d[,2]=(1:n)/n
par(mar=c(5,5,4,2))
par(mfrow=c(1,2))
plot(d,xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",cex=3,xlab="x",ylab = "distribution function",main="random",cex.lab=2,cex.main=3,cex.axis=2)
segments(0,0,x[1],0,col=2,lwd=2)
for(i in 1:(n-1)) segments(x[i],i/n,x[i+1],i/n,col=2,lwd=2)
segments(x[n],1,1,1,col=2,lwd=2)
abline(0,1,col=3,lwd=2,lty=2)


x=((1:n)-.5)/n
d=matrix(0,nrow=n,ncol=2)
d[,1]=x
d[,2]=(1:n)/n
plot(d,xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",cex=3,xlab="x",ylab = "distribution function",main="optimal",cex.lab=2,cex.main=3,cex.axis=2)
segments(0,0,x[1],0,col=2,lwd=2)
for(i in 1:(n-1)) segments(x[i],i/n,x[i+1],i/n,col=2,lwd=2)
segments(x[n],1,1,1,col=2,lwd=2)
abline(0,1,col=3,lwd=2,lty=2)
par(mfrow=c(1,1))

#Figure 5.2

library(spacefillr)
library(SFDesign)
p=2
n=20
set.seed(6)
D = uniform.optim(uniformLHD(n,p)$design)$design
R=matrix(runif(n*2),ncol=2)
S=generate_sobol_set(n,p,seed=sample(1:10000,1))
c(uniform.crit(R),uniform.crit(S),uniform.crit(D))
par(mar=c(5,5,4,2))
par(mfrow=c(1,3))
plot(R,pch=16,col=4,cex=4,main="Monte Carlo",cex.lab=3,cex.main=4,cex.axis=2,xlab=expression(x[1]),ylab=expression(x[2]))
plot(S,pch=16,col=4,cex=4,main="Sobol Sequence",cex.lab=3,cex.main=4,cex.axis=2,xlab=expression(x[1]),ylab=expression(x[2]))
plot(D,pch=16,col=4,cex=4,main="Uniform Design",cex.lab=3,cex.main=4,cex.axis=2,xlab=expression(x[1]),ylab=expression(x[2]))
par(mfrow=c(1,1))


#Figure 5.3

h=function(x) exp(-sum((x-.5)^2))
mean(apply(S,1,h))
mean(apply(D,1,h))
true=pi*(pnorm(.5*sqrt(2))-pnorm(-.5*sqrt(2)))^2

n=seq(10,100,by=10)
randw=sobw=uniw=matrix(0,nrow=10,ncol=30)
rand=sob=uni=matrix(0,nrow=10,ncol=30)
set.seed(1)#this simulation will take time!!!
for(j in 1:10){
  for(i in 1:30){
    R=matrix(runif(n[j]*2),ncol=2)
    S=generate_sobol_set(n[j],2,seed=sample(1:10000,1))
    D = uniform.optim(uniformLHD(n[j],p)$design)$design
    
    randw[j,i]=uniform.crit(R)
    sobw[j,i]=uniform.crit(S)
    uniw[j,i]=uniform.crit(D)
    
    rand[j,i]=mean(apply(R,1,h))
    sob[j,i]=mean(apply(S,1,h))
    uni[j,i]=mean(apply(D,1,h))
  }
}
design=expand.grid(n=seq(10,100,by=10),method=c("MC","Sobol","Uniform"))
integration=cbind(design,rbind(rand,sob,uni))
wd=cbind(design,rbind(randw,sobw,uniw))
names(integration)[3:32]=paste("V",sep="",1:30)
names(wd)[3:32]=paste("V",sep="",1:30)

library(tidyverse)
library(ggplot2)
data <- wd %>% 
  pivot_longer(V1:V30, names_to = "run", values_to = "value") %>% 
  group_by(method, n) %>% 
  summarize(l = quantile(value, 0.05),
            m = quantile(value, 0.5),
            u = quantile(value, 0.95))


ggplot(data) +
  geom_line(aes(x = n, y = m, color = method), linewidth = .75) +
  geom_point(aes(x = n, y = m, shape = method), size = 1) +
  geom_ribbon(aes(x = n, ymin = l, ymax = u, fill = method), alpha = .5) +
  scale_color_discrete(name = "method") +
  scale_shape_discrete(name = "method") +
  scale_fill_discrete(name = "method") +
  scale_x_continuous(breaks = seq(10,100,by=10)) +
  labs(x = "n", y = "discrepancy",title = "Discrepancy") +
  theme_classic() +
  theme(
    legend.title = element_blank(),plot.title = element_text(size=20, hjust = 0.5)
  )
ggsave("wd.pdf", width = 4, height = 3)

data <- integration %>% 
  pivot_longer(V1:V30, names_to = "run", values_to = "value") %>% 
  group_by(method, n) %>% 
  summarize(l = quantile(value, 0.05),
            m = quantile(value, 0.5),
            u = quantile(value, 0.95))


ggplot(data) +
  geom_line(aes(x = n, y = m, color = method), linewidth = .75) +
  geom_point(aes(x = n, y = m, shape = method), size = 1) +
  geom_ribbon(aes(x = n, ymin = l, ymax = u, fill = method), alpha = .5) +
  scale_color_discrete(name = "method") +
  scale_shape_discrete(name = "method") +
  scale_fill_discrete(name = "method") +
  scale_x_continuous(breaks = seq(10,100,by=10)) +
  labs(x = "n",y = "integral",title = "Integration") +
  theme_classic() +
  theme(
    legend.title = element_blank(),plot.title = element_text(size=20, hjust = 0.5)
  )
ggsave("integration.pdf", width = 4, height = 3)

#####################################
#Section 5.2 Nonuniform distributions
#####################################

######################################
#Section 5.2.1 Transformation methods
######################################

#Figure 5.4

p=2
n=50
set.seed(6)
library(SFDesign)
d = uniform.optim(uniformLHD(n, p)$design)$design
par(mfrow=c(1,2))
D=qnorm(d)
f=function(x) exp(-.5*sum(x^2))/(2*pi)
N.plot=300
p1 <- seq(-3, 3, length.out = N.plot)
p2 <- seq(-3, 3, length.out = N.plot)
fc=matrix(apply(expand.grid(p1,p2),1,f),N.plot,N.plot)
image(p1,p2,fc, xlab=expression(x[1]), ylab=expression(x[2]),col=cm.colors(5),main="Normal distribution",cex.main=3,cex.lab=2,cex.axis=2)
contour(p1,p2,fc,add=TRUE,nlevels = 10)
points(D,pch=16,col="blue",cex=3)

D=qgamma(d,2,1)
f=function(x) prod(dgamma(x,2,1))
N.plot=300
p1 <- seq(0, 6, length.out = N.plot)
p2 <- seq(0, 6, length.out = N.plot)
fc=matrix(apply(expand.grid(p1,p2),1,f),N.plot,N.plot)
image(p1,p2,fc, xlab=expression(x[1]), ylab=expression(x[2]),col=cm.colors(5),main="Gamma distribution",cex.main=3,cex.lab=2,cex.axis=2)
contour(p1,p2,fc,add=TRUE,nlevels = 10)
points(D,pch=16,col="blue",cex=3)
par(mfrow=c(1,1))


#Figure 5.5 (Rosenblatt transformation)
library(mnormt)
x1=qnorm(d[,1])
x2=qnorm(d[,2],mean=.5*x1,sd=sqrt(1-.5^2))
D=cbind(x1,x2)
f=function(x) dmnorm(x,mean=c(0,0),varcov = matrix(c(1,.5,.5,1),2,2))
N.plot=300
p1 <- seq(-3, 3, length.out = N.plot)
p2 <- seq(-3, 3, length.out = N.plot)
fc=matrix(apply(expand.grid(p1,p2),1,f),N.plot,N.plot)
image(p1,p2,fc, xlab=expression(x[1]), ylab=expression(x[2]),col=cm.colors(5),main="Correlated normal distribution",cex.main=2,cex.lab=2,cex.axis=2)
contour(p1,p2,fc,add=TRUE,nlevels = 10)
points(D,pch=16,col="blue",cex=3)

##############################
#Section 5.2.2 Support points
#############################

#Figure 5.6

set.seed(1)
library(support)
library(mnormt)
par(mfrow=c(1,2))
D=sp(n,p,dist.str = rep("normal",p))$sp
f=function(x) exp(-.5*sum(x^2))/(2*pi)
N.plot=300
p1 <- seq(-3, 3, length.out = N.plot)
p2 <- seq(-3, 3, length.out = N.plot)
fc=matrix(apply(expand.grid(p1,p2),1,f),N.plot,N.plot)
image(p1,p2,fc, xlab=expression(x[1]), ylab=expression(x[2]),col=cm.colors(5),main="Normal distribution",cex.main=3,cex.lab=2,cex.axis=2)
contour(p1,p2,fc,add=TRUE,nlevels = 10)
points(D,pch=16,col="blue",cex=3)

dist.param <- vector("list",p)
for (l in 1:p){
  dist.param[[l]] <- c(2,1)
}
D=sp(n,p,dist.str = rep("gamma",p),dist.param = dist.param)$sp
f=function(x) prod(dgamma(x,2,1))
N.plot=300
p1 <- seq(0, 6, length.out = N.plot)
p2 <- seq(0, 6, length.out = N.plot)
fc=matrix(apply(expand.grid(p1,p2),1,f),N.plot,N.plot)
image(p1,p2,fc, xlab=expression(x[1]), ylab=expression(x[2]),col=cm.colors(5),main="Gamma distribution",cex.main=3,cex.lab=2,cex.axis=2)
contour(p1,p2,fc,add=TRUE,nlevels = 10)
points(D,pch=16,col="blue",cex=3)
par(mfrow=c(1,1))

#Figure 5.7

set.seed(1)
par(mfrow=c(1,2))
N=10000
X=rmnorm(N,mean=c(0,0),varcov = matrix(c(1,.5,.5,1),2,2))
D=sp(n,p,dist.samp = X)$sp
f=function(x) dmnorm(x,mean=c(0,0),varcov = matrix(c(1,.5,.5,1),2,2))
N.plot=300
p1 <- seq(-3, 3, length.out = N.plot)
p2 <- seq(-3, 3, length.out = N.plot)
fc=matrix(apply(expand.grid(p1,p2),1,f),N.plot,N.plot)
image(p1,p2,fc, xlab=expression(x[1]), ylab=expression(x[2]),col=cm.colors(5),main="Correlated normal distribution",cex.main=2.5,cex.lab=2,cex.axis=2)
contour(p1,p2,fc,add=TRUE,nlevels = 10)
points(D,pch=16,col="blue",cex=3)

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
fc=matrix(exp(apply(expand.grid(p1,p2),1,logf)),N.plot,N.plot)
library(adaptMCMC)
a=MCMC(logf,n=100000, init=c(.5,.8), acc.rate = 0.4)
X=a$samples
D=sp(n,p,dist.samp = X)$sp
image(p1,p2,fc, xlab=expression(theta[1]), ylab=expression(theta[2]),col=cm.colors(5),main="Banana-shaped density",cex.main=2.5,cex.lab=2,cex.axis=2)
contour(p1,p2,fc,add=TRUE,nlevels = 10)
points(D,pch=16,col="blue",cex=2)
par(mfrow=c(1,1))

######################################
#Section 5.2.3 Uncertainty propagation
######################################

#Figure 5.8

library(support)
f=function(x) {
  val=2*pi*x[3]*(x[4]-x[6])/(log(x[2]/x[1])*(1+2*x[7]*x[3]/(log(x[2]/x[1])*x[1]^2*x[8])+x[3]/x[5]))
  return(val)
}

p=8

lower<-c(0.05,100,63070,990,63.1,700,1120,9855)
upper<-c(0.15,50000,115600,1110,116,820,1680,12045)

set.seed(1)
N=100000
X1=0.1+0.01618*rnorm(N)
X2=exp(7.71+1.0056*rnorm(N))
X=matrix(0,nrow=N,ncol=p)
X[,1:2]=cbind(X1,X2)
X[,3:p]=matrix(runif(N*6),nrow=N)
X[,3:p]=sweep(X[,3:p],2,upper[3:p]-lower[3:p],"*")
X[,3:p]=sweep(X[,3:p],2,lower[3:p],"+")
true=apply(X,1,f)



n=100
x=X[sample(1:N,n),]
mc=apply(x,1,f)
plot(density(true),lty=3,lwd=3,xlab="y",main="output distribution", cex.main=3,cex.lab=2,cex.axis=2,xlim=c(0,200))
lines(density(mc),col=2,lty=2,lwd=3)

mu=apply(X,2,mean)
sigma=apply(X,2,sd)
Z=sweep(X,2,mu,"-")
Z=sweep(Z,2,sigma,"/")
D=sp(n,p,dist.samp = Z)$sp
D=sweep(D,2,sigma,"*")
D=sweep(D,2,mu,"+")

y=apply(D,1,f)
lines(density(y),col=3,lwd=3)
legend("topright",legend = c("truth","MC","support points"),col=c(1,2,3),lwd=3,lty=c(3,2,1),bty="n",cex=2)


###########################################
#Section 5.2.4 Population Quasi-Monte Carlo
###########################################

#from https://github.com/BillHuang01/PQMC

#Figure 5.9

set.seed(20210816)
setwd("~/R commands/Book/Ch5/PQMC-master")
source("scripts/lib.R")
library(mvtnorm)
library(randtoolbox)

logmixture <- function(x){
  v <- matrix(NA, nrow = 5, ncol = 2)
  v[1,] <- c(-10,-10)
  v[2,] <- c(0,16)
  v[3,] <- c(13,8)
  v[4,] <- c(-9,7)
  v[5,] <- c(14,-14)
  v <- (v + 20) / 40
  sigma <- array(NA, dim = c(5,2,2))
  sigma[1,,] <- matrix(c(2,0.6,0.6,1), nrow = 2) / 40^2
  sigma[2,,] <- matrix(c(2,-0.4,-0.4,2), nrow = 2) / 40^2
  sigma[3,,] <- matrix(c(2,0.8,0.8,2), nrow = 2) / 40^2
  sigma[4,,] <- matrix(c(3,0,0,0.5), nrow = 2) / 40^2
  sigma[5,,] <- matrix(c(2,-0.1,-0.1,2), nrow = 2) / 40^2
  logf <- rep(0,5)
  for (i in 1:5) logf[i] <- dmvnorm(x, mean = v[i,], sigma = sigma[i,,], log = T)
  logf <- logaddexp(logf) - log(5)
  return (logf)
}
x1 <- x2 <- seq(0, 1, length.out = 101)
x.grid <- expand.grid(x1, x2)
density <- matrix(exp(apply(x.grid, 1, logmixture)), 101, 101)

set.seed(1)
Z <- 1 # normalizing constant
p <- 2 # dimensions
K <- 20 # number of propopsals
J <- 10 # number of samples drawn from each proposal
steps <- 4 # number of PMC iterations
# initilization
ini <- sobol(K, p) # initial centers
# initial covariance
sigma <- 0.2 # initial covariance

#PMC (Multinomial)
# setting
sampling <- "random"
resampling <- "multinomial"
sigma.adapt <- FALSE
# call pmc function
layout(matrix(c(1:9), nrow = 3, byrow = T))
pmc.mn <- pmc(logmixture, K, J, steps, ini, 
              sampling = sampling, resampling = resampling, 
              sigma = sigma, sigma.adapt = sigma.adapt,
              visualization = F,output = "weighted samples")


par(mfrow=c(1,2))
contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 10, main = "Population Monte Carlo",cex.main=2, xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2,cex.axis=2)
points(pmc.mn$center.all[(K*steps+1):(K*steps+K),],col="blue",pch=16,cex=2)
points(pmc.mn$samp.all,col=2,pch=4,cex=.5)

#PQMC (ISP + Lookback)
# setting
sampling <- "qmc"
resampling <- "sp"
sigma.adapt <- TRUE # lookback covariance adaptation
# call pmc function
#layout(matrix(c(1:9), nrow = 3, byrow = T))
pqmc.isp <- pmc(logmixture, K, J, steps, ini, 
                sampling = sampling, resampling = resampling, 
                sigma = sigma, sigma.adapt = sigma.adapt,
                visualization = F,output = "weighted samples")


contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 10, main = "Population Quasi-Monte Carlo",cex.main=2, xlab=expression(x[1]), ylab=expression(x[2]),cex.lab=2,cex.axis=2)
points(pqmc.isp$center.all[(K*steps+1):(K*steps+K),],col="blue",pch=16,cex=2)
points(pqmc.isp$samp.all,col=2,pch=4,cex=.5)
par(mfrow=c(1,1))
