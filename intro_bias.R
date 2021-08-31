#####################################################################
# Illustration Bias of Synthetic Control
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# DISCLAIMER: This software is provided "as is" 
# without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
#####################################################################

#####################################################################
# Preliminaries
#####################################################################

rm(list = ls())

library(Synth)
library(xtable)
library(scinference)
library(limSolve)

set.seed(12345)

setwd("/Users/kasparwuthrich/Dropbox/research/SC/SC with Victor and Yinchu/Asymptotics Paper/ttest_sc")

###################################################################
# Functions
###################################################################

source("common_functions.R")

sim.one.sample <- function(T0,T1,J,K.vec,rho.u,mu){
    
  # Preliminaries
  alpha.sig <- 0.1

  # Generate T01
  T01       <- T0+T1
  
  # Weights
  w <- rep(0,J)
  w[1:3] <- 1/3

  # Linear model
  Y0        <- matrix(rnorm(T01*J),T01,J)
  Y0[,1:3]  <- Y0[,1:3] + 2 # Heterogeneity induces bias (i.e., makes selection mistakes induce bias)
  Y1        <- mu + Y0 %*% w + generate.AR.series(T01,rho.u)

  # SC without debiasing
  Y1.pre  <- Y1[1:T0]
  Y0.pre  <- Y0[1:T0,]
  Y1.post <- Y1[(T0+1):T01]
  Y0.post <- Y0[(T0+1):T01,]
  
  w.hat         <- sc(Y1.pre,Y0.pre)
  tau.hat.nodb  <- mean(Y1.post -Y0.post %*% w.hat)
  
  # Results
  tau.vec.db   <- rep(NA,length(K.vec))
  for (k in 1:length(K.vec)){
    tau.vec.db[k] <- scinference(Y1,Y0,T1=T1,T0=T0,inference_method="ttest",K=K.vec[k])$att
  }
  
  return(c(tau.hat.nodb,tau.vec.db))
  
}


#######################################################
# Simulations
#######################################################

### General specification

nreps   <- 20000
K.vec   <- c(2,3)
T0      <- 15
T1      <- 28
J       <- 16

###  Persistent data

rho.u <- 0.8 # note the sampling distribution depends on the LRV which depends on rho.u

result.all.estimate.sc.per <- matrix(NA,nreps,3)
for (r in 1:nreps){
  result.all.estimate.sc.per[r,] <- sim.one.sample(T0,T1,J,K.vec,rho.u,0)
}

result.all.estimate.misspec.sc.per <- matrix(NA,nreps,3)
for (r in 1:nreps){
  result.all.estimate.misspec.sc.per[r,] <- sim.one.sample(T0,T1,J,K.vec,rho.u,2)
}

#######################################################
# Graphics
#######################################################

sc.original.results.per <- result.all.estimate.sc.per[,1]
sc.debiased.results.per <- result.all.estimate.sc.per[,3] # Results for K=3

sc.original.results.misspec.per <- result.all.estimate.misspec.sc.per[,1]
sc.debiased.results.misspec.per <- result.all.estimate.misspec.sc.per[,3] # Results for K=3

xl <- -2.5
xu <- 2.5

yl <- 0
yu <- 1

sc.original.results.per <- sc.original.results.per[xl < sc.original.results.per & sc.original.results.per<xu]
sc.debiased.results.per <- sc.debiased.results.per[xl < sc.debiased.results.per & sc.debiased.results.per<xu]

graphics.off()
pdf("graphics/hist_sc_per.pdf",pointsize=16,width=8.0,height=6.0)
h <- hist(sc.original.results.per, breaks=c(seq(xl,xu,length=40)), col="gray", xlab="", main="Correct specification: SC estimator w/o debiasing",freq=F,xlim=c(xl,xu),ylim=c(yl,yu))
xfit <- seq(min(xl),max(xu),length=1000) 
yfit <- dnorm(xfit,mean=0,sd=sd(sc.original.results.per)) 
lines(xfit, yfit, col="black", lwd=5)
legend("topleft",legend=c("Simulation","Normal curve"), col=c("gray","black"),fill=NA,border=NA, lty=c(1,1),lwd=c(15,5),merge=T,bty="n")
dev.off()

graphics.off()
pdf("graphics/hist_sc_debiased_per.pdf",pointsize=16,width=8.0,height=6.0)
h <- hist(sc.debiased.results.per, breaks=c(seq(xl,xu,length=40)), col="gray", xlab="", main="Correct specification: De-biased SC estimator",freq=F,xlim=c(xl,xu),ylim=c(yl,yu)) 
xfit <- seq(min(xl),max(xu),length=1000) 
yfit <- dnorm(xfit,mean=0,sd=sd(sc.debiased.results.per)) 
lines(xfit, yfit, col="black", lwd=5)
legend("topleft",legend=c("Simulation","Normal curve"), col=c("gray","black"),fill=NA,border=NA, lty=c(1,1),lwd=c(15,5),merge=T,bty="n")
dev.off()

xl <- -3
xu <- 5

yl <- 0
yu <- 1

sc.original.results.misspec.per <- sc.original.results.misspec.per[xl < sc.original.results.misspec.per & sc.original.results.misspec.per<xu]
sc.debiased.results.misspec.per <- sc.debiased.results.misspec.per[xl < sc.debiased.results.misspec.per & sc.debiased.results.misspec.per<xu]

graphics.off()
pdf("graphics/hist_intercept_sc_per.pdf",pointsize=16,width=8.0,height=6.0)
h <- hist(sc.original.results.misspec.per, breaks=c(seq(xl,xu,length=40)), col="gray", xlab="", main="Misspecification: SC estimator w/o debiasing",freq=F,xlim=c(xl,xu),ylim=c(yl,yu))
xfit <- seq(min(xl),max(xu),length=1000) 
yfit <- dnorm(xfit,mean=0,sd=sd(sc.original.results.misspec.per)) 
lines(xfit, yfit, col="black", lwd=5)
legend("topleft",legend=c("Simulation","Normal curve"), col=c("gray","black"),fill=NA,border=NA, lty=c(1,1),lwd=c(10,3),merge=T,bty="n")
dev.off()

graphics.off()
pdf("graphics/hist_intercept_sc_debiased_per.pdf",pointsize=16,width=8.0,height=6.0)
h <- hist(sc.debiased.results.misspec.per, breaks=c(seq(xl,xu,length=40)), col="gray", xlab="", main="Misspecification: De-biased SC estimator",freq=F,xlim=c(xl,xu),ylim=c(yl,yu)) 
xfit <- seq(min(xl),max(xu),length=1000) 
yfit <- dnorm(xfit,mean=0,sd=sd(sc.debiased.results.misspec.per)) 
lines(xfit, yfit, col="black", lwd=5)
legend("topleft",legend=c("Simulation","Normal curve"), col=c("gray","black"),fill=NA,border=NA, lty=c(1,1),lwd=c(10,3),merge=T,bty="n")
dev.off()
