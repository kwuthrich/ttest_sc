#####################################################################
# Illustration size distortions Newey-West
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# DISCLAIMER: This software is provided "as is" 
# without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
#####################################################################

#####################################################################
# Preliminaries
#####################################################################

rm(list = ls())

library(sandwich)

set.seed(12345)

setwd("/Users/kasparwuthrich/Dropbox/research/SC/SC with Victor and Yinchu/Asymptotics Paper")

###################################################################
# Functions
###################################################################

source("Code/common_functions_001.R")

sim <- function(T0,T1,rho,K,alpha.sig){

  T01   <- T0 + T1
  u <- generate.AR.series(T01,rho)
  
  r     <- min(floor(T0/K),T1)
  u.pre  <- u[1:T0]
  u.post <- u[(T0+1):T01]

  tau.mat <- matrix(NA,K,1)
  for (k in 1:K){
    Hk            <- (T0-(r*K))+seq((k-1)*r+1,k*r,1) 
    tau.mat[k,1]  <- mean(u.post) - mean(u.pre[Hk])
  }
  
  tau.hat <- mean(tau.mat)
  se.hat  <- sqrt(1+((K*r)/T1))*sd(tau.mat)/sqrt(K)
  t.hat   <- tau.hat/se.hat # this has a t_{K-1} distribution
  
  rej.t   <- (abs(t.hat)> qt(1-alpha.sig/2,df=K-1))
  
  # estimated lrv
  m       <- lm(u.pre~ 1)
  #lrv    <- NeweyWest(m, prewhite = FALSE)
  lrv     <- NeweyWest(m, prewhite = TRUE)
  #lrv     <- ((1-rho^2)/(1-rho)^2)/T0
  c0      <- T0/T1
  g.c0.K  <- K*(c0<1)+K/c0*(1<=c0&c0<=K)+1*(c0>K)
  std.lrv <- sqrt(lrv*T0*(min(c0,1)+g.c0.K/K))/sqrt(min(T0,T1))
  z.hat   <- tau.hat/std.lrv
  rej.z   <- (abs(z.hat)> qnorm(1-alpha.sig/2))
  
  # true lrv
  lrv         <- ((1-rho^2)/(1-rho)^2)/T0
  c0          <- T0/T1
  g.c0.K      <- K*(c0<1)+K/c0*(1<=c0&c0<=K)+1*(c0>K)
  std.lrv     <- sqrt(lrv*T0*(min(c0,1)+g.c0.K/K))/sqrt(min(T0,T1))
  z.hat       <- tau.hat/std.lrv
  rej.z.true  <- (abs(z.hat)> qnorm(1-alpha.sig/2))
  
  return(c(rej.t,rej.z,rej.z.true))
  
}

###################################################################
# Simulations
###################################################################

nreps <- 20000

T0 <- 15
T1 <- 28

rho.vec <- c(seq(0,0.9,0.1))
alpha.sig <- 0.1

# K=2
K <- 2
results.mat.2 <- matrix(NA,length(rho.vec),3)
for (r in 1:length(rho.vec)){
  temp.mat <- matrix(NA,nreps,3)
  for (rep in 1:nreps){
    temp.mat[rep,] <- sim(T0,T1,rho.vec[r],K,alpha.sig)
  }
  results.mat.2[r,] <- colMeans(temp.mat)
}

# K=3
K <- 3
results.mat.3 <- matrix(NA,length(rho.vec),3)
for (r in 1:length(rho.vec)){
  temp.mat <- matrix(NA,nreps,3)
  for (rep in 1:nreps){
    temp.mat[rep,] <- sim(T0,T1,rho.vec[r],K,alpha.sig)
  }
  results.mat.3[r,] <- colMeans(temp.mat)
}

###################################################################
# Simulations
###################################################################

graphics.off()
pdf("Paper/Graphics/sizedistortion_lrv.pdf",pointsize=16,width=8.0,height=6.0)
plot(1, xlab="AR(1) coefficient", ylab="Size", main="", col="blue", pch=".", xlim = c(0,1), ylim=c(0,0.6))
lines(rho.vec,results.mat.2[,1],col="black",lwd=5, lty=1,type="b",pch=1)
lines(rho.vec,results.mat.3[,1],col="black",lwd=5, lty=1,type="b",pch=2)
lines(rho.vec,results.mat.3[,2],col="black",lwd=5, lty=1,type="b",pch=3)
abline(h=0.1,lty=1,lwd=1,col="gray")
legend("topleft",legend=c("T-test (K=2)","T-test (K=3)","Newey-West standard errors"), seg.len=2, col=c("black","black"),fill=NA,border=NA, lty=c(1,1,1),lwd=c(5,5,5), pch=c(1,2,3), merge=T,bty="n")
dev.off()

graphics.off()
pdf("Paper/Graphics/sizedistortion_lrv_nw_only.pdf",pointsize=16,width=8.0,height=6.0)
plot(1, xlab="AR(1) coefficient", ylab="Size", main="", col="blue", pch=".", xlim = c(0,1), ylim=c(0,0.6))
lines(rho.vec,results.mat.3[,2],col="black",lwd=5, lty=1,type="b",pch=3)
abline(h=0.1,lty=1,lwd=1,col="gray")
legend("topleft",legend=c("Newey-West standard errors"), seg.len=2, col=c("black","black"),fill=NA,border=NA, lty=c(1),lwd=c(5), pch=c(3), merge=T,bty="n")
dev.off()


