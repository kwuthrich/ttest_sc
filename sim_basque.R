#####################################################################
# Empirical Monte Carlo based on Abadie and Gardeazabal (2003)
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

library(rgl)
set.seed(12345)

###################################################################
# Functions
###################################################################

source("common_functions.R")

sim.one.sample <- function(DGP,T0,T1,J,K,Lambda,rho.u,var.u,var.factors,rho.vec,var.epsl.vec,w0.sc){

  # Preliminaries
  alpha.sig <- 0.1

  # DGPs

  if (DGP==1){ # SC (estimated)
    w       <- w0.sc
    nonstat <- matrix(0,(T0+T1),J)
  }
  if (DGP==2){ # DID 
    w       <- matrix(1,J,1)/J
    nonstat <- matrix(0,(T0+T1),J)
  }
  if (DGP==3){ # Full misspecification
    w       <- -matrix(seq(1,J,1),J,1)/J
    nonstat <- matrix(0,(T0+T1),J)
  }
  if (DGP==4){ # Common deterministic trend
    w       <- w0.sc
    theta   <- c(1:(T0+T1))
    nonstat <- matrix(rep(theta,J),(T0+T1),J)
  }
  if (DGP==5){ # Common random walk
    w       <- w0.sc
    theta   <- cumsum(rnorm((T0+T1)))
    nonstat <- matrix(rep(theta,J),(T0+T1),J)
  }
  if (DGP==6){ # Common deterministic trend + sparse deviation
    w           <- w0.sc
    theta       <- c(1:(T0+T1))
    nonstat     <- matrix(rep(theta,J),(T0+T1),J)
    deviation   <- c(1:(T0+T1))
    nonstat[,1] <- nonstat[,1] + deviation # note that the first element of w0.sc is zero.
  }
  if (DGP==7){ # Common random walk + sparse deviation
    w           <- w0.sc
    theta       <- cumsum(rnorm((T0+T1)))
    nonstat     <- matrix(rep(theta,J),(T0+T1),J)
    deviation   <- cumsum(rnorm((T0+T1)))
    nonstat[,1] <- nonstat[,1] + deviation # note that the first element of w0.sc is zero.
  }
  if (DGP==8){ # Heterogenous deterministic trends
    w       <- w0.sc
    nonstat <- matrix(NA,(T0+T1),J)
    for (j in 1:J){
      nonstat[,j] <- j+j*c(1:(T0+T1))
    }
  }
  if (DGP==9){ # Random walks with heterogenous drifts
    w       <- w0.sc
    nonstat <- matrix(NA,(T0+T1),J)
    for (j in 1:J){
      theta     <- rep(NA,(T0+T1))
      wn        <- rnorm((T0+T1))
      theta[1]  <- wn[1]
      for (t in 2:(T0+T1)){
        theta[t] <- j+theta[(t-1)] + wn[t]
      }
      nonstat[,j] <- theta
    }
  }

  # Generate T01
  T01 <- (T0+T1)

  # Generate errors as AR(1)
  u.sim <- generate.AR.series.var.input(T01,rho.u,var.u)

  # Generate iid factors
  Factors.sim <- matrix(NA,T01,num.factor)
  for (i in 1:num.factor){
    Factors.sim[,i] <- rnorm(T01)*sqrt(var.factors[i])
  }

  # Generate control errors as AR(1)
  epsl.sim <- matrix(NA,T01,J)
  for (i in 1:J){
    epsl.sim[,i] <- generate.AR.series.var.input(T01,rho.vec[i],var.epsl.vec[i])
  }

  # Linear model (under the null of no effect)
  Y0 <- nonstat + Factors.sim %*% t(Lambda) + epsl.sim
  Y1 <- Y0 %*% w + u.sim

  # Results
  r.sc.sim    <- scinference(Y1,Y0,T1=T1,T0=T0,inference_method="ttest",estimation_method="sc",K=K)
  tau.hat.sc  <- r.sc.sim$att
  se.hat.sc   <- r.sc.sim$se
  cov.sc.vec  <- (abs(tau.hat.sc/se.hat.sc)<=qt(1-alpha.sig/2,df=K-1))
  leng.sc.vec <- 2*qt(1-alpha.sig/2,df=K-1)*se.hat.sc

  r.did.sim     <- scinference(Y1,Y0,T1=T1,T0=T0,inference_method="ttest",estimation_method="did",K=K)
  tau.hat.did   <- r.did.sim$att
  se.hat.did    <- r.did.sim$se
  cov.did.vec   <- (abs(tau.hat.did/se.hat.did)<=qt(1-alpha.sig/2,df=K-1))
  leng.did.vec  <- 2*qt(1-alpha.sig/2,df=K-1)*se.hat.did

  w.true          <- w
  r.oracle.sim    <- oracle.cf(Y1,Y0,T0,T1,w.true,K)
  tau.hat.oracle  <- r.oracle.sim$tau.hat
  se.hat.oracle   <- r.oracle.sim$se.hat
  cov.oracle.vec  <- (abs(tau.hat.oracle/se.hat.oracle)<=qt(1-alpha.sig/2,df=K-1))
  leng.oracle.vec <- 2*qt(1-alpha.sig/2,df=K-1)*se.hat.oracle

  length.all    <- c(leng.sc.vec,leng.did.vec,leng.oracle.vec)
  coverage.all  <- c(cov.sc.vec,cov.did.vec,cov.oracle.vec)
  result        <- c(coverage.all,length.all)

  return(result)

}


###################################################################
# Calibration to application
###################################################################
data("basque")
attach(basque)
data.long <- as.data.frame(cbind(regionno,year,gdpcap))
data.wide <- reshape(data.long,dir="wide",idvar="regionno",timevar="year")
detach(basque)

# Raw data

Y1 <- unlist(data.wide[17,-1])
Y0 <- t(matrix(unlist(data.wide[-c(1,17),-1]),16,43))

J <- dim(Y0)[2]

### Pre-treatment period: 1955-1969; post-treatment period: 1970-1997

T0 <- 15
T1 <- length(Y1)-T0
T01 <- T0+T1

### De-trend data

bX <- rowMeans(Y0)

Y1dt <- Y1-bX
Y0dt <- apply(Y0,2,function(x) x - bX)

### Fit factor model

w0.sc         <- as.matrix(scinference:::sc(Y1dt[1:T0],Y0dt[1:T0,],lsei_type = 1)$w.hat)
u.hat         <- Y1dt[1:T0] - Y0dt[1:T0,] %*% w0.sc
ar.obj.u.hat  <- ar(u.hat, order.max=1)
rho.u         <- ar.obj.u.hat$ar
var.u         <- ar.obj.u.hat$var.pred # this is the variance of the AR error term

result <- svd(Y0dt)
max(abs(result$u %*% diag(result$d) %*% t(result$v)-Y0dt))

num.factor            <- 4
dd                    <- result$d
dd[(num.factor+1):J]  <- 0
Y0.PCA.est            <- result$u %*% diag(dd) %*% t(result$v)
Y0.resid              <- Y0dt-Y0.PCA.est

Factors <- result$u[,1:num.factor]
Lambda  <- result$v %*% diag(dd)
Lambda  <- Lambda[,1:num.factor]

var.factors <- matrix(NA,num.factor,1)
for (j in 1:num.factor){
  var.factors[j]  <- var(Factors[,j])
}

rho.vec       <- matrix(NA,J,1)
var.epsl.vec  <- matrix(NA,J,1)
for (j in 1:J){
  ar.obj          <- ar(Y0.resid[,j], order.max=1)
  rho.vec[j,1]    <- ar.obj$ar
  var.epsl.vec[j] <- ar.obj$var.pred
}

round(w0.sc,digits=2)
round(rho.u,digits=2)
round(c(min(rho.vec),max(rho.vec),median(rho.vec)),digits=2)

#######################################################
# Simulations
#######################################################

### Simulations for different DGPs

nreps <- 10000

for (DGP in 1:9){

  results_K2_15 <- results_K3_15 <- results_K2_150 <- results_K3_150 <- matrix(NA,nreps,6)

  for (r in 1:nreps){
    results_K2_15[r,]   <- sim.one.sample(DGP,15,T1,J,2,Lambda,rho.u,var.u,var.factors,rho.vec,var.epsl.vec,w0.sc)
    results_K3_15[r,]   <- sim.one.sample(DGP,15,T1,J,3,Lambda,rho.u,var.u,var.factors,rho.vec,var.epsl.vec,w0.sc)
    results_K2_150[r,]  <- sim.one.sample(DGP,150,T1,J,2,Lambda,rho.u,var.u,var.factors,rho.vec,var.epsl.vec,w0.sc)
    results_K3_150[r,]  <- sim.one.sample(DGP,150,T1,J,3,Lambda,rho.u,var.u,var.factors,rho.vec,var.epsl.vec,w0.sc)
  }

  overall_K2_15   <- colMeans(results_K2_15)
  overall_K3_15   <- colMeans(results_K3_15)
  overall_K2_150  <- colMeans(results_K2_150)
  overall_K3_150  <- colMeans(results_K3_150)

  overall_K2 <- c(overall_K2_15,overall_K2_150)
  overall_K3 <- c(overall_K3_15,overall_K3_150)

  print(DGP)

  print(round(c(overall_K2_15[4]/overall_K2_15[5],overall_K3_15[4]/overall_K3_15[5]),digits=2))
  print(round(c(overall_K2_150[4]/overall_K2_150[5],overall_K3_150[4]/overall_K3_150[5]),digits=2))

  print(xtable(cbind(c(2,3),rbind(overall_K2,overall_K3)),digits=c(0,0,rep(2,12))),include.rownames = F)

}


### Illustrating trade-off when choosing K

Ks  <- 2:5
T0s <- c(20,40,60,80,100)

overall_cov_Ks <- overall_leng_Ks<- matrix(NA,length(Ks),length(T0s))

for (k in 1:length(Ks)){
  
  for (t in 1:length(T0s)){
    
    res_temp <- matrix(NA,nreps,6)
    
    for (r in 1:nreps){
      res_temp[r,]   <- sim.one.sample(1,T0s[t],T1,J,Ks[k],Lambda,rho.u,var.u,var.factors,rho.vec,var.epsl.vec,w0.sc)
    }
    
    overall_cov_Ks[k,t] <- colMeans(res_temp)[1]
    overall_leng_Ks[k,t] <- colMeans(res_temp)[4]
    
  }
  
}

#######################################################
# Figures trade-off
#######################################################

# The coloring is based on the example here: https://rdrr.io/r/graphics/persp.html

cov_axis <- c(0.6,1)
leng_axis <- c(0,0.6)

nrowz <- nrow(overall_cov_Ks)
ncolz <- ncol(overall_cov_Ks)
jet.colors <- colorRampPalette( c("orange", "blue") )
ncol <- 100
color <- jet.colors(ncol)

midpoints_cov <- overall_cov_Ks[-1, -1] + overall_cov_Ks[-1, -ncolz] + overall_cov_Ks[-nrowz, -1] + overall_cov_Ks[-nrowz, -ncolz]
pdf("graphics/3d_cov.pdf",pointsize=18,width=10.0,height=8.0)
pmat <- persp(x=Ks,y=T0s,z=overall_cov_Ks,col = color[cut(midpoints_cov,ncol)], xlim=range(Ks),ylim=range(T0s),
      zlim=range(cov_axis),xlab="K",ylab="Number of pre-treatment periods",zlab="Coverage",ticktype="detailed",expand=0.8,theta = 50,phi=20,d=10,nticks=4,main="Coverage 90% confidence intervals")
dev.off()

jet.colors <- colorRampPalette( c("blue", "orange") )
ncol <- 100
color <- jet.colors(ncol)
midpoints_leng <- overall_leng_Ks[-1, -1] + overall_leng_Ks[-1, -ncolz] + overall_leng_Ks[-nrowz, -1] + overall_leng_Ks[-nrowz, -ncolz]
pdf("graphics/3d_leng.pdf",pointsize=18,width=10.0,height=8.0)
pmat <- persp(x=Ks,y=T0s,z=overall_leng_Ks,col = color[cut(midpoints_leng, ncol)],xlim=range(Ks),ylim=range(T0s),
      zlim=range(leng_axis),xlab="K",ylab="Number of pre-treatment periods",zlab="Average length",ticktype="detailed",expand=0.8,theta = 50,phi=20,d=10,nticks=4,main="Average length 90% confidence intervals")
dev.off()

