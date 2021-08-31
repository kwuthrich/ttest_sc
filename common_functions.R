#####################################################################
# Functions
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# DISCLAIMER: This software is provided "as is" 
# without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
#####################################################################

# AR(1) with given error variance
generate.AR.series.var.input <- function(T01,rho,var.error){
  u     <- matrix(NA,T01,1)
  epsl  <- matrix(rnorm(T01),T01,1)*sqrt(var.error)
  startvalue <- sqrt(var.error/(1-rho^2))
  u[1,] <- rho*startvalue +epsl[1,]
  for (t in 2:T01){
    u[t,] <- rho*u[(t-1),]+epsl[t,]
  }
  return(u)
}

# AR(1)
generate.AR.series <- function(T01,rho){
  u           <- matrix(NA,T01,1)
  epsl        <- matrix(rnorm(T01),T01,1)*sqrt(1-rho^2)
  startvalue  <- rnorm(1)
  u[1,]       <- rho*startvalue + epsl[1,] 
  for (t in 2:T01){
    u[t,] <- rho*u[(t-1),]+epsl[t,]
  }
  return(u)
}

# SC estimator of w 
sc <- function(Y1,Y0){
  E <- matrix(1,1,ncol(Y0))
  F <- 1
  G <- diag(ncol(Y0))
  H <- matrix(0,ncol(Y0),1)
  w.hat <- cbind(lsei(A=Y0, B=Y1, E=E, F=F, G=G, H=H, type=2)$X)
  return(w.hat)
}

# Oracle approach based on known w
oracle.cf <- function(y1,yJ,T0,T1,w.true,K){
  
  T01   <- T0 + T1
  r     <- min(floor(T0/K),T1)
  y1.pre  <- y1[1:T0]
  yJ.pre  <- yJ[1:T0,]
  y1.post <- y1[(T0+1):T01]
  yJ.post <- yJ[(T0+1):T01,]
  
  tau.mat <- matrix(NA,K,1)
  for (k in 1:K){
    Hk            <- (T0-(r*K))+seq((k-1)*r+1,k*r,1) 
    tau.mat[k,1]  <- mean(y1.post-yJ.post%*%w.true) - mean(y1.pre[Hk]-yJ.pre[Hk,]%*%w.true)
  }
  
  tau.hat <- mean(tau.mat)
  se.hat  <- sqrt(1+((K*r)/T1))*sd(tau.mat)/sqrt(K)
  t.hat   <- tau.hat/se.hat # this has a t_{K-1} distribution
  
  return(list(t.hat=t.hat,tau.hat=tau.hat,se.hat=se.hat))
}


