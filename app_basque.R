#######################################################
# Application: Abadie and Gardeazabal (2003)
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# DISCLAIMER: This software is provided "as is" 
# without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
#######################################################

#######################################################
# Preliminaries
#######################################################

rm(list = ls())

library(xtable)
library(Synth)
library(scinference)

set.seed(12345)

setwd("/Users/kasparwuthrich/Dropbox/research/SC/SC with Victor and Yinchu/Asymptotics Paper")

#######################################################
# Data
#######################################################

### Data preparation

detach(basque)
data("basque")
attach(basque)
data.long <- as.data.frame(cbind(regionno,year,gdpcap))
data.wide <- reshape(data.long,dir="wide",idvar="regionno",timevar="year")

Y1 <- unlist(data.wide[17,-1])
Y0 <- t(matrix(unlist(data.wide[-c(1,17),-1]),16,43))

J <- dim(Y0)[2]

### Raw data plot

graphics.off()
time <- c(seq(1955,1997,1))
pdf("Paper/Graphics/basque_data_raw.pdf",pointsize=14,width=9.0,height=6.0)
plot(1, ylab="Per Capita GDP", xlab="Time", main="", col="blue", pch=".", xlim = range(time), ylim=c(1,13))
for (j in 1:dim(Y0)[2]) lines(time,Y0[,j],col="black",lwd=0.5,lty=1)
lines(time,Y1,col="black",lwd=5,lty=1)
legend("topleft",legend=c("Basque Country","Controls"), seg.len=2, col=c("black","black"),fill=NA,border=NA, lty=c(1,1),lwd=c(5,0.5), merge=T,bty="n")
dev.off()

### Pre-treatment period: 1955-1969; post-treatment period: 1970-1997

T0 <- 15
T1 <- length(Y1)-T0
T01 <- T0+T1

#######################################################
# Main results
#######################################################

alpha.sig <- 0.1
K.vec  <- c(2:3)

tau.sc.vec <- se.sc.vec <-  t.sc.vec <- lb.sc.vec <- ub.sc.vec <-  matrix(NA,length(K.vec),1)

for (k in 1:length(K.vec)){

  r.sc.sim          <- scinference(Y1,Y0,T1=T1,T0=T0,inference_method="ttest",K=K.vec[k])
  tau.sc.vec[k,1]   <- r.sc.sim$att
  se.sc.vec[k,1]    <- r.sc.sim$se
  t.sc.vec[k,1]     <- r.sc.sim$att/r.sc.sim$se
  lb.sc.vec[k,1]    <- r.sc.sim$lb  
  ub.sc.vec[k,1]    <- r.sc.sim$ub  
  
}

results.sc  <- cbind(tau.sc.vec,se.sc.vec,t.sc.vec,lb.sc.vec,ub.sc.vec)

round(results.sc,digits=3)
round(qt(1-alpha.sig/2,df=1),digits=2)
round(qt(1-alpha.sig/2,df=2),digits=2)



