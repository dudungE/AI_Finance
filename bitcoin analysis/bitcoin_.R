########################################### 
####R Mid term  
####경제학부    
####20176386    
####정재현     
###########################################


#데이터불러오기
data1 <- read.csv("C:/Users/JJHR/Desktop/R midterm/data_BitcoinHedge.csv", head=T)
data <- data1[,(-9:-13)]
summary(data)


#데이터정제~시계열, 로그수익률
date <- as.POSIXct(data$date)


lag(data$btc)
install.packages("plyr")
library(plyr)
f_lnR <- function(x) (diff(log(x), lag = 1)*100)
lnR <- plyr::colwise(f_lnR)(data[, -1])
lnR





#######################기초통계량 분석  
#####graph그리기
#btc
y1 <- rep(NA,1166)
y1[1] <- log(data$btc)[1]
for (i in 2:1166){
  y1[i] <- (log(data$btc)[i] - log(data$btc)[i-1])*100
}

par(mfrow=c(2,1), cex=0.5, mex=0.5)
plot(date,data$btc,type="l",lty=1,col="blue",lwd=2,ylab="btc",xlab="Time")
plot(date[-1],y1[-1],type="l",lty=1,col="blue",lwd=2,ylab="btc",xlab="Time")
plot(date[-1],lnR$btc,type="l",lty=1,col="blue",lwd=2,ylab="btc",xlab="Time")

#ltc
y2 <- rep(NA,1166)
y2[1] <- log(data$ltc)[1]
for (i in 2:1166){
  y2[i] <- (log(data$ltc)[i] - log(data$ltc)[i-1])*100
}
plot(date,data$ltc,type="l",lty=1,col="blue",lwd=2,ylab="ltc",xlab="Time")
plot(date[-1],lnR$ltc,type="l",lty=1,col="blue",lwd=2,ylab="ltc",xlab="Time")

#xrp  #이것부터 로그수익률 그릴 때 간단하게 lnR데이터만이용
plot(date,data$xrp,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")
plot(date[-1],lnR$xrp,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")

#sp
plot(date,data$sp,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")
plot(date[-1],lnR$sp,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")


#gold
plot(date,data$gold,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")
plot(date[-1],lnR$gold,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")


#bond
plot(date,data$bond,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")
plot(date[-1],lnR$bond,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")

#vix
plot(date,data$vix,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")
plot(date[-1],lnR$vix,type="l",lty=1,col="blue",lwd=2,ylab="xrp",xlab="Time")


#기초통계량
summary(lnR)
install.packages("pastecs")
library(pastecs)
stat.desc(lnR)







#####
#####모형의 추정 및 분석
#####

##회귀분석
#Rc,t = α + β1Rs,t + β2Rs,t(τ) + δ1Rb,t + δ2Rb,t(τ) + γ1Rg,t + γ2Rg,t(τ) + ut

### i) τ = 0.1
#1.btc
spt <- ifelse(lnR$sp <= quantile(lnR$sp,0.1), lnR$sp, 0)
bondt <- ifelse(lnR$bond <= quantile(lnR$bond,0.1), lnR$bond, 0)
goldt <- ifelse(lnR$gold <= quantile(lnR$gold,0.1), lnR$gold, 0)
data_Rt <- cbind(spt,goldt,bondt)
data_lnRt <- cbind(lnR,data_Rt)
attach(data_lnRt)


res1 <- lm(btc~sp+spt+bond+bondt+gold+goldt)
summary(res1)
serial_corr <- lm(res1$resid[-1]~res1$resid[-1165])
summary(serial_corr)

library(sandwich)
sqrt(diag(vcovHAC(res1,type="HAC")))




prais.winsten.lm <- function(mod){
  X <- model.matrix(mod)
  y <- model.response(model.frame(mod))
  e <- residuals(mod)
  n <- length(e)
  names <- colnames(X)
  rho <- sum(e[1:(n-1)]*e[2:n])/sum(e^2)
  y <- c(y[1] * (1 - rho^2)^0.5, y[2:n] - rho * y[1:(n-1)])
  X <- rbind(X[1,] * (1 - rho^2)^0.5, X[2:n,] - rho * X[1:(n-1),])
  mod <- lm(y ~ X - 1)
  result <- list()
  result$coefficients <- coef(mod)
  names(result$coefficients) <- names
  summary <- summary(mod, corr = F)
  result$cov <- (summary$sigma^2) * summary$cov.unscaled
  dimnames(result$cov) <- list(names, names)
  result$sigma <- summary$sigma
  result$rho <- rho
  class(result) <- 'prais.winsten'
  result
}

prais.winsten.lm(res1)
sqrt(diag(prais.winsten.lm(res1)$cov))



#2.ltc
res2 <- lm(ltc~sp+spt+bond+bondt+gold+goldt)
summary(res2)
serial_corr <- lm(res2$resid[-1]~res2$resid[-1165])
summary(serial_corr)

#3.xrp
res3 <- lm(xrp~sp+spt+bond+bondt+gold+goldt)
summary(res3)
serial_corr <- lm(res3$resid[-1]~res3$resid[-1165])
summary(serial_corr)



### ii) τ = 0.05
#1.btc
spt <- ifelse(lnR$sp <= quantile(lnR$sp,0.05), lnR$sp, 0)
bondt <- ifelse(lnR$bond <= quantile(lnR$bond,0.05), lnR$bond, 0)
goldt <- ifelse(lnR$gold <= quantile(lnR$gold,0.05), lnR$gold, 0)
data_Rt <- cbind(spt,goldt,bondt)
data_lnRt <- cbind(lnR,data_Rt)
attach(data_lnRt)

res1 <- lm(btc~sp+spt+bond+bondt+gold+goldt)
summary(res1)
serial_corr <- lm(res1$resid[-1]~res1$resid[-1165])
summary(serial_corr)

#2.ltc
res2 <- lm(ltc~sp+spt+bond+bondt+gold+goldt)
summary(res2)
serial_corr <- lm(res2$resid[-1]~res2$resid[-1165])
summary(serial_corr)


#3.xrp
res3 <- lm(xrp~sp+spt+bond+bondt+gold+goldt)
summary(res3)
serial_corr <- lm(res3$resid[-1]~res3$resid[-1165])
summary(serial_corr)



### iii) τ = 0.01
#1.btc
spt <- ifelse(lnR$sp <= quantile(lnR$sp,0.01), lnR$sp, 0)
bondt <- ifelse(lnR$bond <= quantile(lnR$bond,0.01), lnR$bond, 0)
goldt <- ifelse(lnR$gold <= quantile(lnR$gold,0.01), lnR$gold, 0)
data_Rt <- cbind(spt,goldt,bondt)
data_lnRt <- cbind(lnR,data_Rt)
attach(data_lnRt)

res1 <- lm(btc~sp+spt+bond+bondt+gold+goldt)
summary(res1)
serial_corr <- lm(res1$resid[-1]~res1$resid[-1165])
summary(serial_corr)

#2.ltc
res2 <- lm(ltc~sp+spt+bond+bondt+gold+goldt)
summary(res2)
serial_corr <- lm(res2$resid[-1]~res2$resid[-1165])
summary(serial_corr)


#3.xrp
res3 <- lm(xrp~sp+spt+bond+bondt+gold+goldt)
summary(res3)
serial_corr <- lm(res3$resid[-1]~res3$resid[-1165])
summary(serial_corr)
















#############################분산투자

############1.첫 번째 집합(stock,bond,gold)
####i) Mean-Variance w/o shortsale

#표본으로부터 portfolio weights 구하기
install.packages("quadprog")
require(quadprog)

rt <- lnR[4:6]
dim(rt)

mu.hat <- colMeans(rt)
Sigma <- cov(rt)
mu.g <- mean(mu.hat)

Dmat <- Sigma
dvec <- rep(0,3)
Amat <- t(rbind( mu.hat, c(1,1,1), diag(rep(1,3)) ))
bvec <- c(mu.g, 1, rep(0,3))
meq <- 2
res.mv <- solve.QP(Dmat, dvec, Amat, bvec, meq)
res.mv$solution

#표본외 수익률 구하기(표본내 마지막 한달기준)
rt0 <- rt[1135:1165,]
rt0

rp <- t(res.mv$solution)%*%t(as.matrix(rt0))
mu <- mean(rp)
sigma <- sd(rp)
SR <- mu/sigma
SR

#diversification measure
#허핀달지수
HHI <- 1-sum(res.mv$solution^2)
HHI

#weight entrophy
WE <- -sum(res.mv$solution*log(res.mv$solution), na.rm = T)
WE

#diversification ratio
DR <- sum(res.mv$solution*diag(Sigma))/sigma
DR



#####ii) Minimum-Variance w/o shortsale

#표본으로부터 portfolio weights 구하기
mu.hat <- colMeans(rt)
Sigma <- cov(rt)
mu.g <- mean(mu.hat)

Dmat <- Sigma
dvec <- rep(0,3)
Amat <- t(rbind( c(1,1,1), diag(rep(1,3)) ))
bvec <- c(1, rep(0,3))
meq <- 1
res.mv <- solve.QP(Dmat, dvec, Amat, bvec, meq)
res.mv


#표본외 수익률 구하기
rp <- t(res.mv$solution)%*%t(as.matrix(rt0))
mu <- mean(rp)
sigma <- sd(rp)
SR <- mu/sigma
SR

##diversification measure
#허핀달지수
HHI <- 1-sum(res.mv$solution^2)
HHI

#weight entrophy
WE <- -sum(res.mv$solution*log(res.mv$solution), na.rm = T)
WE

#diversification ratio
DR <- sum(res.mv$solution*diag(Sigma))/sigma
DR


####iii) Mean-Variance w/o shortsale - resampled

install.packages("mvtnorm")
library(mvtnorm)

B <- 500
rs.pi <- matrix(NA, ncol=3, nrow=B)

for (i in 1:B){
  r.rt <-rmvnorm(500, mu.hat, Sigma)
  r.mu.hat <- colMeans(r.rt)
  r.Sigma <- cov(r.rt)
  r.mu.g <- mean(r.mu.hat)
  dvec <- rep(0,3)
  r.Dmat <- r.Sigma
  r.bvec <- c(r.mu.g, 1, rep(0,3))
  r.Amat <- t(rbind(r.mu.hat, c(1,1,1), diag(rep(1,3))))
  meq <- 2
  rs.mv <- solve.QP(r.Dmat, dvec, r.Amat, r.bvec, meq)
  rs.pi[i,] <- rs.mv$solution
}

resampled.mv <- colMeans(rs.pi)
resampled.mv
sum(resampled.mv)
rp <- t(resampled.mv)%*%t(as.matrix(rt0))
mu <- mean(rp)
sigma <- sd(rp)
SR <- mu/sigma
SR

##diversification measure
#허핀달지수
HHI <- 1-sum(resampled.mv^2)
HHI

#weight entrophy
WE <- -sum(resampled.mv*log(resampled.mv), na.rm = T)
WE

#diversification ratio
DR <- sum(resampled.mv*diag(Sigma))/sigma
DR



####iv) Minimum-Variance w/o shortsale - resampled
B <- 500
rs.pi <- matrix(NA, ncol=3, nrow=B)

for (i in 1:B){
  r.rt <-rmvnorm(500, mu.hat, Sigma)
  r.mu.hat <- colMeans(r.rt)
  r.Sigma <- cov(r.rt)
  r.mu.g <- mean(r.mu.hat)
  dvec <- rep(0,3)
  r.Dmat <- r.Sigma
  r.bvec <- c(1, rep(0,3))
  r.Amat <- t(rbind(c(1,1,1), diag(rep(1,3))))
  meq <- 1
  rs.mv <- solve.QP(r.Dmat, dvec, r.Amat, r.bvec, meq)
  rs.pi[i,] <- rs.mv$solution
}

resampled.mv <- colMeans(rs.pi)
resampled.mv
sum(resampled.mv)

rp <- t(resampled.mv)%*%t(as.matrix(rt0))
rp
mu <- mean(rp)
sigma <- sd(rp)
SR <- mu/sigma
SR

##diversification measure
#허핀달지수
HHI <- 1-sum(resampled.mv^2)
HHI

#weight entrophy
WE <- -sum(resampled.mv*log(resampled.mv), na.rm = T)
WE

#diversification ratio
DR <- sum(resampled.mv*diag(Sigma))/sigma
DR



########2.두 번째 집합(stock,bond,gold+btc,ltc,xrp)

####i) Mean-Variance w/o shortsale
#표본으로부터 portfolio weights 구하기
rt <- lnR[1:6]
dim(rt)

mu.hat <- colMeans(rt)
Sigma <- cov(rt)
mu.g <- mean(mu.hat)

Dmat <- Sigma
dvec <- rep(0,6)
Amat <- t(rbind( mu.hat, c(1,1,1,1,1,1), diag(rep(1,6)) ))
bvec <- c(mu.g, 1, rep(0,6))
meq <- 2
res.mv <- solve.QP(Dmat, dvec, Amat, bvec, meq)
res.mv$solution

#표본외 수익률 구하기(표본내 마지막 한달기준)
rt0 <- rt[1135:1165,]
rt0

rp <- t(res.mv$solution)%*%t(as.matrix(rt0))
mu <- mean(rp)
sigma <- sd(rp)
SR <- mu/sigma
SR

##diversification measure
#허핀달지수
HHI <- 1-sum(res.mv$solution^2)
HHI

#weight entrophy
WE <- -sum(res.mv$solution*log(res.mv$solution), na.rm = T)
WE

#diversification ratio
DR <- sum(res.mv$solution*diag(Sigma))/sigma
DR



####ii) Minimum-Variance w/o shortsale

#표본으로부터 portfolio weights 구하기
mu.hat <- colMeans(rt)
Sigma <- cov(rt)
mu.g <- mean(mu.hat)

Dmat <- Sigma
dvec <- rep(0,6)
Amat <- t(rbind( c(1,1,1,1,1,1), diag(rep(1,6)) ))
bvec <- c(1, rep(0,6))
meq <- 1
res.mv <- solve.QP(Dmat, dvec, Amat, bvec, meq)
res.mv


#표본외 수익률 구하기
rp <- t(res.mv$solution)%*%t(as.matrix(rt0))
mu <- mean(rp)
sigma <- sd(rp)
SR <- mu/sigma
SR

##diversification measure
#허핀달지수
HHI <- 1-sum(res.mv$solution^2)
HHI

#weight entrophy
WE <- -sum(res.mv$solution*log(res.mv$solution), na.rm = T)
WE

#diversification ratio
DR <- sum(res.mv$solution*diag(Sigma))/sigma
DR


####iii) Mean-Variance w/o shortsale - resampled

install.packages("mvtnorm")
library(mvtnorm)

B <- 500
rs.pi <- matrix(NA, ncol=6, nrow=B)

for (i in 1:B){
  r.rt <-rmvnorm(500, mu.hat, Sigma)
  r.mu.hat <- colMeans(r.rt)
  r.Sigma <- cov(r.rt)
  r.mu.g <- mean(r.mu.hat)
  dvec <- rep(0,6)
  r.Dmat <- r.Sigma
  r.bvec <- c(r.mu.g, 1, rep(0,6))
  r.Amat <- t(rbind(r.mu.hat, c(1,1,1,1,1,1), diag(rep(1,6))))
  meq <- 2
  rs.mv <- solve.QP(r.Dmat, dvec, r.Amat, r.bvec, meq)
  rs.pi[i,] <- rs.mv$solution
}

resampled.mv <- colMeans(rs.pi)
resampled.mv
sum(resampled.mv)

rp <- t(resampled.mv)%*%t(as.matrix(rt0))
mu <- mean(rp)
sigma <- sd(rp)
SR <- mu/sigma
SR

##diversification measure
#허핀달지수
HHI <- 1-sum(resampled.mv^2)
HHI

#weight entrophy
WE <- -sum(resampled.mv*log(resampled.mv), na.rm = T)
WE

#diversification ratio
DR <- sum(resampled.mv*diag(Sigma))/sigma
DR



####iv) Minimum-Variance w/o shortsale - resampled
B <- 500
rs.pi <- matrix(NA, ncol=6, nrow=B)

for (i in 1:B){
  r.rt <-rmvnorm(500, mu.hat, Sigma)
  r.mu.hat <- colMeans(r.rt)
  r.Sigma <- cov(r.rt)
  r.mu.g <- mean(r.mu.hat)
  dvec <- rep(0,6)
  r.Dmat <- r.Sigma
  r.bvec <- c(1, rep(0,6))
  r.Amat <- t(rbind(c(1,1,1,1,1,1), diag(rep(1,6))))
  meq <- 1
  rs.mv <- solve.QP(r.Dmat, dvec, r.Amat, r.bvec, meq)
  rs.pi[i,] <- rs.mv$solution
}

resampled.mv <- colMeans(rs.pi)
resampled.mv
sum(resampled.mv)

rp <- t(resampled.mv)%*%t(as.matrix(rt0))
rp
mu <- mean(rp)
sigma <- sd(rp)
SR <- mu/sigma
SR

##diversification measure
#허핀달지수
HHI <- 1-sum(resampled.mv^2)
HHI

#weight entrophy
WE <- -sum(resampled.mv*log(resampled.mv), na.rm = T)
WE

#diversification ratio
DR <- sum(resampled.mv*diag(Sigma))/sigma
DR

