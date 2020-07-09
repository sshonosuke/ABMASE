set.seed(2)
source("ABMASE-function.R")
source("GREG-function.R")


# settings
p <- 30    # maximum number of covariates
N <- 10000    # population size
n <- 300    # sample size

# covariates
rho <- 0.2
C <- matrix(NA,p,p)
for(j in 1:p){ 
  for(k in 1:p){ C[k,j] <- rho^(abs(j-k)) }
}
D <- diag( rep(sqrt(2), p) )
X <- mvrnorm(N, rep(1,p), D%*%C%*%D)
mX <- apply(X, 2, mean)    # population mean vector


# true coeffieicnts and varaince 
Beta <- rep(0, p)
Beta[c(1, 4, 7, 10)] <- c(1, -0.5, 1, -0.5)
sig <- 2

# response
Y <- as.vector( X%*%Beta + sig*rnorm(N) )

# true mean
Mu <- mean(Y)   




###----------------------------------------------------------------###
###          Scenario 1 (simple random sampling)                   ###
###----------------------------------------------------------------###
Pi <- rep(n/N, N)
Sel <- sort( sample(1:N, n) )   # simple random sampling
Ind <- rep(0, N)
Ind[Sel] <- 1
sY <- Y[Ind==1]
sX <- X[Ind==1,]
sPi <- Pi[Ind==1]


# HT estimator
Delta <- matrix(NA,n,n)
for(i in 1:n){
  Delta[i,i] <- 1-n/N
  pp2 <- n*(n-1)/(N*(N-1))
  Delta[i,-i] <- 1-(n/N)^2/pp2
}
HT <- sum(sY/sPi)/N
HTV <- as.vector( t(sY/sPi)%*%Delta%*%(sY/sPi)/N^2 )

# GREG estimator
greg1 <- GREG(sY, sX, sPi, mX, shrink="NO", samp="SRS")$Greg
greg2 <- GREG(sY, sX, sPi, mX, shrink="Lasso", samp="SRS")$Greg
greg3 <- GREG(sY, sX, sPi, mX, shrink="Ridge", samp="SRS")$Greg
greg4 <- GREG.forward(sY, sX, sPi, mX, samp="SRS")$Greg
greg5 <- GREG.mixed(sY, sX, sPi, mX, samp="SRS")$Greg

# Proposed methods
fit1 <- ABMA.LM(sY, sX, sPi, mX, N, samp="SRS", prior="default", mc=3000, burn=1000)
fit2 <- ABMA.LM(sY, sX, sPi, mX, N, samp="SRS", prior="Laplace", mc=3000, burn=1000)
fit3 <- ABMA.LM(sY, sX, sPi, mX, N, samp="SRS", prior="HS", mc=3000, burn=1000)


##  point estimates  ##
HT
greg1[1]
greg2[1]
greg3[1]
greg4[1]
greg5[1]
mean(fit1$Mu)
mean(fit2$Mu)
mean(fit3$Mu)
Mu  # true value 


##  95% confidence (credible) intervals  ##
zz <- qnorm(0.975)
HT+c(-1,1)*zz*sqrt(HTV)
greg1[1]+c(-1,1)*zz*sqrt(greg1[2])
greg2[1]+c(-1,1)*zz*sqrt(greg2[2])
greg3[1]+c(-1,1)*zz*sqrt(greg3[2])
greg4[1]+c(-1,1)*zz*sqrt(greg4[2])
greg5[1]+c(-1,1)*zz*sqrt(greg5[2])
quantile(fit1$Mu, prob=c(0.025, 0.975))
quantile(fit2$Mu, prob=c(0.025, 0.975))
quantile(fit3$Mu, prob=c(0.025, 0.975))





###----------------------------------------------------------------###
###    Scenario 2 (probability-proportional-to-size sampling)      ###
###----------------------------------------------------------------###
zz <- log( 1+abs(Y+rexp(N,2)) )   # size variable
zz[zz<1] <- 1
Pi <- n*zz/sum(zz)
Sel <- sort(sample(1:N, n, prob=zz/sum(zz)))   # PPS sampling
Ind <- rep(0,N)
Ind[Sel] <- 1
sY <- Y[Ind==1]
sX <- X[Ind==1,]
sPi <- Pi[Ind==1]


# HT estimator
HT <- sum(sY/sPi)/N
HTV <- n*var(sY/sPi)/N^2

# GREG estimator
greg1 <- GREG(sY, sX, sPi, mX, shrink="NO", samp="PPS")$Greg
greg2 <- GREG(sY, sX, sPi, mX, shrink="Lasso", samp="PPS")$Greg
greg3 <- GREG(sY, sX, sPi, mX, shrink="Ridge", samp="PPS")$Greg
greg4 <- GREG.forward(sY, sX, sPi, mX, samp="PPS")$Greg
greg5 <- GREG.mixed(sY, sX, sPi, mX, samp="PPS")$Greg

# Proposed methods
fit1 <- ABMA.LM(sY, sX, sPi, mX, N, samp="PPS", prior="default", mc=3000, burn=1000)
fit2 <- ABMA.LM(sY, sX, sPi, mX, N, samp="PPS", prior="Laplace", mc=3000, burn=1000)
fit3 <- ABMA.LM(sY, sX, sPi, mX, N, samp="PPS", prior="HS", mc=3000, burn=1000)


##  point estimates  ##
HT
greg1[1]
greg2[1]
greg3[1]
greg4[1]
greg5[1]
mean(fit1$Mu)
mean(fit2$Mu)
mean(fit3$Mu)
Mu  # true value 


##  95% confidence (credible) intervals  ##
zz <- qnorm(0.975)
HT+c(-1,1)*zz*sqrt(HTV)
greg1[1]+c(-1,1)*zz*sqrt(greg1[2])
greg2[1]+c(-1,1)*zz*sqrt(greg2[2])
greg3[1]+c(-1,1)*zz*sqrt(greg3[2])
greg4[1]+c(-1,1)*zz*sqrt(greg4[2])
greg5[1]+c(-1,1)*zz*sqrt(greg5[2])
quantile(fit1$Mu, prob=c(0.025, 0.975))
quantile(fit2$Mu, prob=c(0.025, 0.975))
quantile(fit3$Mu, prob=c(0.025, 0.975))

