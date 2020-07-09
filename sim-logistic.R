set.seed(2)
source("ABMASE-function.R")

logistic <- function(x){ exp(x)/(1+exp(x)) }

# settings
p <- 30    # maximum number of covariates
N <- 10000    # population size
n <- 400    # sample size

# covariates
rho <- 0.2
C <- matrix(NA,p,p)
for(j in 1:p){ 
  for(k in 1:p){ C[k,j] <- rho^(abs(j-k)) }
}
D <- diag( rep(sqrt(2), p) )
X <- mvrnorm(N, rep(1,p), D%*%C%*%D)
tX <- X    # population covariate vector


# true coeffieicnts and varaince 
Beta <- rep(0, p)
Beta[c(1, 4, 7, 10)] <- c(1, -0.5, 1, -0.5)
sig <- 2

# response
pp <- logistic( as.vector(-1+X%*%Beta) )
Y <- rbinom(N, 1, pp)

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


# Proposed methods
fit1 <- ABMA.BIN(sY, sX, sPi, tX, N, samp="SRS", prior="default", mc=3000, burn=1000)
fit2 <- ABMA.BIN(sY, sX, sPi, tX, N, samp="SRS", prior="Laplace", mc=3000, burn=1000)
fit3 <- ABMA.BIN(sY, sX, sPi, tX, N, samp="SRS", prior="HS", mc=3000, burn=1000)


##  point estimates  ##
HT
mean(fit1$Mu)
mean(fit2$Mu)
mean(fit3$Mu)
Mu  # true value 


##  95% confidence (credible) intervals  ##
zz <- qnorm(0.975)
HT+c(-1,1)*zz*sqrt(HTV)
quantile(fit1$Mu, prob=c(0.025, 0.975))
quantile(fit2$Mu, prob=c(0.025, 0.975))
quantile(fit3$Mu, prob=c(0.025, 0.975))





###----------------------------------------------------------------###
###    Scenario 2 (probability-proportional-to-size sampling)      ###
###----------------------------------------------------------------###
zz <- log(1+abs(0.5*Y+rexp(N,3)))   # size variable
zz[ zz<0.5 ] <- 0.5
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

# Proposed methods
fit1 <- ABMA.BIN(sY, sX, sPi, tX, N, samp="PPS", prior="default", mc=3000, burn=1000)
fit2 <- ABMA.BIN(sY, sX, sPi, tX, N, samp="PPS", prior="Laplace", mc=3000, burn=1000)
fit3 <- ABMA.BIN(sY, sX, sPi, tX, N, samp="PPS", prior="HS", mc=3000, burn=1000)


##  point estimates  ##
HT
mean(fit1$Mu)
mean(fit2$Mu)
mean(fit3$Mu)
Mu  # true value 


##  95% confidence (credible) intervals  ##
zz <- qnorm(0.975)
HT+c(-1,1)*zz*sqrt(HTV)
quantile(fit1$Mu, prob=c(0.025, 0.975))
quantile(fit2$Mu, prob=c(0.025, 0.975))
quantile(fit3$Mu, prob=c(0.025, 0.975))

