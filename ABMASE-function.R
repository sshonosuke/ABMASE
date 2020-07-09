## inastall packages
library(statmod)
library(glmnet)
library(MASS)
library(MCMCpack)




###    Approximate Bayesian approach (linear working model)    ###
# sY: sampled response variable
# sX: (n,p)-matrix of covariates without intercept
# sPi: sampling probability
# mX: population mean vector of covariates without intercept
# N: population size
# samp: sampling mechanism ("SRS", "PPS" or "Pois") 
# prior: prior for regression coefficients ("default", "Laplace" or "HS")
# mc: number of Monte Carlo samples
# burn: length of burn-in period in MCMC
# a: shape parameter in the hyperprior of lambda in the Laplace prior

ABMA.LM <- function(sY, sX, sPi, mX, N, samp="SRS", prior="default", mc=5000, burn=1000, a=NULL){
  p <- dim(sX)[2]
  Gam <- diag(1/sPi)
  XX <- cbind(1, sX)
  est <- as.vector( ginv(t(XX)%*%Gam%*%XX)%*%t(XX)%*%Gam%*%sY )   
  halpha <- est[1]  # point estimate of intercept
  hbeta <- est[-1]  # point estimates of regression coefficients
  resid <- as.vector( sY-XX%*%est )   # residual 
  
  ## Delta matrix for each sampling mechanism
  if(samp=="Pois"){ mat1 <- t(XX)%*%diag(resid^2*(1-sPi)/sPi^2)%*%XX }
  if(samp=="SRS"){
    Delta = matrix(NA,n,n)
    for(i in 1:n){
      Delta[i,i] = 1-n/N
      pp2 = n*(n-1)/(N*(N-1))
      Delta[i,-i] = 1-(n/N)^2/pp2
    }
    mat1 <- t(XX*resid/sPi)%*%Delta%*%(XX*resid/sPi)
  }
  if(samp=="PPS"){
    UU <- XX*resid/sPi
    sUU <- scale(t(UU), scale=F)
    mat1 <- n*sUU%*%t(sUU)/(n-1)
  }
  
  mat2 <- ginv( t(XX)%*%Gam%*%XX )
  V <- mat2%*%mat1%*%mat2      # covariance matrix of beta 
  Vb <- V[-1, -1]
  Va <- V[1, 1]
  Vab <- V[1, -1]
  
  ## Posterior of beta (default non-informative prior)
  if(prior=="default"){
    Beta.pos <- mvrnorm(mc, hbeta, Vb)
    Aux <- list()
  }
  
  ## Posterior of beta (Lasso prior)
  if(prior=="Laplace"){
    if(is.null(a)){
      opt.lam <- cv.glmnet(sX,sY,family="gaussian", alpha=1, weights=1/sPi, nfolds=10)$lambda.min
      a <- opt.lam^2    # hyperparameter (shape) in gamma prior for lambda
    }   
    b <- 1      # hyperparameter (scale) in gamma prior for lambda
    mc2 <- mc+burn
    Prec <- solve(Vb)   # precision matrix for beta
    vAlpha <- max(as.vector( Va-t(Vab)%*%solve(Vb)%*%Vab ), 0)   # varinace for alpha
    
    Beta.pos <- matrix(NA, mc2, p)
    Alpha.pos <- c()
    Tau2.pos <- matrix(NA, mc2, p)
    Lam.pos <- c()
    
    # Initial values
    Beta <- hbeta
    Alpha <- halpha
    Tau2 <- hbeta^2
    Lam <- 1
    
    # MCMC 
    for(k in 1:mc2){
      # Beta
      invA <- ginv(Prec + diag(1/Tau2))
      mBeta <- as.vector( invA%*%Prec%*%hbeta )
      Beta <- mvrnorm(1, mBeta, invA)
      Beta.pos[k,] <- Beta
      # Alpha
      mAlpha <- halpha - as.vector( Vab%*%solve(Vb)%*%(hbeta-Beta) )
      Alpha <- rnorm(1, mAlpha, sqrt(vAlpha))
      Alpha.pos[k] <- Alpha
      # Tau and Lam
      mu <- sqrt(Lam^2/Beta^2)
      Tau2 <- 1/rinvgauss(p, mu, Lam^2)
      Lam <- sqrt( rgamma(1, p+a, sum(Tau2)/2+b) )
      Tau2.pos[k,] <- Tau2
      Lam.pos[k] <- Lam
    }
    # result
    om <- 1:burn
    Beta.pos <- Beta.pos[-om,]
    Alpha.pos <- Alpha.pos[-om]
    Tau2.pos <- Tau2.pos[-om,]
    Lam.pos <- Lam.pos[-om]
    Aux <- list(Alpha=Alpha.pos, Tau2=Tau2.pos, Lam=Lam.pos)   # auxiliary results
  }
  
  ##  Posterior of beta (Horseshoe prior)
  if(prior=="HS"){    
    mc2 <- mc+burn
    Prec <- solve(Vb)   # precision matrix for beta
    vAlpha <- max(as.vector( Va-t(Vab)%*%solve(Vb)%*%Vab ), 0)   # varinace for alpha
    
    Beta.pos <- matrix(NA, mc2, p)
    Alpha.pos <- c()    
    U2.pos <- matrix(NA, mc2, p)
    Lam2.pos <- c()
    Xi.pos <- matrix(NA, mc2, p)
    Gam.pos <- c()
    
    # Initial values
    Beta <- hbeta
    Alpha <- halpha
    U2 <- hbeta^2
    Lam2 <- 1
    Xi <- rep(1, p)
    Gam <- 1
    
    # MCMC
    for(k in 1:mc2){
      # Beta
      invA <- ginv( Prec+diag(U2*Lam2) ) 
      mBeta <- invA%*%Prec%*%hbeta
      Beta <- mvrnorm(1, mBeta, invA)
      Beta.pos[k,] <- Beta
      # Alpha
      mAlpha <- halpha - as.vector( Vab%*%solve(Vb)%*%(hbeta-Beta) )
      Alpha <- rnorm(1, mAlpha, sqrt(vAlpha))
      Alpha.pos[k] <- Alpha
      # U2
      bb <- 0.5*Beta^2/Lam2 + 1/Xi
      U2 <- rinvgamma(p, 1, bb)
      U2.pos[k, ] <- U2
      # Xi
      Xi <- rinvgamma(p, 1, 1+1/U2)
      Xi.pos[k,] <- Xi
      # Lam2
      cc <- 0.5*sum(Beta^2/U2) + 1/Gam
      Lam2 <- rinvgamma(1, (p+1)/2, cc)
      Lam2.pos[k] <- Lam2
      # Gam
      Gam <- max(0.001, rinvgamma(1, 1, 1/Lam2))
      Gam.pos[k] <- Gam
    }
    # result
    om <- 1:burn
    Beta.pos <- Beta.pos[-om,]
    Alpha.pos <- Alpha.pos[-om]
    U2.pos <- U2.pos[-om,]
    Lam2.pos <- Lam2.pos[-om]
    Aux <- list(Alpha=Alpha.pos, U2=U2.pos, Lam2=Lam2.pos)   # auxiliary results
  }
  
  ##  Posterior of population mean    
  Resid.pos <- sY-sX%*%t(Beta.pos)
  mm.pos <- Beta.pos%*%mX + apply(Resid.pos/sPi,2,sum)/N
  
  if(samp=="Pois"){ ss.pos <- apply(Resid.pos^2*(1-sPi)/sPi^2/N^2, 2, sum) }
  if(samp=="SRS"){ ss.pos <- diag( t(Resid.pos/sPi)%*%Delta%*%(Resid.pos/sPi)/N^2 ) }
  if(samp=="PPS"){ ss.pos <- n*apply(Resid.pos/sPi, 2, var)/N^2  }
  
  Mu.pos <- rnorm(mc, mm.pos, sqrt(ss.pos))
  
  ## Summary
  Res <- list(Mu=Mu.pos, Beta=Beta.pos, Aux=Aux)
  return(Res)
}




















###    Approximate Bayesian approach (linear working model)    ###
# sY: sampled response variable
# sX: (n,p)-matrix of sampled covariates without intercept
# sPi: sampling probability
# tX: (N,p)-matrix of population covariates without intercept
# samp: sampling mechanism ("SRS", "PPS" or "Pois") 
# prior: prior for regression coefficients ("default", "Laplace" or "HS")
# mc: number of Monte Carlo samples
# burn: length of burn-in period in MCMC
# a: shape parameter in the hyperprior of lambda in the Laplace prior

ABMA.BIN <- function(sY, sX, sPi, tX, samp="SRS", prior="default", mc=5000, burn=1000, a=NULL){
  N <- dim(tX)[1]
  p <- dim(sX)[2]
  Gam <- diag(1/sPi)
  XX <- cbind(1, sX)
  
  ## function for weighted logistic regression  
  NR.wlogit <- function(sY, sX){
    est <- coef( glm(sY~sX, family="binomial") )
    sXX <- cbind(1, sX)
    maxitr <- 100
    for(k in 1:maxitr){
      est0 <- est
      pr <- logistic( as.vector(sXX%*%est) )
      vv <- pr*(1-pr)/sPi
      est <- est + ginv( t(sXX)%*%diag(vv)%*%sXX)%*%t(sXX)%*%((sY-pr)/sPi )
      est <- as.vector(est)
      dd <- sum( abs(est-est0) )
      if( dd < 10^(-10) ){ break() }
    }
    return(est)
  }
  
  ## point estimates
  est <- NR.wlogit(sY, sX)
  mm <- as.vector( logistic(cbind(1,sX)%*%est) )
  resid <- sY-mm
  halpha <- est[1]
  hbeta <- est[-1]
  
  ## Delta matrix for each sampling mechanism
  if(samp=="Pois"){ mat1 <- t(XX)%*%diag(resid^2*(1-sPi)/sPi^2)%*%XX }
  if(samp=="SRS"){
    Delta = matrix(NA,n,n)
    for(i in 1:n){
      Delta[i,i] = 1-n/N
      pp2 = n*(n-1)/(N*(N-1))
      Delta[i,-i] = 1-(n/N)^2/pp2
    }
    mat1 <- t(XX*resid/sPi)%*%Delta%*%(XX*resid/sPi)
  }
  if(samp=="PPS"){
    UU <- XX*resid/sPi
    sUU <- scale(t(UU), scale=F)
    mat1 <- n*sUU%*%t(sUU)/(n-1)
  }
  
  Gam2 <- diag( mm*(1-mm)/sPi )
  mat2 <- ginv(t(XX)%*%Gam2%*%XX)
  V <- mat2%*%mat1%*%mat2     # covariance matrix of beta 
  Vb <- V[-1, -1]
  Va <- V[1, 1]
  Vab <- V[1, -1]
  
  ## Posterior samples (Default prior)
  if(prior=="default"){
    Beta.pos <- mvrnorm(mc, hbeta, Vb)
    Alpha.pos <- rep(halpha, mc)
    Aux <- list()
  }
  
  ## Posterior of beta (Lasso prior)
  if(prior=="Laplace"){
    if(is.null(a)){
      opt.lam <- cv.glmnet(sX,sY,family="binomial", alpha=1, weights=1/sPi, nfolds=10)$lambda.min
      a <- opt.lam^2    # hyperparameter (shape) in gamma prior for lambda
    }   
    b <- 1      # hyperparameter (scale) in gamma prior for lambda
    mc2 <- mc+burn
    Prec <- solve(Vb)   # precision matrix for beta
    vAlpha <- max(as.vector( Va-t(Vab)%*%solve(Vb)%*%Vab ), 0)   # varinace for alpha
    
    Beta.pos <- matrix(NA, mc2, p)
    Alpha.pos <- c()
    Tau2.pos <- matrix(NA, mc2, p)
    Lam.pos <- c()
    
    # Initial values
    Beta <- hbeta
    Alpha <- halpha
    Tau2 <- hbeta^2
    Lam <- 1
    
    # MCMC 
    for(k in 1:mc2){
      # Beta
      invA <- ginv(Prec + diag(1/Tau2))
      mBeta <- as.vector( invA%*%Prec%*%hbeta )
      Beta <- mvrnorm(1, mBeta, invA)
      Beta.pos[k,] <- Beta
      # Alpha
      mAlpha <- halpha - as.vector( Vab%*%solve(Vb)%*%(hbeta-Beta) )
      Alpha <- rnorm(1, mAlpha, sqrt(vAlpha))
      Alpha.pos[k] <- Alpha
      # Tau and Lam
      mu <- sqrt(Lam^2/Beta^2)
      Tau2 <- 1/rinvgauss(p, mu, Lam^2)
      Lam <- sqrt( rgamma(1, p+a, sum(Tau2)/2+b) )
      Tau2.pos[k,] <- Tau2
      Lam.pos[k] <- Lam
    }
    # result
    om <- 1:burn
    Beta.pos <- Beta.pos[-om,]
    Alpha.pos <- Alpha.pos[-om]
    Tau2.pos <- Tau2.pos[-om,]
    Lam.pos <- Lam.pos[-om]
    Aux <- list(Alpha=Alpha.pos, Tau2=Tau2.pos, Lam=Lam.pos)   # auxiliary results
  }
  
  ##  Posterior of beta (Horseshoe prior)
  if(prior=="HS"){    
    mc2 <- mc+burn
    Prec <- solve(Vb)   # precision matrix for beta
    vAlpha <- max(as.vector( Va-t(Vab)%*%solve(Vb)%*%Vab ), 0)   # varinace for alpha
    
    Beta.pos <- matrix(NA, mc2, p)
    Alpha.pos <- c()    
    U2.pos <- matrix(NA, mc2, p)
    Lam2.pos <- c()
    Xi.pos <- matrix(NA, mc2, p)
    Gam.pos <- c()
    
    # Initial values
    Beta <- hbeta
    Alpha <- halpha
    U2 <- hbeta^2
    Lam2 <- 1
    Xi <- rep(1, p)
    Gam <- 1
    
    # MCMC
    for(k in 1:mc2){
      # Beta
      invA <- ginv( Prec+diag(U2*Lam2) ) 
      mBeta <- invA%*%Prec%*%hbeta
      Beta <- mvrnorm(1, mBeta, invA)
      Beta.pos[k,] <- Beta
      # Alpha
      mAlpha <- halpha - as.vector( Vab%*%solve(Vb)%*%(hbeta-Beta) )
      Alpha <- rnorm(1, mAlpha, sqrt(vAlpha))
      Alpha.pos[k] <- Alpha
      # U2
      bb <- 0.5*Beta^2/Lam2 + 1/Xi
      U2 <- rinvgamma(p, 1, bb)
      U2.pos[k, ] <- U2
      # Xi
      Xi <- rinvgamma(p, 1, 1+1/U2)
      Xi.pos[k,] <- Xi
      # Lam2
      cc <- 0.5*sum(Beta^2/U2) + 1/Gam
      Lam2 <- rinvgamma(1, (p+1)/2, cc)
      Lam2.pos[k] <- Lam2
      # Gam
      Gam <- max(0.001, rinvgamma(1, 1, 1/Lam2))
      Gam.pos[k] <- Gam
    }
    # result
    om <- 1:burn
    Beta.pos <- Beta.pos[-om,]
    Alpha.pos <- Alpha.pos[-om]
    U2.pos <- U2.pos[-om,]
    Lam2.pos <- Lam2.pos[-om]
    Aux <- list(Alpha=Alpha.pos, U2=U2.pos, Lam2=Lam2.pos)   # auxiliary results
  }
  
  ##   Posterior of population mean    
  Resid.pos <- sY-logistic( halpha + sX%*%t(Beta.pos) )
  mm.pos <- apply(logistic(halpha + tX%*%t(Beta.pos)), 2, mean) + apply(Resid.pos/sPi,2,sum)/N
  
  if(samp=="Pois"){ ss.pos <- apply(Resid.pos^2*(1-sPi)/sPi^2/N^2, 2, sum) }
  if(samp=="SRS"){ ss.pos <- diag( t(Resid.pos/sPi)%*%Delta%*%(Resid.pos/sPi)/N^2 ) }
  if(samp=="PPS"){ ss.pos <- n*apply(Resid.pos/sPi, 2, var)/N^2  }
  
  Mu.pos <- rnorm(mc, mm.pos, sqrt(ss.pos))
  
  ## Summary
  Res <- list(Mu=Mu.pos, Beta=Beta.pos, Aux=Aux)
  return(Res)
}



