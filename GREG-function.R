## inastall packages
library(statmod)
library(glmnet)
library(MASS)
library(MCMCpack)
library(olsrr)




Tr <- function(A){ sum(diag(A)) }


GREG <- function(sY, sX, sPi, mX, shrink="NO", samp="SRS"){
  pp = dim(sX)[2]
  nf = 10
  if(shrink=="NO"){
    fit = lm(sY~-1+sX,weights=1/sPi)
    hbeta = coef(fit)
    Resid = as.vector(sY-predict(fit))
    est = sum(na.omit(mX*hbeta))+sum(Resid/sPi)/N
    Aux = hbeta
    names(Aux) = paste0("Beta",1:pp)
  }
  if(shrink=="Lasso"){
    opt.lam = cv.glmnet(sX,sY,family="gaussian",alpha=1,weights=1/sPi,nfolds=nf)$lambda.min
    Lasso = glmnet(sX,sY,family="gaussian",alpha=1,weights=1/sPi,lambda=opt.lam)
    hbeta = coef(Lasso)
    Resid = as.vector(sY-cbind(1,sX)%*%hbeta)
    est = as.vector(c(1,mX)%*%hbeta+sum(Resid/sPi)/N)
    Aux = c(opt.lam,as.vector(hbeta))
    names(Aux) = c("Lam",paste0("Beta",0:pp))
  }
  if(shrink=="Ridge"){
    opt.lam = cv.glmnet(sX,sY,family="gaussian",alpha=0,weights=1/sPi,nfolds=nf)$lambda.min
    Lasso = glmnet(sX,sY,family="gaussian",alpha=0,weights=1/sPi,lambda=opt.lam)
    hbeta = coef(Lasso)
    Resid = as.vector(sY-cbind(1,sX)%*%hbeta)
    est = as.vector(c(1,mX)%*%hbeta+sum(Resid/sPi)/N)
    Aux = c(opt.lam,as.vector(hbeta))
    names(Aux) = c("Lam", paste0("Beta",0:pp))
  }
  
  # variance estimation
  if(samp=="Pois"){  VV = sum(Resid^2*(1-sPi)/sPi^2)/N^2 }
  if(samp=="SRS"){ 
    Delta = matrix(NA,n,n)
    for(i in 1:n){
      Delta[i,i] = 1-n/N
      pp2 = n*(n-1)/(N*(N-1))
      Delta[i,-i] = 1-(n/N)^2/pp2
    }
    VV = t(Resid/sPi)%*%Delta%*%(Resid/sPi)/N^2
  }
  if(samp=="PPS"){
    VV = n*var(Resid/sPi)/N^2
  }
  
  # results
  result=list(Greg=c(est,VV), Aux=Aux)
  return(result)
}





# GREG: forward selection 
GREG.forward <- function(sY, sX, sPi, mX, samp="SRS"){
  p <- dim(sX)[2]
  model <- lm(sY~., data=data.frame(sY, sX), weights=1/sPi)
  VS <- ols_step_forward_p(model)
  Ind <- sort( as.numeric( factor(VS$predictors, levels=paste0("X",1:p)) ) )
  hbeta <- rep(0, p)
  hbeta[Ind] <- coef(lm(sY~-1+sX[,Ind], weights=1/sPi))
  pred <- as.vector(sX%*%hbeta)
  Resid <- sY-pred
  est <- sum(na.omit(mX*hbeta))+sum(Resid/sPi)/N
  
  # variance estimation
  if(samp=="Pois"){  VV = sum(Resid^2*(1-sPi)/sPi^2)/N^2 }
  if(samp=="SRS"){ 
    Delta = matrix(NA,n,n)
    for(i in 1:n){
      Delta[i,i] = 1-n/N
      pp2 = n*(n-1)/(N*(N-1))
      Delta[i,-i] = 1-(n/N)^2/pp2
    }
    VV = t(Resid/sPi)%*%Delta%*%(Resid/sPi)/N^2
  }
  if(samp=="PPS"){
    VV = n*var(Resid/sPi)/N^2
  }
  
  # results
  result <- list(Greg=c(est,VV), Aux=hbeta)
  return(result)
}




# GREG: mixed model approach
GREG.mixed <- function(sY, sX, sPi, mX, samp="SRS"){
  p <- dim(sX)[2]
  DD <- diag(sPi)
  invDD <- diag(1/sPi)
  mat <- sX%*%t(sX)
  
  Q <- function(vv){
    psi <- vv[1]
    sig2 <- vv[2]
    Sig <- psi*mat + sig2*DD
    val <- sum(log( eigen(Sig)$value )) + as.vector( t(sY)%*%solve(Sig)%*%sY )
    return(val)
  }
  
  est.vv <- optim(par=c(0.1, mean(1/sPi)), fn=Q, method="L-BFGS-B", lower=rep(0.001, 2), upper=rep(1000, 2))$par
  hpsi <- est.vv[1]/est.vv[2]
  mat <- solve( t(sX)%*%invDD%*%sX+hpsi*diag(p) )
  hbeta <- as.vector( mat%*%t(sX)%*%invDD%*%sY )
  
  # point estimate
  pred <- as.vector(sX%*%hbeta)
  Resid <- sY-pred
  est <- sum(mX*hbeta) + sum(Resid/sPi)/N
  
  # variance estimation 
  H <- mat%*%t(sX)%*%invDD%*%sX
  alpha <- 1/(N*sPi)
  XHT <- apply(sX/sPi, 2, sum)/N
  ww <- alpha + as.vector( t(mX-XHT)%*%H%*%solve(t(sX)%*%invDD%*%sX)%*%t(sX) )/sPi
  
  if(samp=="SRS"){ 
    Delta = matrix(NA,n,n)
    for(i in 1:n){
      Delta[i,i] = 1-n/N
      pp2 = n*(n-1)/(N*(N-1))
      Delta[i,-i] = 1-(n/N)^2/pp2
    }
    V1 <- t(Resid*ww)%*%Delta%*%(Resid*ww)
    mat1 <- t(sX*Resid*ww)%*%Delta%*%(sX*Resid*ww)*N^2
    VXH <- t(sX/sPi)%*%Delta%*%(sX/sPi)/N^2
  }
  
  if(samp=="PPS"){
    V1 = n*var(Resid*ww)
    UU <- sX*Resid*ww*N
    sUU <- scale(t(UU), scale=F)
    mat1 <- n*sUU%*%t(sUU)/(n-1)
    VXH <- n*var(sX/sPi)/N^2
  }
  
  mat2 <- ginv( t(sX)%*%invDD%*%sX )
  V.beta <- mat2%*%mat1%*%mat2
  GG <- solve(t(sX)%*%invDD%*%sX + diag(p)/hpsi)%*%t(sX)%*%invDD%*%sX - diag(p)
  V2 <- Tr( (hbeta%*%t(hbeta) - V.beta)%*%t(GG)%*%VXH%*%GG )

  VV <- V1 + V2
  
  # results
  result <- list(Greg=c(est,VV), Aux=hbeta)
  return(result)
}


