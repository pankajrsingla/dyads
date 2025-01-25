p2 <- function (net, sender = NULL, receiver = NULL , density = NULL, reciprocity = NULL, burnin = NULL, sample = NULL, adapt= NULL, seed = NULL) 
{
  model <- p2model(net, sender, receiver, density, reciprocity)
  y <- model$y
  X <- model$X
  X1 <- model$X1
  X2 <- model$X2
  X3 <- model$X3
  X4 <- model$X4
  nact <- model$nact
  ns <- model$ns
  nre <- model$nre
  nd <- model$nd
  nr <- model$nr
  # number of parameters
  nb <- ns+ nre+ nd
  nrand <- nact*2
  npar <- nrand + nb
  nvarpar <- 2
  nvarcovpar <- nvarpar*(nvarpar+1)/2
  # sampling parameters
  if(!is.null(adapt)){
    Nadapt <- adapt
  } else {
    Nadapt <- 100
  } 
  Sadapt <- 125
  gacc <- round(Sadapt/3)
  if(!is.null(burnin)){
    Nburn <- burnin
  } else {
    Nburn <- 10000
  } 
  if(!is.null(burnin)){
    Nsamp <- sample
  } else {
    Nsamp <- 80000
  } 
  if(!is.null(seed)){
    set.seed(seed)  
  } else {
    set.seed(1)
  } 
  # prior distributions
  pSDb <- 10
  sigmab <- c(if (ns>0) {pSDb/apply(X1, 2, sd)}, if (nre>0) {pSDb/apply(X2, 2, sd)}, pSDb, if (nd>1) {pSDb/apply(X3[,2:nd, drop=F], 2, sd)})
  sigmaR <- 1
  sigmapar <- c(sigmab, rep(sigmaR, nrand))
  pmb <- as.vector(rep(0, npar))
  pVb <- diag(sigmapar^2)
  pmr <- as.vector(rep(0, nr))
  pSDg4 <- 10
  pScaler <- c(pSDg4, if (nr>1) {pSDg4/apply(X4[,2:nr, drop=F], 2, sd)})
  pVr <- diag(pScaler^2, ncol = nr)
  pmbr <- c(pmb[1:nb], pmr)
  pVbr <- diag(c(sigmapar[1:nb]^2, pScaler^2))
  pdfR <- nvarpar +1
  pVR <- diag(rep(nvarpar))
  # create vectors and matrices
  beta <- as.vector(rep(0, npar))
  g4 <- as.vector(rep(0, nr))
  varsimsAD <- matrix(rep(NA,nvarpar^2*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=nvarpar^2)
  bsimsAD <- matrix(rep(NA, npar*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=npar)
  rsimsAD <- matrix(rep(NA, nr*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=nr)
  varsimsBI <- matrix(rep(NA,nvarpar^2*Nburn), nrow= Nburn, ncol=nvarpar^2)
  bsimsBI <- matrix(rep(NA, npar*Nburn), nrow= Nburn, ncol=npar)
  rsimsBI <- matrix(rep(NA, nr*Nburn), nrow= Nburn, ncol=nr)
  var <- matrix(rep(NA,nvarpar^2*Nsamp), nrow= Nsamp, ncol=nvarpar^2)
  sims <- matrix(rep(NA, npar*Nsamp), nrow= Nsamp, ncol=npar)
  rsims <- matrix(rep(NA, nr*Nsamp), nrow= Nsamp, ncol=nr)
  Ml1 <- matrix(rep(0, nact*nact), ncol=nact)
  My <- matrix(rep(0, nact*nact), ncol=nact)
  My2 <- matrix(rep(0, nact*nact), ncol=nact)
  M <- matrix(rep(0, nact*nact), ncol=nact)
  M2 <- matrix(rep(0, nact*nact), ncol=nact)
  R <- matrix(rep(0, nact*nact), ncol=nact)
  R2 <- matrix(rep(0, nact*nact), ncol=nact)
  # calculations only needed once
  tX <- t(X)
  postdfR <- pdfR + nact
  varAD <- pVR
  IpVb <- diag(1/(sigmab^2), npar)
  rInd <- net*t(net)
  covRWADr <- pVr/1000
  covRWADb <- pVb/1000
  covRWADC <- pVR/100
  # Adaptation
  ll1 <- llp2(y, X, X4, beta, g4, M, My, R, rInd)
  c <- cbind(beta[c((nb[1]+1):(nb[1]+nact))], beta[c((nb[1]+nact+1):npar)])
  Sb <- 1
  Sc <- 1
  Sr <- 1
  #print(paste("adaptive sequence"))
  #pb <- txtProgressBar(min = 0, max = Nadapt, style = 3)
  for (i in 1:Nadapt){
    accb <- 0
    accc <- 0
    accr <- 0
    for (j in 1:Sadapt){ 
      # parameters
      beta2 <- beta 
      beta2[1:nb] <- beta[1:nb] + as.vector(Rfast::rmvnorm(1, pmb[1:nb], as.matrix(covRWADb[1:nb, 1:nb])))
      g4s <- g4 + as.vector(Rfast::rmvnorm(1, pmr, as.matrix(covRWADr[1:nr, 1:nr])))
      ll2 <- llp2(y, X, X4, beta2, g4s, M, My, R, rInd)
      ll1br <- ll1 + mvtnorm::dmvnorm(t(c(beta[1:nb], g4)), mean= pmbr, sigma= pVbr)    
      ll2br <- ll2 + mvtnorm::dmvnorm(t(c(beta2[1:nb], g4s)), mean= pmbr, sigma= pVbr) 
      if (runif(1, min = 0, max = 1) <  min(1, exp(ll2br-ll1br))){
        bsimsAD[((i-1)*Sadapt + j), ] <- beta2
        beta <- beta2
        rsimsAD[((i-1)*Sadapt + j),] <- g4s
        g4 <- g4s
        ll1 <- ll2
        accb <- accb + 1
        accr <- accr + 1
      } else {
        bsimsAD[((i-1)*Sadapt + j), ] <- beta
        rsimsAD[((i-1)*Sadapt + j),] <- g4
      }   
      # random effects
      beta2 <- beta
      beta2[(nb+1):npar] <- beta[(nb+1):npar] + as.vector(Rfast::rmvnorm(nact, c(0,0), covRWADC))
      ll2 <- llp2(y, X, X4, beta2, g4, M, My, R, rInd)
      ll1C <- ll1 + sum(log(mvtnorm::dmvnorm(c, mean= c(0,0), sigma= varAD)))   
      c2 <- cbind(beta2[c((nb[1]+1):(nb[1]+nact))], beta2[c((nb[1]+nact+1):npar)])
      ll2C <- ll2 + sum(log(mvtnorm::dmvnorm(c2, mean= c(0,0), sigma= varAD)))
      if (runif(1, min = 0, max = 1) <  min(1, exp(ll2C-ll1C))){
        bsimsAD[((i-1)*Sadapt + j), ] <- beta2
        beta <- beta2
        ll1 <- ll2  
        accc <- accc +1
        c <- c2
      } else {
        bsimsAD[((i-1)*Sadapt + j), ] <- beta
      }
      # Sigma
      varAD <- MASS::ginv(matrix(stats::rWishart(1, postdfR, MASS::ginv(t(c)%*%c + pVR)), ncol=nvarpar))
      varsimsAD[i,] <- as.vector(varAD)
    }  
    if (accb > gacc){
      Sb <- Sb*(1+(1-(Sadapt-accb)/(Sadapt-gacc)))  
    } else {
      Sb <- Sb/(1+(1-(accb/gacc)))
    }
    covRWADb <- Sb*cov(bsimsAD, use= "complete.obs")
    if (accc > gacc){
      Sc <- Sc*(1+(1-(Sadapt-accc)/(Sadapt-gacc)))  
    } else {
      Sc <- Sc/(1+(1-(accc/gacc)))
    }
    covRWADC <- Sc*matrix(colMeans(varsimsAD, na.rm = T), ncol=2)
    if (accr > gacc){
      Sr <- Sr*(1+(1-(Sadapt-accr)/(Sadapt-gacc)))  
    } else {
      Sr <- Sr/(1+(1-(accr/gacc)))
    }
    covRWADr <- Sr*cov(rsimsAD, use= "complete.obs")
    # update progress bar
    #Sys.sleep(0.1)
    #setTxtProgressBar(pb, i)
  } 
  #close(pb)
  # Burn in
  bsimsBI[1,] <- beta
  rsimsBI[1,] <- g4
  c <- cbind(beta[c((nb[1]+1):(nb[1]+nact))], beta[c((nb[1]+nact+1):npar)])
  varsimsBI[1,] <- as.vector(varAD)
  varBI <- varAD
  covRWb <- covRWADb
  covRWC <- covRWADC
  covRWr <- covRWADr
  #print(paste("burn-in"))
  #pb <- txtProgressBar(min = 0, max = Nburn, style = 3)
  for (i in 2:Nburn){
    # parameters 
    beta2 <- beta 
    beta2[1:nb] <- beta[1:nb] + as.vector(Rfast::rmvnorm(1, pmb[1:nb], as.matrix(covRWADb[1:nb, 1:nb])))
    g4s <- g4 + as.vector(Rfast::rmvnorm(1, pmr, as.matrix(covRWADr[1:nr, 1:nr])))
    ll2 <- llp2(y, X, X4, beta2, g4s, M, My, R, rInd)
    ll1br <- ll1 + log(mvtnorm::dmvnorm(t(c(beta[1:nb], g4)), mean= pmbr, sigma= pVbr))    
    ll2br <- ll2 + log(mvtnorm::dmvnorm(t(c(beta2[1:nb], g4s)), mean= pmbr, sigma= pVbr))    
    if (runif(1, min = 0, max = 1) <  min(1, exp(ll2br-ll1br))){
      bsimsBI[i,] <- beta2
      beta <- beta2
      rsimsBI[i,] <- g4s
      g4 <- g4s
      ll1 <- ll2
    } else {
      bsimsBI[i,] <- beta
      rsimsBI[i,] <- g4
    }   
    # random effects
    beta2 <- beta
    beta2[(nb+1):npar] <- beta[(nb+1):npar] + as.vector(Rfast::rmvnorm(nact, c(0,0), covRWADC))
    ll2 <- llp2(y, X, X4, beta2, g4, M, My, R, rInd)
    ll1C <- ll1 + sum(log(mvtnorm::dmvnorm(c, mean= c(0,0), sigma= varBI)))   
    c2 <- cbind(beta2[c((nb[1]+1):(nb[1]+nact))], beta2[c((nb[1]+nact+1):npar)])
    ll2C <- ll2 + sum(log(mvtnorm::dmvnorm(c2, mean= c(0,0), sigma= varBI)))
    if (runif(1, min = 0, max = 1) <  min(1, exp(ll2C-ll1C))){
      bsimsBI[i,] <- beta2
      beta <- beta2
      ll1 <- ll2  
      c <- c2
    } else {
      bsimsBI[i,] <- beta
    } 
    # Sigma
    varBI <- MASS::ginv(matrix(stats::rWishart(1, postdfR, MASS::ginv(t(c)%*%c + pVR)), ncol=nvarpar))
    varsimsBI[i,] <- as.vector(varBI)
    # update progress bar
    #Sys.sleep(0.1)
    #setTxtProgressBar(pb, i)
  } 
  #close(pb)
  # Sample
  sims[1,] <- beta
  rsims[1,] <- g4
  c <- cbind(beta[c((nb[1]+1):(nb[1]+nact))], beta[c((nb[1]+nact+1):npar)])
  var[1,] <- as.vector(varBI)
  vartmp <- varBI
  #print(paste("sample"))
  #pb <- txtProgressBar(min = 0, max = Nsamp, style = 3)
  for (i in 2:Nsamp){
    # parameters
    beta2 <- beta 
    beta2[1:nb] <- beta[1:nb] + as.vector(Rfast::rmvnorm(1, pmb[1:nb], as.matrix(covRWADb[1:nb, 1:nb])))
    g4s <- g4 + as.vector(Rfast::rmvnorm(1, pmr, as.matrix(covRWADr[1:nr, 1:nr])))
    ll2 <- llp2(y, X, X4, beta2, g4s, M, My, R, rInd)
    ll1br <- ll1 + log(mvtnorm::dmvnorm(t(c(beta[1:nb], g4)), mean= pmbr, sigma= pVbr))    
    ll2br <- ll2 + log(mvtnorm::dmvnorm(t(c(beta2[1:nb], g4s)), mean= pmbr, sigma= pVbr))    
    if (runif(1, min = 0, max = 1) <  min(1, exp(ll2br-ll1br))){
      sims[i,] <- beta2
      beta <- beta2
      rsims[i,] <- g4s
      g4 <- g4s
      ll1 <- ll2
    } else {
      sims[i,] <- beta
      rsims[i,] <- g4
    }   
    # random effects
    beta2 <- beta
    beta2[(nb+1):npar] <- beta[(nb+1):npar] + as.vector(Rfast::rmvnorm(nact, c(0,0), covRWADC))
    ll2 <- llp2(y, X, X4, beta2, g4, M, My, R, rInd)
    ll1C <- ll1 + sum(log(mvtnorm::dmvnorm(c, mean= c(0,0), sigma= vartmp)))   
    c2 <- cbind(beta2[c((nb[1]+1):(nb[1]+nact))], beta2[c((nb[1]+nact+1):npar)])
    ll2C <- ll2 + sum(log(mvtnorm::dmvnorm(c2, mean= c(0,0), sigma= vartmp)))
    if (runif(1, min = 0, max = 1) <  min(1, exp(ll2C-ll1C))){
      sims[i,] <- beta2
      beta <- beta2
      ll1 <- ll2  
      c <- c2
    } else {
      sims[i,] <- beta
    } 
    # Sigma
    vartmp <- MASS::ginv(matrix(stats::rWishart(1, postdfR, MASS::ginv(t(c)%*%c + pVR)), ncol=nvarpar))
    var[i,] <- as.vector(vartmp)
    # update progress bar
    #Sys.sleep(0.1)
    #setTxtProgressBar(pb, i)
  } 
  #close(pb) 
  output.matrix <- matrix(NA, nvarcovpar+nb+nr, 10)
  collabels <- c("Estimate", "SE", "Q.05", "Q2.5", "Q25", "Q50", "Q75", "Q97.5", "Q99.5", "Neff")
  colnames(output.matrix) <- collabels
  rowlabels <- c("sender variance", "sender receiver covariance", "receiver variance", all.vars(sender),
                 all.vars(receiver), "density", all.vars(density), "reciprocity", all.vars(reciprocity))
  rownames(output.matrix) <- rowlabels
  output.matrix[,1] <- c(colMeans(var)[c(1,2,4)],if (nb>1) {colMeans(sims[,1:nb])} else {mean(sims[,1])}, if (nr>1) {colMeans(rsims[, 1:nr])} else {mean(rsims)}) 
  output.matrix[,2] <- c(sqrt(diag(cov(var)))[c(1,2,4)], if (nb>1) {sqrt(diag(cov(sims[, 1:nb])))} else {sd(sims)}, if (nr>1) {sqrt(diag(cov(rsims)))} else {sd(rsims)}) 
  output.matrix[1,3:9] <- quantile(var[,1],  probs = c(0.5, 2.5, 25, 50, 75, 97.5, 99.5)/100)
  output.matrix[,3:9] <- t(apply(cbind(var[, c(1,2,4)], sims[,1:nb], rsims), 2, quantile, probs = c(0.5, 2.5, 25, 50, 75, 97.5, 99.5)/100))
  output.matrix[,10] <- t(apply(cbind(var[, c(1,2,4)], sims[,1:nb], rsims), 2, effectiveEst))
  return(round(output.matrix, digits=3))
}