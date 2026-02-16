p2ML <- function (nets, sender = NULL, receiver = NULL, density =~ 1, reciprocity =~ 1, adapt = NULL, burnin = NULL, center = NULL, separateSigma = NULL, seed = NULL, densVar = NULL, recVar = NULL) 
{
  # sampling parameters
  if(!is.null(adapt)){
    Nadapt <- adapt
  } else {
    Nadapt <- 250
  } 
  nccyc <- 1
  Sadapt <- 125
  gacc <- round(Sadapt*.234)
  gaccC <- round(Sadapt*.234)
  gaccM <- round(Sadapt*.44)
  gaccMt <- round(Sadapt*.234)
  gaccRt <- round(Sadapt*.234)
  gaccR <- round(Sadapt*.44)
  if(!is.null(burnin)){
    NburnR <- burnin
  } else {
    NburnR <- 2500
  }
  if(NburnR < 2){
    NburnR <- 1
  }
  burn <- -1*seq(NburnR)
  if(!is.null(center)){
    center <- center
  } else {
    center <- TRUE
  } 
  if(!is.null(separateSigma)){
    separateSigma <- separateSigma
  } else {
    separateSigma <- FALSE
  } 
  if(!is.null(densVar)){
    densVar <- densVar
  } else {
    densVar <- TRUE
  } 
  if(!is.null(recVar)){
    recVar <- recVar
  } else {
    recVar <- TRUE
  } 
  if(!is.null(seed)){
    set.seed(seed)
    RcppZiggurat::zsetseed(seed)
  } else {
    set.seed(1)
    RcppZiggurat::zsetseed(1)
  } 
  # obtain model
  nnets <- length(nets)
  model <- modelMLB19lb(nets, sender, receiver, density, reciprocity, center, densVar, recVar)
  yl <- model$yl
  Xl <- model$Xl
  X1 <- model$X1
  X <- model$X
  X2 <- model$X2
  X3 <- model$X3
  X3l <- model$X3l
  X4 <- model$X4
  X4l <- model$X4l
  XS <- model$XS
  XRE <- model$XRE
  nact <- model$nact
  ns <- model$ns
  nre <- model$nre
  nd <- model$nd
  nr <- model$nr
  densVar <- model$densVar
  recVar <- model$recVar
  nrandd <- model$nrandd
  nrandr <- model$nrandr
  netnums <- model$netnums
  nactnets <- model$nactnets
  actnums <- unlist(lapply(1:nnets, function(k) {seq(nactnets[k])}))
  nrows <- model$nrows
  netrows <- unlist(lapply(1:nnets, function(k){rep(k, nrows[k])}))
  if (nrandd>0){drandd <- 1} else {drandd <- 0}
  if (nrandr>0){drandr <- 1} else {drandr <- 0}
  if (densVar == TRUE | nnets == 1){overalld <- 1} else {overalld <- 0}
  if (recVar == TRUE | nnets == 1){overallr <- 1} else {overallr <- 0}
  # number of parameters
  nb <- ns+ nre+ nd
  nrand <- nact*2
  npar <- nrand + nb
  nvarpar <- 2
  nvarcovpar <- nvarpar*(nvarpar+1)/2
  subC1 <- c((nb[1]+1):(nb[1]+nact))
  subC2 <- c((nb[1]+nact+1):npar)
  subC <- c(subC1, subC2)
  subCnets <- c(rep(0, nb[1]), netnums, netnums)
  subCHCnets <- c(netnums, netnums + nnets)
  subCnets2 <- c(rep(0, nb[1]), netnums, netnums+nnets)
  subS1nets <- rep(1:nnets, each=nvarpar^2)
  if ((nnets > 1)){
    subb <- -c(drandd*((ns+nre+1+overalld):(ns+nre+nnets+overalld)),(nb+1):npar)
    subM <- c((ns+nre+overalld+1):(ns+nre+nrandd+overalld))
  } else {
    subb <- c((1):nb)
  }
  if ((nnets > 1)){
    subr <- -c(drandr*((overallr+1):(overallr+nrandr)))
    subR <- c((overallr+1):(overallr+nrandr))
  } else {
    subr <- c(1:nr)
  }
  subm <- c(ns+nre+1)
  subd <- c((ns+nre+1): nd)
  subs <- c(if((ns)>0) 1:(ns))
  subre <- c(if((nre)>0) (1+ns):(ns+nre))
  # prior distributions
  pSDb <- 10
  pSDbsq <- pSDb^2
  sigmab <- c(if (ns>0) {pSDb/apply(X1, 2, sd)}, if (nre>0) {pSDb/apply(X2, 2, sd)}, if (densVar == TRUE | nnets == 1) {pSDb}, if (nrandd>0) {rep(pSDb*sqrt(nnets), nnets)}, if ((nd-(as.numeric(nnets>1))*nnets-overalld)>0) {pSDb/apply(X3[,((as.numeric(nnets>1))*nnets+overalld+1):nd, drop=F], 2, sd)})
  sigmab[is.infinite(sigmab)] <- pSDb
  sigmaR <- 1
  sigmapar <- c(sigmab, rep(sigmaR, nrand))
  pmb <- as.vector(rep(0, npar))
  pVb <- diag(sigmapar^2)
  pmr <- as.vector(rep(0, nr))
  pSDg4 <- 10
  pSDg4sq <- pSDg4^2
  pScaler <- c(if (recVar == TRUE | nnets == 1) {pSDg4}, if (nrandr>0) {rep(pSDg4*sqrt(nnets), nnets)}, if (nr-((as.numeric(nnets>1))*nnets+overallr)>0) {pSDg4/apply(X4[,((as.numeric(nnets>1))*nnets+overallr+1):nr, drop=F], 2, sd)})
  pScaler[is.infinite(pScaler)] <- pSDg4
  pVr <- diag(pScaler^2, nrow=nr)
  pmbr <- c(pmb[c(subb)], pmr[subr])
  pVbr <- diag(c(sigmapar[c(subb)]^2, pScaler[subr]^2))
  invrho <- 1/pSDg4^2
  Xrho <- rep(1, nnets)
  pdfC <- 3
  CPC <- diag(rep(nvarpar))*1
  pVC <- diag(2)
  pdfM <- 2
  CPM <- 0.5
  pVM <- 1
  pdfR <- 2
  CPR <- 0.5
  pVR <- 1
  npC <- 0
  npM <- 0
  npR <- 0
  # create vectors and matrices
  Nsamp <- Nadapt*Sadapt
  beta <- as.vector(rep(0, npar))
  g4 <- as.vector(rep(0, nr))
  bsimsAD <- matrix(rep(NA, npar*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=npar)
  rsimsAD <- matrix(rep(NA, nr*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=nr)
  varCAD <- matrix(rep(NA,nvarpar^2*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=nvarpar^2)
  varMAD <- matrix(rep(NA, 1*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=1)
  sqrtvarMAD <- matrix(rep(NA, 1*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=1)
  varRAD <- matrix(rep(NA, 1*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=1)
  sqrtvarRAD <- matrix(rep(NA, 1*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=1)
  varCADb <- matrix(rep(NA,nvarpar^2*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=nvarpar^2)
  varMADb <- matrix(rep(NA, 1*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=1)
  varRADb <- matrix(rep(NA, 1*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=1)
  varRBI <- matrix(rep(NA, 1*NburnR), nrow= NburnR, ncol=1)
  varC <- matrix(rep(NA,nvarpar^2*Nsamp), nrow= Nsamp, ncol=nvarpar^2)
  varM <- matrix(rep(NA, 1*Nsamp), nrow= Nsamp, ncol=1)
  varR <- matrix(rep(NA, 1*Nsamp), nrow= Nsamp, ncol=1)
  bsims <- matrix(rep(NA, npar*Nsamp), nrow= Nsamp, ncol=npar)
  rsimsBI <- matrix(rep(NA, nr*NburnR), nrow= NburnR, ncol=nr)
  overallrBI <- matrix(rep(NA, 1*NburnR), nrow= NburnR, ncol=1)
  rsims <- matrix(rep(NA, nr*Nsamp), nrow= Nsamp, ncol=nr)
  varS1 <- list()
  varS1W <- list()
  varS1Mat <- list()
  prvarS1 <- list()
  sqrtvarS1 <- list()
  corvarS1 <- list()
  varS1u <- list()
  EvarS1 <- list()
  Ml1l <- list()
  Myl <- list()
  Ml <- list()
  Rl <- list()
  cSignl <- list()
  alphaC <- list()
  alphaM <- list()
  alphaR <- list()
  varc <- list()
  pVCtmp <- list()
  c2Tot <- list()
  rInd <- list() 
  for (k in 1:nnets){
    varS1[[k]] <- pVC
    varS1W[[k]] <- pVC
    varS1Mat[[k]] <- diag(nvarpar*2)
    prvarS1[[k]] <- pVC
    sqrtvarS1[[k]] <- pVC
    corvarS1[[k]] <- 0
    varS1u[[k]] <- pVC
    EvarS1[[k]] <- pVC
    Ml1l[[k]] <- matrix(rep(0, nactnets[k]*nactnets[k]), ncol=nactnets[k])
    Myl[[k]] <- matrix(rep(0, nactnets[k]*nactnets[k]), ncol=nactnets[k])
    Ml[[k]] <- matrix(rep(0, nactnets[k]*nactnets[k]), ncol=nactnets[k])
    Rl[[k]] <- matrix(rep(0, nactnets[k]*nactnets[k]), ncol=nactnets[k])
    cSignl[[k]] <- (2*nets[[k]]-1)*(2*t(nets[[k]])-1)
    alphaC[[k]] <- abs(log(0.1/(nactnets[k] -1 -0.1))) 
    alphaR[[k]] <- abs(2*log(0.25/((nactnets[k] -1)/2 -0.25))) 
    varc[[k]] <- matrix(rep(0, 4), ncol= nvarpar)
    pVCtmp[[k]] <- diag(rep(pSDbsq*nnets*nactnets[k]/2,2))
    c2Tot[[k]] <- diag(nvarpar)*nactnets[k]
    rInd[[k]] <- nets[[k]]*t(nets[[k]])
  }
  alphaM <- alphaC
  varS1tmp <- varS1
  ll1Ct <- as.vector(rep(NA, nnets))
  ll2Ct <- as.vector(rep(NA, nnets))
  lltmp <- as.vector(rep(NA, nnets))
  tmplpC <- as.vector(rep(NA, nnets))
  ll1Mt <- as.vector(rep(NA, nnets))
  ll2Mt <- as.vector(rep(NA, nnets))
  ll1Rt <- as.vector(rep(NA, nnets))
  ll2Rt <- as.vector(rep(NA, nnets))
  lpCgntmp <- as.vector(rep(NA, nact))
  varS1AD <- matrix(rep(NA,(nvarpar^2)*nnets*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=(nvarpar^2)*nnets)
  sqrtvarS1AD <- matrix(rep(NA,(nvarpar^2)*nnets*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=(nvarpar^2)*nnets)
  corvarS1AD <- matrix(rep(NA, nnets*Nadapt*Sadapt), nrow= Nadapt*Sadapt, ncol=nnets)
  betaC2tmp <- matrix(rep(NA, nact*nvarpar), ncol=nvarpar)
  betaRE <- as.vector(rep(0, nnets))
  betaRE2 <- as.vector(rep(0, nnets))
  # calculations only needed once
  tX <- t(X)
  if (ns>0){
    tXS  <- t(XS)
  }
  if (nre>0){
    tXRE  <- t(XRE)
  }
  postdfC <- pdfC + nactnets
  postdfCtot <- pdfC + nact
  postdfM <- pdfM + nnets
  if (nrandd>0) {DM <- diag(nnets)} else {DM <- diag(2)}
  if (nrandr>0) {DR <- diag(nnets)} else {DR <- diag(2)}
  postdfR <- pdfR + nnets
  IpVb <- diag(1/(sigmapar^2), npar)
  net <- blockMatrixDiagonal(nets)
  nett <- t(net)
  cSign <- (2*net-1)*(2*nett-1)
  Dn <- diag(rep(1, nact))
  # Adaptive sequence m and r
  covRWADr <- sqrt(pVr[c(subr), c(subr)])/(nact)
  covRWADb <- (pVb[c(subb), c(subb)])/(nact^2)
  covRWADC <- list()
  covRWADS1 <- list()
  for (k in 1:nnets){
    covRWADC[[k]] <- pVC/nactnets[k]^2
    covRWADS1[[k]] <- diag(1)/nactnets[k]^2
  }
  covRWADMt <- pVM/nact
  covRWADRt <- pVR/nact
  covRWADM <- list()
  covRWADR <- list()
  for (k in 1:nnets){
    covRWADM[[k]] <- pVM/nact
    covRWADR[[k]] <- pVR/nact
  }
  covRWADS2 <- pVR/nact
  covRWADS3 <- pVR/nact
  varS2 <- pVM
  varS2W <- pVM
  EvarS2 <- pVM
  varS3 <- pVR
  varS3W <- pVR
  EvarS3 <- pVR
  VRM <- DM * varS2
  VRR <- DR * varS3
  beta[c(subC)] <- Ccenter(Rfast::rmvnorm(nact, mu= c(0,0), sigma= pVC), netnums)
  if ((nnets > 1)){
    beta[c(subM)] <- rep(c(0), nrandd)
    rwM <- rep(c(0), nrandd)
    m <- beta[c(subM)]
    mtmp <- beta[c(subM)]
    m2Tot <- nnets
    mTot <- m
  } else {
    m <- rep(0, nnets)
    mtmp <- as.vector(rep(0, nnets))
    m2Tot <- 1
  }
  if ((nnets > 1)){
    g4[c(subR)] <-  rep(c(0), nrandr)
    rwR <-  rep(c(0), nrandr)
    r <- g4[c(subR)]
    rtmp <- g4[c(subR)] 
    r2Tot <- nnets
    rTot <- r
  } else {
    r <- rep(0, nnets)
    rtmp <- rep(0, nnets) 
    r2Tot <- 1
  }
  ll1C <- lapply(1:nnets, function(k){llp2ML(yl[[k]], nets[[k]], Xl[[k]], X4l[[k]], c(beta[1:nb], beta[subCnets==k]), g4, Ml[[k]], Myl[[k]], Rl[[k]], rInd[[k]])})
  ll1 <- unlist(lapply(1:nnets, function(k){sum(ll1C[[k]])}))
  c <- cbind(beta[c(subC1)], beta[c(subC2)])
  IC <- matrix(rep(1, nact*2), ncol=2)
  Ibeta <- rep(0, nb)
  Ibeta[subm] <- 1
  ctmp <- c
  mc <- 0
  Sb <- 1
  Sba <- 1
  Sbb <- 1
  Sr <- 1
  Sra <- 1
  Srb <- 1
  SC <- 1/(2*nactnets)
  SvarS1 <- rep(1, nnets)
  SvarS2 <- 1
  SvarS3 <- 1
  SMt <- 1
  SM <- rep(1, nnets)
  SMb <- rep(1, nnets)
  SRt <- 1/nnets
  SR <- rep(1, nnets)
  SRa <- rep(1, nnets)
  SRb <- rep(1, nnets)
  accb <- 0
  sumaccb <- 0
  accr <- 0
  sumaccr <- 0
  accC <- rep(0, nnets)
  accMt <- 0
  accRt <- 0
  accvarS1 <- rep(0, nnets)
  accvarS2 <- 0
  accvarS3 <- 0
  accM <- rep(0, nnets)
  accR <- rep(0, nnets)
  sumaccR <- rep(0, nnets)
  Sbt <- rep(NA, Nadapt)
  Srt <- rep(NA, Nadapt)
  SCt <- matrix(rep(NA, nnets*Nadapt), ncol = Nadapt)
  varm <- 0
  varr <- 0
  MSwM <- 0
  MSwR <- 0
  for (i in 1:Nadapt){
    accb <- 0
    accr <- 0
    accC <- rep(0, nnets)
    accvarS1 <- rep(0, nnets)
    accvarS2 <- 0
    accvarS3 <- 0
    accMt <- 0
    accM <- rep(0, nnets)
    accRt <- 0
    accR <- rep(0, nnets)
    for (j in 1:Sadapt){ 
      num <- ((i-1)*Sadapt + j)
      beta2 <- beta 
      beta2[c(subb)] <- beta[c(subb)] +as.vector(Rfast::rmvnorm(1,pmb[c(subb)], covRWADb))
      g4s <- g4
      g4s[c(subr)] <- g4[c(subr)] + as.vector(Rfast::rmvnorm(1, pmr[c(subr)], covRWADr))
      ll2C <- lapply(1:nnets, function(k){llp2ML(yl[[k]], nets[[k]], Xl[[k]], X4l[[k]], c(beta2[1:nb], beta2[subCnets==k]), g4s, Ml[[k]], Myl[[k]], Rl[[k]], rInd[[k]])})
      ll2 <- unlist(lapply(1:nnets, function(k){sum(ll2C[[k]])}))
      ll1br <- sum(ll1) + Rfast::dmvnorm(t(c(beta[c(subb)], g4[subr])), mu= pmbr, sigma= pVbr, logged = TRUE)
      ll2br <- sum(ll2) + Rfast::dmvnorm(t(c(beta2[c(subb)], g4s[subr])), mu= pmbr, sigma= pVbr, logged = TRUE)
      if (runif(1, min = 0, max = 1) <  min(1, exp(ll2br-ll1br))){
        bsimsAD[((i-1)*Sadapt + j), ] <- beta2
        beta <- beta2
        rsimsAD[((i-1)*Sadapt + j),] <- g4s
        g4 <- g4s
        ll1 <- ll2
        ll1C <- ll2C
        accb <- accb + 1
        accr <- accr + 1
      } else {
        bsimsAD[((i-1)*Sadapt + j), ] <- beta
        rsimsAD[((i-1)*Sadapt + j),] <- g4
      }
      # random actor effects
      for (ccyc in 1:nccyc){
        beta2 <- beta
        cOld <- c
        for (k in 1:nnets){
          betaC2tmp[netnums==k,] <-  Rfast::rmvnorm(nactnets[k], mu= c(0,0), sigma= covRWADC[[k]])
        }
        c2 <- c + betaC2tmp
        beta2[c(subC)] <- beta[c(subC)] + as.vector(betaC2tmp)
        ll2C <- lapply(1:nnets, function(k){llp2ML(yl[[k]], nets[[k]], Xl[[k]], X4l[[k]], c(beta2[1:nb], beta2[subCnets==k]), g4, Ml[[k]], Myl[[k]], Rl[[k]], rInd[[k]])})
        ll1CpC <- unlist(lapply(1:nnets, function(k){sum(ll1C[[k]])})) + lpC(c, varS1, netnums, tmplpC)
        ll2CpC <-  unlist(lapply(1:nnets, function(k){sum(ll2C[[k]])})) + lpC(c2, varS1, netnums, tmplpC) 
        for (k in 1:nnets){
          if (runif(1, min = 0, max = 1) <  min(1, exp(ll2CpC[k]-ll1CpC[k]))){
            ll1C[[k]] <- ll2C[[k]]  
            accC[k] <- accC[k] +1
            c[netnums==k,] <- c2[netnums==k,]
          } 
        }
        bsimsAD[((i-1)*Sadapt + j), c(subC)] <- as.vector(c)
        beta[subC] <- as.vector(c)
        ll1 <- unlist(lapply(1:nnets, function(k){sum(ll1C[[k]])}))
        if (separateSigma == FALSE){
          varS1tot <-  CholWishart::rInvWishart(1, postdfCtot, crossprod(c) + CPC)[,,1]
          for (k in 1:nnets){
            varS1[[k]] <- varS1tot
          }
        } else{
          for (k in 1:nnets){
            varS1[[k]] <- CholWishart::rInvWishart(1, postdfC[k], crossprod(c[netnums==k,]) + CPC)[,,1]
          }
        }
        varS1AD[((i-1)*Sadapt + j),] <- unlist(varS1)
        varCAD[((i-1)*Sadapt + j),] <- apply(simplify2array(varS1), 1:2, mean) 
      }
      # random effects M and R
      if ((nnets > 1)){
        #
        beta2 <- beta
        for (k in 1:nnets){
          rwM[k] <-  Rfast::rmvnorm(1, mu= c(0), sigma= covRWADM[[k]])
        }
        m2 <-  beta[subM] + rwM
        beta2[subM] <- m2
        ll2C <- lapply(1:nnets, function(k){llp2ML(yl[[k]], nets[[k]], Xl[[k]], X4l[[k]], c(beta2[1:nb], beta[subCnets==k]), g4, Ml[[k]], Myl[[k]], Rl[[k]], rInd[[k]])})
        ll2 <- unlist(lapply(1:nnets, function(k){sum(ll2C[[k]])}))
        if (densVar == TRUE) {
          ll1M <- ll1 +  dnorm(m, mean= c(0), sd= sqrt(varS2), log = TRUE)
          ll2M <- ll2 + dnorm(m2, mean= c(0), sd= sqrt(varS2), log = TRUE)
        } else {
          ll1M <- ll1 + dnorm(m, mean= c(0), sd= pSDb, log = TRUE)
          ll2M <- ll2 + dnorm(m2, mean= c(0), sd= pSDb, log = TRUE)
        }
        for (k in 1:nnets) {
          if (runif(1, min = 0, max = 1) <  min(1, exp(ll2M[k]-ll1M[k]))){
            ll1C[[k]] <- ll2C[[k]]
            accM[k] <- accM[k] +1
            m[k] <- m2[k]
          }
        }
        #
        ll1 <- unlist(lapply(1:nnets, function(k){sum(ll1C[[k]])}))
        mMat <-  cbind(m)
        beta[subM] <- m
        bsimsAD[((i-1)*Sadapt + j), ] <- beta
        if (densVar == TRUE) {
          varS2 <- CholWishart::rInvWishart(1, postdfM, crossprod(mMat) + CPM)[,,1]
          VRM <- DM * varS2
          varMAD[((i-1)*Sadapt + j),] <- varS2
        }
        #
      }
      if ((nnets > 1)){
        rOld <- r
        g4s <- g4
        for (k in 1:nnets){
          rwR[k] <-  Rfast::rmvnorm(1, mu= c(0), sigma= covRWADR[[k]])
        }
        r2 <-  g4[subR] + rwR
        g4s[subR] <- r2
        ll2C <- lapply(1:nnets, function(k){llp2ML(yl[[k]], nets[[k]], Xl[[k]], X4l[[k]], c(beta[1:nb], beta[subCnets==k]), g4s, Ml[[k]], Myl[[k]], Rl[[k]], rInd[[k]])})
        ll2 <- unlist(lapply(1:nnets, function(k){sum(ll2C[[k]])}))
        if (recVar == TRUE) { 
          ll1R <- ll1 +  dnorm(r, mean= c(0), sd= sqrt(varS3), log = TRUE) 
          ll2R <- ll2 + dnorm(r2, mean= c(0), sd= sqrt(varS3), log = TRUE) 
        } else {
          ll1R <- ll1 + dnorm(r, mean= c(0), sd= pSDg4, log = TRUE)
          ll2R <- ll2 + dnorm(r2, mean= c(0), sd= pSDg4, log = TRUE)
        }
        for (k in 1:nnets) {
          if (runif(1, min = 0, max = 1) <  min(1, exp(ll2R[k]-ll1R[k]))){
            ll1C[[k]] <- ll2C[[k]]  
            accR[k] <- accR[k] +1
            r[k] <- r2[k]
          } 
        }
        g4[subR] <- r
        ll1 <- unlist(lapply(1:nnets, function(k){sum(ll1C[[k]])}))
        rMat <- cbind(r) 
        rsimsAD[((i-1)*Sadapt + j), ] <- g4
        rsimsAD[((i-1)*Sadapt + j), subR] <- r
        if (recVar == TRUE){
          varS3 <- CholWishart::rInvWishart(1, postdfR, crossprod(rMat) + CPR)[,,1]
          VRR <- DR * varS3
          varRAD[((i-1)*Sadapt + j),] <- varS3
        }
      }
      callback()
    }
    sumSadapt <- Sadapt*i
    sumgacc <- gacc*i
    sumaccb <- sumaccb + accb
    fc <- 1/sqrt(i)
    if (sumaccb > sumgacc){
      Sba <- Sba*(1+(1-(sumSadapt-sumaccb)/(sumSadapt-sumgacc)))  
    } else {
      Sba <- Sba/(1+(1-(sumaccb/sumgacc)))
    }
    if (accb > gacc){
      Sbb <- Sbb*(1+fc*(1-(Sadapt-accb)/(Sadapt-gacc)))  
    } else {
      Sbb <- Sbb/(1+fc*(1-(accb/gacc)))
    }
    Sb <- Sbb
    Sbt[i] <- Sb
    covRWADb <- Sb*cov(as.matrix(bsimsAD[,subb]), use= "complete.obs")
    if (length(eigen(covRWADb)$values[eigen(covRWADb)$values >1*10^-15]) < length(subb)){diag(covRWADb) <- diag(covRWADb) + 1*10^-10}
    sumaccr <- sumaccr + accr
    if (sumaccr > sumgacc){
      Sra <- Sra*(1+(1-(sumSadapt-sumaccr)/(sumSadapt-sumgacc)))  
    } else {
      Sra <- Sra/(1+(1-(sumaccr/sumgacc)))
    }
    if (accr > gacc){
      Srb <- Srb*(1+fc*(1-(Sadapt-accr)/(Sadapt-gacc)))  
    } else {
      Srb <- Srb/(1+fc*(1-(accr/gacc)))
    }
    Sr <- Srb
    Srt[i] <- Sr
    covRWADr <- Sr*cov(as.matrix(rsimsAD[,subr]), use= "complete.obs")
    if (length(eigen(covRWADr)$values[eigen(covRWADr)$values >1*10^-15]) < length(subr)){diag(covRWADr) <- diag(covRWADr) + 1*10^-8}
    #
    for (k in 1:nnets){
      if (accC[k] > gaccC){
        SC[k] <- SC[k]*(1+fc*(1-(Sadapt-accC[k])/(Sadapt-gaccC)))  
      } else {
        SC[k] <- SC[k]/(1+fc*(1-(accC[k]/gaccC)))
      }
      covRWADC[[k]] <- SC[k]*matrix(colMeans(varS1AD[,subS1nets==k], na.rm = T), ncol=2)
      if (length(eigen(covRWADC[[k]])$values[eigen(covRWADC[[k]])$values >1*10^-15]) < nvarpar){diag(covRWADC[[k]]) <- diag(covRWADC[[k]]) + 1*10^-10}
      SCt[k,i] <-  SC[k] 
    }
    if ((nnets > 1)){
      for (k in 1:nnets){
        if (accM[k] > gaccM){
          SMb[k] <- SMb[k]*(1+fc*(1-(Sadapt-accM[k])/(Sadapt-gaccM)))
        } else {
          SMb[k] <- SMb[k]/(1+fc*(1-(accM[k]/gaccM)))
        }
        SM[k] <- SMb[k]
        if (densVar == TRUE){
          covRWADM[[k]] <- SM[k]*matrix(colMeans(varMAD, na.rm = T), ncol=1)
        } else {
          covRWADM[[k]] <- SM[k]*cov(as.matrix(bsimsAD[,subM[k]]), use= "complete.obs")
          if (length(eigen(covRWADM[[k]])$values[eigen(covRWADM[[k]])$values >1*10^-15]) < 1){covRWADM[[k]] <- 1*10^-8}
        }
      }
      for (k in 1:nnets){
        if (accR[k] > gaccR){
          SRb[k] <- SRb[k]*(1+fc*(1-(Sadapt-accR[k])/(Sadapt-gaccR)))
        } else {
          SRb[k] <- SRb[k]/(1+fc*(1-(accR[k]/gaccR)))
        }
        SR[k] <- SRb[k]
        if (recVar == TRUE){
          covRWADR[[k]] <- SR[k]*matrix(colMeans(varRAD, na.rm = T), ncol=1)
        } else {
          covRWADR[[k]] <- SR[k]*cov(as.matrix(rsimsAD[,subR[k]]), use= "complete.obs")
          if (length(eigen(covRWADR[[k]])$values[eigen(covRWADR[[k]])$values >1*10^-15]) < 1){covRWADR[[k]] <- 1*10^-8}
        }
      }
    }
  } 
  MCMCsims <- cbind(varCAD[burn,c(1,2,4)], if (densVar == TRUE) varMAD[burn,1], if (recVar == TRUE) varRAD[burn,1], bsimsAD[burn,1:nb], rsimsAD[burn,1:nr]) 
  colnames(MCMCsims) <- c("sender variance", "sender receiver covariance","receiver variance", if (densVar == TRUE & nnets > 1) c("density variance"), 
                          if (recVar == TRUE & nnets > 1) c("reciprocity variance"), all.vars(sender), all.vars(receiver), 
                          if (overalld == 1) {"density"}, if (nnets > 1) { if (overalld == 0) {sprintf("density net %d",seq(1:nnets))} else { if (densVar == TRUE) sprintf("net %d",seq(1:nnets))}}, 
                          all.vars(density), if (overallr == 1) {"reciprocity"}, 
                          if (nnets > 1) { if (overallr == 0) {sprintf("reciprocity net %d",seq(1:nnets))} else { if (recVar == TRUE) sprintf("net %d",seq(1:nnets))}}, 
                          all.vars(reciprocity))
  z <- list(MCMCsims = MCMCsims, y = net, densVar=densVar, recVar=recVar, drandd=drandd, drandr=drandr, nrandd=nrandd, nrandr=nrandr, ns=ns, nre=nre, nd=nd-overalld-nrandd, nr=nr-overallr-nrandr, acc=c(accvarS1,accvarS2,accvarS3), nnets=nnets)
  class(z) <- c("2MLb")
  return(z)
}


