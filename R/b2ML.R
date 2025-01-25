b2ML <- function (nets, actor = NULL, density = NULL, adapt = NULL, burnin = NULL, center = NULL, separate= NULL, densVar = NULL, seed = NULL) 
{
  # sampling parameters
  if(!is.null(adapt)){
    Nadapt <- adapt
  } else {
    Nadapt <- 100
  } 
  nccyc <- 1
  Sadapt <- 125
  gacc <- round(Sadapt*.234)
  gaccC <- round(Sadapt*.234)
  gaccMt <- round(Sadapt*.234)
  if(!is.null(burnin)){
    NburnR <- burnin
  } else {
    NburnR <- 5000
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
  if(!is.null(separate)){
    separate <- separate
  } else {
    separate <- FALSE
  } 
  if(!is.null(densVar)){
    densVar <- densVar
  } else {
    densVar <- TRUE
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
  model <- b2modelB2ML1l(nets, actor, density, center, separate, densVar)
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
  nrandd <- model$nrandd
  nrandr <- model$nrandr
  netnums <- model$netnums
  nactnets <- model$nactnets
  actnums <- unlist(lapply(1:nnets, function(k) {seq(nactnets[k])}))
  nrows <- model$nrows
  netrows <- unlist(lapply(1:nnets, function(k){rep(k, nrows[k])}))
  if (nrandd>0){drandd <- 1} else {drandd <- 0}
  if (nrandr>0){drandr <- 1} else {drandr <- 0}
  if (separate == FALSE){overalld <- 1; overallr <- 1} else {overalld <- 0; overallr <- 0}
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
  if ((separate == FALSE) & (nnets > 1) & (densVar == TRUE)){
    subb <- -c(if ((ns+nre)>0) 1:(ns+nre), drandd*((ns+nre+1):(ns+nre+nnets+1)),(nb+1):npar)
    if((ns+nre)==0){
      subb <- -c(if ((ns+nre)>0) 1:(ns+nre), drandd*((ns+nre+1):(ns+nre+nnets+1)),(nb+1):npar)
    }
    subr <- -c(drandr*((overallr):(overallr+nrandr)))
    if((nr-nrandr)==1){
      subr <- -c(drandr*((overallr+1):(overallr+nrandr)))
    }
    subm <- c(ns+nre+1)
    subd <- c((ns+nre+1): nd)
    subM <- c((ns+nre+2):(ns+nre+nrandd+1))
    subR <- c((overallr+1):(overallr+nrandr))
    subs <- c(if((ns)>0) 1:(ns))
    subre <- c(if((nre)>0) (1+ns):(ns+nre))
  } else {
    subb <- c((ns+nre+2):nb)
    subr <- c(1:nr)
    subm <- c(ns+nre+1)
    subd <- c((ns+nre+1): nd)
    subM <- NA
    subs <- c(if((ns)>0) 1:(ns))
    subre <- c(if((nre)>0) (1+ns):(ns+nre))
  }
  # prior distributions
  pSDb <- 3
  pSDbsq <- pSDb^2
  sigmab <- c(if (ns>0) {pSDb/apply(X1, 2, sd)}, if (nre>0) {pSDb/apply(X2, 2, sd)}, if (separate == FALSE) {pSDb}, if (nrandd>0) {rep(pSDb*sqrt(nnets), nnets)}, if ((nd-(as.numeric(nnets>1))*nnets-overalld)>0) {pSDb/apply(X3[,((as.numeric(nnets>1))*nnets+overalld+1):nd, drop=F], 2, sd)})
  sigmab[is.infinite(sigmab)] <- pSDb
  sigmaR <- 1
  sigmapar <- c(sigmab, rep(sigmaR, nrand))
  pmb <- as.vector(rep(0, npar))
  pVb <- diag(sigmapar^2)
  pmr <- as.vector(rep(0, nr))
  invmu <- 1/pSDb^2
  invCmu <- 1/((0.5*pSDb)^2)
  Xmu <- rep(1, nnets)
  XCmu <- rep(1, nact)
  if ((separate == FALSE) & (nnets > 1) & (densVar == TRUE)){
    XCAREmu <- cbind(rep(1, nact),matrix(unlist(lapply(1:(nnets), function(k){as.numeric(netnums==k)})), ncol=nnets))
  } else {
    XCAREmu <- rep(1, nact)
  }
  if (nre > 0){
    XCAREmu <- cbind(XCAREmu, XRE)
  }
  varHCmu <- diag(nact*2)
  varHCAmu <- diag(nact)
  XCm <- matrix(unlist(lapply(1:nnets, function(k){as.numeric(netnums==k)})), ncol=nnets)
  invCovRe <- list()
  if (nre > 0){
    for (k in 1:nnets){
      invCovRe[[k]] <- diag((1/(pSDb/apply(as.matrix(XRE[netnums==k,], ncol= nre), 2, sd))^2), ncol= nre)
      for (l in 1:nre){
        if (is.infinite(invCovRe[[k]][l,l])){invCovRe[[k]][l,l] <- pSDb}
      }
    }
  }
  if (nre>0){
    invCovReTot <- diag((1/(pSDb/apply(as.matrix(XRE, ncol= nre), 2, sd))^2), ncol= nre)
  }
  invCovmM <- diag(1/(sigmab[(ns+nre+1):(ns+nre+overalld+nrandd)]))
  pmbr <- c(pmb[c(subb)])
  pVbr <- diag(as.numeric(sigmapar[c(subb)]^2), nrow = length(sigmapar[c(subb)]))
  Xrho <- rep(1, nnets)
  pdfC <- 3
  CPC <- diag(rep(nvarpar))*1
  pVC <- diag(2)
  pdfM <- 2
  CPM <- 0.5
  pVM <- 1
  varS1Mat <- CPC
  varS2Mat <- CPM
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
  if (nre>0){
    tXRE  <- t(XRE)
  }
  postdfC <- pdfC + nactnets
  postdfCtot <- pdfC + nact
  postdfM <- pdfM + nnets
  if (nrandd>0) {DM <- diag(nnets)} else {DM <- diag(2)}
  if (nrandr>0) {DR <- diag(nnets)} else {DR <- diag(2)}
  IpVb <- diag(1/(sigmapar^2), npar)
  net <- blockMatrixDiagonal(nets)
  nett <- t(net)
  cSign <- (2*net-1)*(2*nett-1)
  Dn <- diag(rep(1, nact))
  # Adaptive sequence m and r
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
  # intialiseer random effecten op waarden ongelijk aan nul
  beta[c(subC)] <- rep(Ccenter(Rfast::rmvnorm(nact, mu= c(0), sigma= pVC[1,1]), netnums), 2)
  if ((separate == FALSE) & (nnets > 1) & (densVar == TRUE)){
    beta[c(subM)] <- scale(as.vector(Rfast::rmvnorm(nrandd, c(0), pVM)), center=TRUE, scale=FALSE)
    g4[c(subR)] <- scale(as.vector(Rfast::rmvnorm(nrandr, c(0), pVR)), center=TRUE, scale=FALSE) 
    m <- beta[c(subM)]
    mtmp <- beta[c(subM)]
    r <- g4[c(subR)]
    rtmp <- g4[c(subR)] 
    m2Tot <- nnets
    r2Tot <- nnets
    mTot <- m
    rTot <- r
  } else {
    m <- rep(0, nnets)
    mtmp <- as.vector(rep(0, nnets))
    r <- rep(0, nnets)
    rtmp <- rep(0, nnets) 
    m2Tot <- 1
    r2Tot <- 1
  }
  ll1C <- lapply(1:nnets, function(k){llb2MLC(yl[[k]], nets[[k]], Xl[[k]], c(beta[1:nb], beta[subCnets==k]), Ml[[k]], Myl[[k]])})
  ll1 <- unlist(lapply(1:nnets, function(k){sum(ll1C[[k]])/2}))
  c <- cbind(beta[c(subC1)], beta[c(subC2)])
  IC <- matrix(rep(1, nact*2), ncol=2)
  Ibeta <- rep(0, nb)
  Ibeta[subm] <- 1
  ctmp <- c
  cTot <- c
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
    accRt <- 0
    accR <- rep(0, nnets)
    for (j in 1:Sadapt){ 
      num <- ((i-1)*Sadapt + j)
      if (length(beta[c(subb)]) > 0){
        beta2 <- beta 
        beta2[c(subb)] <- beta[c(subb)] +as.vector(Rfast::rmvnorm(1,pmb[c(subb)], covRWADb))
        ll2C <- lapply(1:nnets, function(k){llb2MLC(yl[[k]], nets[[k]], Xl[[k]], c(beta2[1:nb], beta2[subCnets==k]), Ml[[k]], Myl[[k]])})
        ll2 <- unlist(lapply(1:nnets, function(k){sum(ll2C[[k]])/2}))
        ll1br <- sum(ll1) + Rfast::dmvt(t(c(beta[c(subb)])), mu= pmbr, sigma= pVbr, nu=7, logged = TRUE)
        ll2br <- sum(ll2) + Rfast::dmvt(t(c(beta2[c(subb)])), mu= pmbr, sigma= pVbr, nu=7, logged = TRUE)
        if (runif(1, min = 0, max = 1) <  min(1, exp(ll2br-ll1br))){
          bsimsAD[((i-1)*Sadapt + j), ] <- beta2
          beta <- beta2
          ll1 <- ll2
          ll1C <- ll2C
          accb <- accb + 1
          accr <- accr + 1
        } else {
          bsimsAD[((i-1)*Sadapt + j), ] <- beta
        }
      }
      # random actor effects
      for (ccyc in 1:nccyc){
        beta2 <- beta
        for (k in 1:nnets){
          betaC2tmp[netnums==k,] <-  rep(Rfast::rmvnorm(nactnets[k], mu= c(0), sigma= covRWADC[[k]][1,1]), 2)
        }
        c2 <- c + betaC2tmp
        beta2[c(subC)] <- beta[c(subC)] + as.vector(betaC2tmp)
        ll2C <- lapply(1:nnets, function(k){llb2MLC(yl[[k]], nets[[k]], Xl[[k]], c(beta2[1:nb], beta2[subCnets==k]),  Ml[[k]], Myl[[k]])})
        ll1CpC <- unlist(lapply(1:nnets, function(k){sum(ll1C[[k]])/2})) + lpC(c, varS1, netnums, tmplpC)/2 #+ lpCm2(ctmp, alphaC, netnums, tmplpC) 
        ll2CpC <- unlist(lapply(1:nnets, function(k){sum(ll2C[[k]])/2})) + lpC(c2, varS1, netnums, tmplpC)/2 #+ lpCm2(c2tmp, alphaC, netnums, tmplpC) 
        for (k in 1:nnets){
          if (runif(1, min = 0, max = 1) <  min(1, exp(ll2CpC[k]-ll1CpC[k]))){
            ll1C[[k]] <- ll2C[[k]]  
            accC[k] <- accC[k] +1
            c[netnums==k,] <- c2[netnums==k,]
          } 
        }
        bsimsAD[((i-1)*Sadapt + j), c(subC)] <- as.vector(c)
        beta[subC] <- as.vector(c)
        ll1 <- unlist(lapply(1:nnets, function(k){sum(ll1C[[k]])/2}))
        varS1tot <-  diag(rep(CholWishart::rInvWishart(1, postdfCtot, diag(diag(crossprod(c))) + CPC)[,,1][1,1], 2))
        for (k in 1:nnets){
          varS1[[k]] <- varS1tot
        }
        varS1AD[((i-1)*Sadapt + j),] <- unlist(varS1)
        varCAD[((i-1)*Sadapt + j),] <- apply(simplify2array(varS1), 1:2, mean)
      }
      if ((separate == FALSE) & (nnets > 1) & (densVar == TRUE)){
        # draw M
        invVarS1 <- t(sapply(lapply(1:nnets, function(k){chol2inv(chol(varS1[[k]]))}), '['))
        tXCAREmuvarHCmu <- t(XCAREmu*(c(rep(invVarS1[,1], nactnets)))) 
        if (nre>0){
          covHCAmu <- chol2inv(chol(tXCAREmuvarHCmu%*%XCAREmu + diag(c(invCmu, rep(4/varS2, nnets), diag(invCovReTot)))))
        } else {
          covHCAmu <- chol2inv(chol(tXCAREmuvarHCmu%*%XCAREmu + diag(c(invCmu, rep(4/varS2, nnets)))))
        }
        ustarA <-  as.vector(c[,2])  + XCAREmu%*%c(beta[subm]/2, m/2, beta[(subre)])
        mHCAmu <- covHCAmu%*%(tXCAREmuvarHCmu%*%ustarA)
        munewA <- as.vector(Rfast::rmvnorm(1,mHCAmu, covHCAmu))
        cA <- ustarA - XCAREmu%*%munewA
        c <- matrix(rep(cA,2), ncol=2)
        beta[subC]  <- as.vector(c)
        beta[subm] <-  munewA[c(1)]*2
        beta[subM] <-  munewA[2:(nnets+1)]*2
        m <- beta[subM]
        if (nre>0){
          beta[subs] <- munewA[(1+nnets+1):(1+nnets+ns)]
          beta[subre] <- beta[subs]
        }
        bsimsAD[((i-1)*Sadapt + j), ] <- beta
      }   
      if ((nnets == 1) | (densVar == FALSE)){
        k <- 1
        # draw m
        invVarS1 <- t(sapply(lapply(1:nnets, function(k){chol2inv(chol(varS1[[k]]))}), '['))
        tXCAREmuvarHCmu <- t(XCAREmu*(c(rep(invVarS1[,1], nactnets)))) 
        if (nre>0){
          covHCAmu <- chol2inv(chol(tXCAREmuvarHCmu%*%XCAREmu + diag(c(invCmu, diag(invCovReTot)))))
        } else {
          covHCAmu <- chol2inv(chol(tXCAREmuvarHCmu%*%XCAREmu + invCmu))
        }
        ustarA <-  as.vector(c[,2])  + XCAREmu%*%c(beta[subm]/2, beta[(subre)])
        mHCAmu <- covHCAmu%*%(tXCAREmuvarHCmu%*%ustarA)
        munewA <- as.vector(Rfast::rmvnorm(1,mHCAmu, covHCAmu))
        cA <- ustarA - XCAREmu%*%munewA
        c <- matrix(rep(cA,2), ncol=2)
        beta[subC]  <- as.vector(c)
        beta[subm] <-  munewA[c(1)]*2
        if (nre>0){
          beta[subs] <- munewA[2:(1+ns)]
          beta[subre] <- beta[subs]
        }
        bsimsAD[((i-1)*Sadapt + j), ] <- beta
      }  
      # random effects M 
      if ((separate == FALSE) & (nnets > 1) & (densVar == TRUE)){
        mMat <- m
        varS2 <- CholWishart::rInvWishart(1, postdfM, crossprod(mMat) + CPM)[,,1]
        VRM <- DM * varS2
        varMAD[((i-1)*Sadapt + j),] <- varS2 
      }
      callback()
    }
    sumSadapt <- Sadapt*i
    sumgacc <- gacc*i
    sumaccb <- sumaccb + accb
    fc <- 1/sqrt(i)
    if (length(beta[c(subb)]) > 0){
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
    }
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
    if ((separate == FALSE) & (nnets > 1) & (densVar == TRUE)){
      if (accMt > gaccMt){
        SMt <- SMt*(1+fc*(1-(Sadapt-accMt)/(Sadapt-gaccMt)))  
      } else {
        SMt <- SMt/(1+fc*(1-(accMt/gaccMt)))
      }
      covRWADMt <- SMt*matrix(colMeans(varMAD, na.rm = T), ncol=1)
    }
  }
  selFpar  <- 1:nb 
  MCMCsims <- cbind(varCAD[burn,c(1)], if (drandd>0) varMAD[burn,1], bsimsAD[burn,1:nb]) 
  colnames(MCMCsims) <- c("actor variance", if (separate == FALSE & nnets > 1 & densVar == TRUE) c("density variance"), 
                          all.vars(actor), all.vars(actor),
                          if (separate == FALSE) {"density"}, if (nnets > 1) { if (separate == TRUE) {sprintf("density net %d",seq(1:nnets))} else { if (densVar == TRUE) sprintf("net %d",seq(1:nnets))}}, 
                          all.vars(density))
  z <- list(MCMCsims = MCMCsims, y = net, separate= separate, drandd=drandd, nrandd=nrandd, ns=ns, nre=nre, nd=nd-overalld-nrandd, acc=c(accvarS1,accvarS2,accvarS3))
  class(z) <- c("b2ML")
  return(z)
}