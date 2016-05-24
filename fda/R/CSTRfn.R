CSTRfn <- function(parvec, datstruct, fitstruct,
                   CSTRbasis, lambda, gradwrd=TRUE){
#function [res, Dres, fitstruct, df, gcv] =
#    CSTRfn(parvec, datstruct, fitstruct, CSTRbasis, lambda, gradwrd)

# Last modified with R port 2007.07.21 by Spencer Graves
#%  Previously modified 29 July 2005

  max.log.betaCC <- (log(.Machine$double.xmax)/3)
# For certain values of 'coef',
# naive computation of betaCC will return +/-Inf,
# which generates NAs in Dres.
# Avoid this by clipping betaCC
#
# log(.Machine$double.xmax)/2 is too big,
# because a multiple of it is squared in CSTRfn ...
#

#if nargin < 6, gradwrd = 1;  end
##
## 1.  Modify fitstruct for use with CSTRfitLS
##
#  cat("CSTRfn: parvec =", parvec, "\n")
#
  eps <- .Machine$double.eps
  eps1 <- (1+2*eps)
  fit <- fitstruct$fit
  if(is.null(fit))
    stop("fitstruct has no 'fit' component")

#%  load parameters
#[kref, EoverR, a, b] = par2vals(parvec, fitstruct)
#%  set up fitstruct
#fitstruct.kref   = kref
#fitstruct.EoverR = EoverR
#fitstruct.a      = a
#fitstruct.b      = b

  estimate <- fitstruct$estimate
  if(is.null(estimate))
    stop("fitstruct has no 'estimate' component")
#  Compare estimate with parvec
  if(sum(estimate) != length(parvec)){
    cat("ERROR in CSTRfn:  sum(estimate) != length(parvec)\n")
    cat("parvec = ", parvec, "; \n")
    stop("fitstruct$estimate = ",
         paste(estimate, collapse=" "))
  }
  m <- 0
  {
#  1.1.  Estimate kref starting from parvec or use fitstruct?
    if(estimate[1]){
      m <- m+1
      kref <- parvec[m]
      fitstruct$kref <- kref
    }
    else
      kref <- fitstruct$kref
  }
  {
#  1.2.  Estimate EoverR starting from parvec or use fitstruct?
    if(estimate[2]){
      m <- m+1
      EoverR <- parvec[m]
      fitstruct$EoverR <- EoverR
    }
    else
      EoverR<- fitstruct$EoverR
  }
  {
#  1.3.  Estimate 'a' starting from parvec or use fitstruct?
    if(estimate[3]){
      m <- m+1
      a <- parvec[m]
      fitstruct$a <- a
    }
    else
      a <- fitstruct$a
  }
  {
#  1.4.  Estimate b starting from parvec or use fitstruct?
    if(estimate[4]){
      m <- m+1
      b <- parvec[m]
      fitstruct$b <- b
    }
    else
      b<- fitstruct$b
  }
##
## 2.  Set up inner optimization:  optimize fit with respect to coef
##
#  2.1.  Set up
  tolval <- 1e-10
  itermax <- 10
  coef0 <- as.matrix(fitstruct$coef0)
  dim.coef0 <- dim(coef0)
  if(length(dim.coef0)>1){
    coef0 <- as.vector(coef0)
    if(length(dim.coef0)==2 && (dim.coef0[2]==2)){
      coefNames <- outer(c("Conc", "Temp"), 1:dim.coef0[1], paste,
                         sep="")
      names(coef0) <- t(coefNames)
    }
  }
  ncoef <- length(coef0)
#  2.2.  initial value
#[res0, Dres] = CSTRfitLS(coef0, datstruct, fitstruct, lambda, 1)
  CSTR0 <- CSTRfitLS(coef0, datstruct, fitstruct, lambda, 1)

#  2.3.  Reshape for easier manipulation
  res0 <- with(CSTR0$res, c(Sres, Lres))
#
  N <- dim(CSTR0$res$Sres)[1]
  k12 <- dim(CSTR0$res$Sres)[2]
  nquad <- dim(CSTR0$res$Lres)[1]
  n.Sres <- N*k12
  n.Lres <- nquad*2
  Nval <- n.Sres + n.Lres
#
  nbasis <- (dim(CSTR0$Dres$DSres)[2]/2)
  Dres0 <- with(CSTR0$Dres, rbind(DSres, DLres))
# res0.mat <- readMat("CSTRfitLSres0.mat")# OK
# d.res0 <- (res0-res0.mat$res0)
# quantile(d.res0)
#      0%      25%      50%      75%     100%
#-7.1e-06 -1.4e-07  1.3e-08  2.1e-07  6.0e-06
# res0 good.

# Dres.mat <- read.csv("CSTRfitLSDres.csv", header=FALSE)
# d.Dres <- (Dres0-as.matrix(Dres.mat))
# quantile(d.Dres)
#           0%           25%           50%           75%          100%
#-0.0006900152  0.0000000000  0.0000000000  0.0000000000  0.0006814540
# sqrt(mean(Dres0^2)) # 0.747
# sqrt(mean(d.Dres^2)) # 9.8e-6
# Good.

# F0 = mean(res0.^2)
  F0 <- mean(res0^2)
  F1 <- 0
  iter <- 0
  fundif <- 1
# gradnorm0 = mean(res0'*Dres)
  r.D <- crossprod(res0, Dres0)
  gradnorm0 <- mean(r.D)
  if(is.na(gradnorm0))
    stop("Initial call to CSTRfitLS returned NAs with parvec = ",
         paste(parvec, collapse=", "), ";  sum(is.na(res0)) = ",
         sum(is.na(res0)), "; sum(is.na(Dres0)) = ", sum(is.na(Dres0)))
#
  gradnorm1 <- 1
##
## 3.  Gauss-Newton optimization loop:
##     optimize fit with respect to coef
##
  while(( gradnorm1 > tolval) | (fundif > tolval)){
#
    iter <- iter + 1
    if( iter > itermax) break
#
#    Dcoef = Dres\res0
    nNA.res0 <- sum(is.na(res0))
    nNA.Dres0 <- sum(is.na(Dres0))
    if(nNA.res0 | nNA.Dres0) {
      dump("parvec", "parvecError.R")
      cat("Error: ", nNA.res0, "and", nNA.Dres0,
          "NAs found in res0 and Dres0 with parvec =",
          parvec, "; parvec dumped to parvecError.R\n")
    }
#
#    Dcoef <- (lm.fit(Dres0, res0)$coefficients)
#
    Dres.svd <- svd2(Dres0)
    ikeep <- with(Dres.svd, which(d > eps*max(d)))
    Dres.rank <- length(ikeep)
    if(Dres.rank < min(dim(Dres0)))
      warning("Dres0 has rank ", Dres.rank, " < dim(Dres0) = ",
              paste(dim(Dres0), collapse=", "),
              " in iteration ", iter,
              ".  Ignoring singular subspace.  ",
              "First (rank+1) singular values = ",
              paste(Dres.svd$d[1:(Dres.rank+1)], collapse=", "))
    Dcoef <- with(Dres.svd, v %*% ((t(u) %*% res0) / d))
#
#    %  initial step:  alpha = 1
    coef1 <- coef0 - Dcoef
#
#[res1, Dres] = CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
    CSTR1 <- CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
    res1 <- with(CSTR1$res, c(Sres, Lres))
    Dres1 <- with(CSTR1$Dres, rbind(DSres, DLres))
#
    F1 <- mean(res1^2)
    alpha <- 1
    fundif <- abs(F0-F1)/abs(F0)

#%     %  smaller steps as required, halving alpha each time
#%     while F1 >= F0*(1+2*eps)
    while(is.na(F1) || F1>= (eps1*F0) || any(is.na(Dres1))){
      alpha <- alpha/2
      if(is.na(F1)){
        n.na <- sum(is.na(res1))
        attr(coef1, "n.na.in.res1") <- n.na
#        ..CSTRfn.coef1.gen.NA <<- list(parvec=parvec, coef1=coef1)
        warning(n.na, " NAs in res1.")
#        warning(n.na, " NAs in res1;  coef1 saved in ",
#                "'..CSTRfn.coef1.gen.NA'")
      }
      if(alpha< 1e-4){
        pv <- paste(signif(parvec, 3), collapse=", ")
#        ..CSTRfn.coef1.gen.5 <<- list(parvec=parvec, coef1=coef1)
#        warning("Stepsize reduced below the minimum with parvec = ",
#             pv, " on iteration ", iter,
#             " in trying to optimize ", length(coef0),
#             " coefficients;  using suboptimal coefficients;  ",
#                "saved in '..CSTRfn.coef1.gen.5'")
        warning("Stepsize reduced below the minimum with parvec = ",
             pv, " on iteration ", iter,
             " in trying to optimize ", length(coef0),
             " coefficients;  using suboptimal coefficients. ")
        break
      }
#
      coef1 <- coef0 - alpha*Dcoef
#%    [res1, Dres] = CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
      CSTR1 <- CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
      res1 <- with(CSTR1$res, c(Sres, Lres))
      Dres1 <- with(CSTR1$Dres, rbind(DSres, DLres))
      F1 <- mean(res1^2)
      fundif <- abs(F0-F1)/abs(F0)
    }
#    gradnorm1 = mean(res1'*Dres)
    gradnorm1 <- mean(crossprod(res1, Dres1))
#%     disp(num2str([iter, F1, fundif, gradnorm1]))
    coef0 <- coef1
    res0 <- res1
    Dres0 <- Dres1
    F0 <- F1
    gradnorm0 <- gradnorm1
# end while(( gradnorm1 > tolval) | (fundif > tolval)){
  }
##
##% 4.  update coef
##
  coef <- coef0
  fitstruct$coef0 <- coef
#
# DresMat <- read.csv("CSTRfnDres.csv", header=FALSE)
# d.Dres <- Dres - as.matrix(DresMat)
# quantile(d.Dres)
#           0%           25%           50%           75%          100%
#-0.0006207360  0.0000000000  0.0000000000  0.0000000000  0.0006317204
# sqrt(mean(Dres^2)) # 0.747
# sqrt(mean(d.Dres^2))# 1.0e-5
#
##
##% 5.  compute df and gcv
##
  Zmat <- Dres0[1:ncoef,]
#  Nval = length(res1)
  Rfac <- Dres0[(ncoef+1):Nval,]
#  Smat = Zmat*inv(Zmat'*Zmat + Rfac'*Rfac)*Zmat'
# Use singular value decomposition so we never have to worry about
# ill conditioning.
  Zsvd <- svd2(Zmat)
# Zmat = with(Zsvd, u %*% diag(d) %*% t(v))
# so Z'Z+R'R = v d^2 v' + R'R
#          = v %*% (d^2 + (R%*%v)'(R%*%v))
  Rv <- (Rfac %*% Zsvd$v)
  d2.vR.Rv <- (diag(Zsvd$d^2)+crossprod(Rv))
  d2.eig <- eigen(d2.vR.Rv, symmetric=TRUE)
# Check for ill conditioning
  d2.ev <- d2.eig$values
  ZR.rank <- sum(d2.ev>0)
  {
    if(ZR.rank<1){
      warning("Z'Z+R'R = 0")
      Smat <- array(0, dim=c(ncoef, ncoef))
    }
    else{
      if(ZR.rank<ncoef)
        warning("Z'Z+R'R is not positive definite.  rank = ",
                ZR.rank, ";  dim = ", ncoef, "; increasing the ",
                ncoef-ZR.rank, " smallest eigenvalues ",
                "to improve numeric stability.")
      d2.evMin <- eps*d2.ev[1]
      ZR.rank1 <- sum(d2.ev >= d2.evMin)
      if(ZR.rank1 < ZR.rank)
        warning("Z'Z+R'R is ill conditioned.  Increasing the ",
                ncoef-ZR.rank1, " smallest eigenvalues ",
                "to improve numeric stability.")
#  Z'(inv(Z'Z+R'R)Z
#     = (u%*%d%*%w) solve(LAM) (udw)'
      udw <- ((Zsvd$u * rep(Zsvd$d, each=ncoef)) %*% d2.eig$vectors)
# or
#udw1<-((Zsvd$u[,1:ZR.rank1, drop=FALSE]*rep(Zsvd$d[1:ZR.rank1],e=ncoef))%*%d2.eig$vectors[1:ZR.rank1,,drop=FALSE])
      d2.ev2 <- pmax(d2.ev, d2.evMin)
      Smat <- ((udw / rep(d2.ev2, each=ncoef)) %*% t(udw))
# or
#Smat1<-((udw1/rep(d2.ev2, each=ncoef)) %*% t(udw1))
    }
  }
  df. <- sum(diag(Smat))
  dfe <- ncoef - df.
  gcv <- (ncoef/dfe)*sum(res1[1:ncoef]^2)/dfe
##
##  6.  compute fits and residuals
##
# [N, nbasis] = size(datstruct.basismat)
  ind1 <- 1:nbasis
  ind2 <- (nbasis+1):(2*nbasis)
  Ccoef <- coef[ind1]
  Tcoef <- coef[ind2]
#
  phimat <- datstruct$basismat
#
  Chat0 <- phimat %*% Ccoef
  That0 <- phimat %*% Tcoef
#
  yobs <- datstruct$y
#
  Cwt <- as.vector(datstruct$Cwt)
  Twt <- as.vector(datstruct$Twt)
#
#  res = []
#  if fit(1) res = [res; (yobs(:,1) - Chat0)./sqrt(Cwt)]
#  if fit(2)res = [res; (yobs(:,2) - That0)./sqrt(Twt)]
#
  yNames <- c("Conc", "Temp")
  fitNames <- yNames[as.logical(fit)]
  fit12 <- length(fitNames)
  basisNames <- dimnames(phimat)
  res <- array(NA, dim=c(N, fit12), dimnames=
                list(basisNames[[1]], fitNames))
  if(fit[1])res[, 1] <- ((yobs[,1] - Chat0)/sqrt(Cwt))
#  if(fit[2])res[, fit12] <- ((yobs[,fit12] - That0)/sqrt(Twt))
  if(fit[2])res[, fit12] <- ((yobs[,2] - That0)/sqrt(Twt))
#  res[1:5] matches Matlab 2007.05.29
##
## 7.  Derivatives?
##
  Dres <- Dres0
  if( gradwrd){
#
#  7.1.  set up basis and basis derivatve matrices
#
    quadmat <- datstruct$quadbasismat
    Dquadmat <- datstruct$Dquadbasismat
#
#    [nquad, nbasis] = size(quadmat)
#    onesb = ones(1,nbasis)
#    onesq = ones(nquad, 1)
    onesb <- rep(1, nbasis)
    onesq <- rep(1, nquad)
#
#  7.2.  set up some constants that are required
#
    V      <- fitstruct$V
    rho    <- fitstruct$rho
    rhoc   <- fitstruct$rhoc
    delH   <- fitstruct$delH
    Cp     <- fitstruct$Cp
    Cpc    <- fitstruct$Cpc
    Tref   <- fitstruct$Tref
#
#  7.3.  Set up input arrays
#
    F.  <- datstruct$F.
    CA0 <- datstruct$CA0
    T0  <- datstruct$T0
    Tc  <- datstruct$Tcin
    Fc  <- datstruct$Fc
#
#  7.4.  C and T values at fine grid
#
    Chat  <- as.vector(quadmat%*%Ccoef)
    That  <- as.vector(quadmat%*%Tcoef)
    DChat <- as.vector(Dquadmat%*%Ccoef)
    DThat <- as.vector(Dquadmat%*%Tcoef)
#
#  7.5.  betaCC and betaTC depend on kref and Eover R
    Tdif   <- 1/That - 1/Tref
#    temp   = exp(-1e4*EoverR*Tdif)
    log.temp <- (-1e4*EoverR*Tdif)
    oops <- (log.temp > max.log.betaCC)
    if(any(oops)){
      warning(sum(oops), " of ", length(log.temp),
              " values of (-1e4*EoverR*Tdif) exceed the max = ",
              max.log.betaCC, ";  thresholding.")
      log.temp[oops] <- max.log.betaCC
    }
    temp <- exp(log.temp)
#
    betaCC <- kref*temp
    TCfac  <- (-delH/(rho*Cp))
    betaTC <- TCfac*betaCC
#
#  7.6.  betaTT depends on a and b
    Fc2b    <- Fc^b
    aFc2b   <- a*Fc2b
    K1      <- V*rho*Cp
    K2      <- 1./(2.*rhoc*Cpc)
    betaTT  <- Fc*aFc2b/(K1*(Fc + K2*aFc2b))
    betaTT0 <- F./V
#
#  7.7.  compute derivatives of residuals
#
#    %  L values
#
    LC <- DChat + (betaTT0 + betaCC)*Chat - betaTT0*CA0*onesq
    LT <- DThat + ((betaTT0 + betaTT)*That - betaTC*Chat -
                 (betaTT0*T0 + betaTT*Tc))
#
#    %  first order derivatives of L values
#    %  derivatives of L values with respect to
#    %  coefficient vectors c and t
#
    DtbetaCC <- (1e4*EoverR/That^2)*betaCC
    DtbetaTC <- TCfac*DtbetaCC
#
    DcDChat  <- (-(betaCC + betaTT0))
    DtDChat  <- (-DtbetaCC*Chat)
    DcDThat  <-  betaTC
    DtDThat  <- (-(betaTT+betaTT0) + DtbetaTC*Chat)
#
    DcLC   <- (Dquadmat - outer(DcDChat, onesb)*quadmat)
    DtLC   <-          - outer(DtDChat, onesb)*quadmat
    DcLT   <-          - outer(DcDThat, onesb)*quadmat
    DtLT   <- (Dquadmat - outer(DtDThat, onesb)*quadmat)
#
    quadwts    <- datstruct$quadwts
    rootwts    <- sqrt(quadwts)
    quadwtsmat <- outer(quadwts, onesb)
#
#    %  k derivatives
#
    lamC <- lambda[1]
    lamT <- lambda[2]

#  7.8.  assemble the Jacobian matrix

#    DLC <- sqrt(lamC/Cwt).*[DcLC, DtLC]
#    DLT <- sqrt(lamT/Twt).*[DcLT, DtLT]
    DLC <- sqrt(lamC/Cwt)*cbind(DcLC, DtLC)
    DLT <- sqrt(lamT/Twt)*cbind(DcLT, DtLT)
#
#    Jacobian <- [DLC; DLT]
    Jacobian <- rbind(DLC, DLT)
#
#  7.9.  compute derivatives with respect to parameters
#
#    %  set up right hand side of equation D2GDc
#
#    D2GDc <- []
    D2GDc. <- vector('list', 4)
    names(D2GDc.) <- c("kref", "EoverR", "a", "b")
#
#    %  kref
#
    if( estimate[1]) {
#
#        %  first derivative of L values
#
      DkbetaCC <- temp
      DkbetaTC <- TCfac*DkbetaCC
#
      DkLC <-  DkbetaCC*Chat
      DkLT <- -DkbetaTC*Chat
#
#        %  second derivative of L values
#
      DktbetaCC <- (1e4*EoverR/That^2)*temp
      DktbetaTC <- TCfac*DktbetaCC
#
      DkcLC <- outer(  DkbetaCC,   onesb)*quadmat
      DkcLT <- outer( -DkbetaTC, onesb)*quadmat
      DktLC <- outer( DktbetaCC*Chat, onesb)*quadmat
      DktLT <- outer(-DktbetaTC*Chat, onesb)*quadmat
#
#      D2GDck <- zeros(2*nbasis,1)
#
      D2GDck <- rep(0, 2*nbasis)
#        D2GDck(ind1,1) <- (lamC/Cwt).* ...
#            (DcLC'*(DkLC.*quadwts) + DkcLC'*(LC.*quadwts)) + ...
#            (lamT/Twt).* ...
#            (DcLT'*(DkLT.*quadwts) + DkcLT'*(LT.*quadwts))
      D2GDck[ind1] <- ((lamC/Cwt)*
            (t(DcLC)%*%(DkLC*quadwts) + t(DkcLC)%*%(LC*quadwts)) +
            (lamT/Twt)*
            (t(DcLT)%*%(DkLT*quadwts) + t(DkcLT)%*%(LT*quadwts)) )
#        D2GDck(ind2,1) <- (lamC/Cwt).* ...
#            (DtLC'*(DkLC.*quadwts) + DktLC'*(LC.*quadwts)) + ...
#            (lamT/Twt).* ...
#            (DtLT'*(DkLT.*quadwts) + DktLT'*(LT.*quadwts))
      D2GDck[ind2] <- ((lamC/Cwt)*
            (t(DtLC)%*%(DkLC*quadwts) + t(DktLC)%*%(LC*quadwts)) +
            (lamT/Twt)*
            (t(DtLT)%*%(DkLT*quadwts) + t(DktLT)%*%(LT*quadwts)) )
#
#        D2GDc <- [D2GDc, D2GDck]
      D2GDc.$kref <- D2GDck
#
#    end kref
    }

#    %  EoverR
    if(estimate[2]){
#
#        %  first derivative of L values
#
      Dtemp    <- (-1e4*kref*Tdif*temp)
      DEbetaCC <- Dtemp
      DEbetaTC <- TCfac*DEbetaCC
#
      DELC  <-  DEbetaCC*Chat
      DELT  <- (-DEbetaTC*Chat)
#
#      DEtbetaCC <- (1e4.*kref  ./That.^2).* ...
#            (1 - 1e4.*EoverR.*Tdif).*temp
      DEtbetaCC <- ((1e4*kref/That^2)*
            (1 - 1e4*EoverR*Tdif)*temp)
      DEtbetaTC <- TCfac*DEtbetaCC
#
#        DEcLC <- (  DEbetaCC        *onesb).*quadmat
#        DEcLT <- ( -DEbetaTC        *onesb).*quadmat
#        DEtLC <- (( DEtbetaCC.*Chat)*onesb).*quadmat
#        DEtLT <- ((-DEtbetaTC.*Chat)*onesb).*quadmat
      DEcLC <- outer(  DEbetaCC, onesb)*quadmat
      DEcLT <- outer( -DEbetaTC, onesb)*quadmat
      DEtLC <- outer(( DEtbetaCC*Chat), onesb)*quadmat
      DEtLT <- outer((-DEtbetaTC*Chat), onesb)*quadmat
#
#        D2GDcE <- zeros(2*nbasis,1)
      D2GDcE <- rep(0, 2*nbasis)
#        D2GDcE(ind1,1) <- (lamC/Cwt).* ...
#            (DcLC'*(DELC.*quadwts) + DEcLC'*(LC.*quadwts)) + ...
#            (lamT/Twt).* ...
#            (DcLT'*(DELT.*quadwts) + DEcLT'*(LT.*quadwts))
      DcLC.. <- (crossprod(DcLC, DELC*quadwts) +
                  crossprod(DEcLC, LC*quadwts))
      DcLT.. <- (crossprod(DcLT, DELT*quadwts) +
                  crossprod(DEcLT, LT*quadwts))
      D2GDcE[ind1] <- ((lamC/Cwt)* DcLC.. + (lamT/Twt)* DcLT..)
#        D2GDcE(ind2,1) <- (lamC./Cwt).* ...
#            (DtLC'*(DELC.*quadwts) + DEtLC'*(LC.*quadwts)) + ...
#            (lamT./Twt).* ...
#            (DtLT'*(DELT.*quadwts) + DEtLT'*(LT.*quadwts))
      DtLC.. <- (crossprod(DtLC, DELC*quadwts) +
                  crossprod(DEtLC, LC*quadwts))
      DtLT.. <- (crossprod(DtLT, DELT*quadwts) +
                  crossprod(DEtLT, LT*quadwts))
      D2GDcE[ind2] <- ((lamC/Cwt)* DtLC.. + (lamT/Twt)* DtLT..)
#        D2GDc <- [D2GDc, D2GDcE]
      D2GDc.$EoverR <- D2GDcE
#
#    end EoverR
    }

#    %  a
    if(estimate[3]){
#
#        %  first derivative of L values
#
#        DhbetaTT <- (betaTT./aFc2b).*(1 - K1.*K2.*betaTT./Fc)
#        DabetaTT <- DhbetaTT.*Fc2b
      DhbetaTT <- ((betaTT/aFc2b)*(1 - K1*K2*betaTT/Fc))
      DabetaTT <- DhbetaTT*Fc2b
#
      DaLT <- (DabetaTT*(That - Tc))
#
      DatLT <- outer(DabetaTT, onesb)*quadmat
#
#        D2GDca <- zeros(2*nbasis,1)
      D2GDca <- rep(0, 2*nbasis)
#
#        D2GDca(ind1,1) <- (lamT/Twt).*(DcLT'*(DaLT.*quadwts))
#        D2GDca(ind2,1) <- (lamT./Twt).* ...
#            (DtLT'*(DaLT.*quadwts) + DatLT'*(LT.*quadwts))
      D2GDca[ind1] <- ((lamT/Twt)*(t(DcLT) %*%(DaLT*quadwts)))
      D2GDca[ind2] <- ((lamT/Twt)*
            (t(DtLT)%*%(DaLT*quadwts) + t(DatLT)%*%(LT*quadwts)) )
#
#        D2GDc <- [D2GDc, D2GDca]
      D2GDc.$a <- D2GDca
#
#    end 'a'
    }

    if( estimate[4]){
#
#        %  b derivative of L values
#
#        DhbetaTT <- (betaTT./aFc2b).*(1 - K1.*K2.*betaTT./Fc)
#        DbbetaTT <- DhbetaTT.*b.*aFc2b./Fc
      DhbetaTT <- (betaTT/aFc2b)*(1 - K1*K2*betaTT/Fc)
      DbbetaTT <- DhbetaTT*b*aFc2b/Fc
#
      DbLT <- DbbetaTT*(That - Tc)
#
#        DbtLT <- (DbbetaTT*onesb).*quadmat
      DbtLT <- outer(DbbetaTT, onesb)*quadmat
#
#        D2GDcb <- zeros(2*nbasis,1)
#        D2GDcb(ind1,1) <- (lamT/Twt).*(DcLT'*(DbLT.*quadwts))
#        D2GDcb(ind2,1) <- (lamT./Twt).* ...
#            (DtLT'*(DbLT.*quadwts) + DbtLT'*(LT.*quadwts))
#
      D2GDcb <- rep(0, 2*nbasis)
      D2GDcb[ind1] <- ((lamT/Twt)*(t(DcLT)%*%(DbLT*quadwts)) )
      D2GDcb[ind2] <- ((lamT/Twt)*
            (t(DtLT)%*%(DbLT*quadwts) + t(DbtLT)%*%(LT*quadwts)) )
#
#        D2GDc <- [D2GDc, D2GDcb]
      D2GDc.$b <- D2GDcb
#
#    end 'b'
    }
#   Convert from a list to a matrix,
#   dropping columns not in estimate
    D2GDc <- do.call(cbind, D2GDc.)
##
## 8.  Construct D2GDc2
##
#  8.1.  First part
#    Wmat <- [quadwtsmat, quadwtsmat; quadwtsmat, quadwtsmat]
    W.5 <- cbind(quadwtsmat, quadwtsmat)
    Wmat <- rbind(W.5, W.5)
#
#    D2GDc2  <- (Jacobian.*Wmat)'*Jacobian
    D2GDc2  <- crossprod(Jacobian*Wmat, Jacobian)
#
#    ZtZmat <- phimat'*phimat
    ZtZmat <- crossprod(phimat)

    if( fit[1])
        D2GDc2[ind1,ind1] <- D2GDc2[ind1,ind1] + ZtZmat/Cwt
#    end
    if( fit[2])
        D2GDc2[ind2,ind2] <- D2GDc2[ind2,ind2] + ZtZmat/Twt
#    end

#  8.2.  Add second derivative information

    DttbetaCC <- ((1e4*kref*EoverR/That^2)*
                 (1e4*EoverR/That^2 - 2/That)*temp)
    DttbetaTC <- TCfac*DttbetaCC
#
#    DctLC <- sparse(zeros(nbasis,nbasis))
#    DttLC <- sparse(zeros(nbasis,nbasis))
#    DctLT <- sparse(zeros(nbasis,nbasis))
#    DttLT <- sparse(zeros(nbasis,nbasis))
    DctLC <- array(0, dim=c(nbasis, nbasis))
    DttLT <- DctLT <- DttLC <- DctLC

#    norder <- nbasis - length(getbasispar(CSTRbasis))
## *** 'getbasispar' <- interior knots
    norder <- nbasis - length(CSTRbasis$params)
#
    for( i in 1:nbasis){
      jstart <- max(c(1,i-norder+1))
      for( j in jstart:i) {
        qijvec <- quadmat[,i]*quadmat[,j]*quadwts
        DctLC[i,j] <- sum(qijvec*LC*DtbetaCC)
        DttLC[i,j] <- sum(qijvec*LC*DttbetaCC*Chat)
        DctLT[i,j] <- sum(qijvec*LT*DtbetaTC)
        DttLT[i,j] <- sum(qijvec*LT*DttbetaTC*Chat)
        if( i != j){
          DctLC[j,i] <- DctLC[i,j]
          DttLC[j,i] <- DttLC[i,j]
          DctLT[j,i] <- DctLT[i,j]
          DttLT[j,i] <- DttLT[i,j]
#       end if(i!=j)
        }
#    end for(j in jstart:i)
      }
#   end for(i in 1:nbasis)
    }
#    DctL <- lamC.*DctLC./Cwt + lamT.*DctLT./Twt
    DctL <- lamC*DctLC/Cwt + lamT*DctLT/Twt
#      DttL <- lamC.*DttLC./Cwt + lamT.*DttLT./Twt
    DttL <- lamC*DttLC/Cwt + lamT*DttLT/Twt

#  8.3.  modify D2GDc2

    D2GDc2[ind1,ind2] <- D2GDc2[ind1,ind2] + DctL
    D2GDc2[ind2,ind1] <- D2GDc2[ind2,ind1] + t(DctL)
    D2GDc2[ind2,ind2] <- D2GDc2[ind2,ind2] + DttL

#  8.4.  compute (D2GDc2)^{-1} D2GDc

#    DcDtheta <- D2GDc2\D2GDc
    DcDtheta <- try(solve(D2GDc2, D2GDc))
    if(class(DcDtheta)=="try-error"){
#
      D2GDc2.eig <- eigen(D2GDc2, symmetric=TRUE)
      Dc.ev <- D2GDc2.eig$values
      Dc.rank <- sum(Dc.ev>0)
#
      if(Dc.rank<ncoef)
        warning("D2GDc2 has reduced rank ", Dc.rank, "; using ginverse.")
      Dc.rank1 <- sum(Dc.ev > (eps*Dc.ev[1]))
      if(Dc.rank1 < Dc.rank)
        warning("D2GDc2 is ill conditioned.  Reducing rank to ",
                Dc.rank1, " from ", Dc.rank)
      jrank <- 1:Dc.rank1
      DcDtheta <- with(D2GDc2.eig, (vectors[, jrank] /
          rep(Dc.ev[jrank], each=ncoef)) %*% crossprod(vectors[, jrank], D2GDc))
    }
#
#  8.5.  set up Dres
#
#    Dres <- []
    Dres <- NULL
#
    if( fit[1])
      Dres <- phimat%*%DcDtheta[ind1,]/sqrt(Cwt)
#       end
    if( fit[2])
      Dres <- rbind(Dres, phimat%*%DcDtheta[ind2,]/sqrt(Twt))
  }
##
## 9.  Done
##
#  list(res=res, Dres=Dres, fitstruct=fitstruct, df=df, lambda=lambda, gradwrd=gradwrd)
  list(res=res, Dres=Dres, fitstruct=fitstruct, df=df., gcv=gcv)
}

