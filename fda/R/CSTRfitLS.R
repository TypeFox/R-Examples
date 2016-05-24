CSTRfitLS <- function(coef, datstruct, fitstruct, 
                                 lambda, gradwrd=FALSE){
#
#function [res, Dres] = CSTRfitLS(coef, datstruct, fitstruct, ...
#                                 lambda, gradwrd)

#  Last modified 2007.05.10 by Spencer Graves
#% previously modified 9 May 2005
##
## 1.  Set up   
##
  max.log.betaCC <- (log(.Machine$double.xmax)/3)
# For certain values of 'coef',
# naive computation of betaCC will return +/-Inf,
# which generates NAs in Dres.
# Avoid this by clipping betaCC
#
# log(.Machine$double.xmax)/2 is too big,
# because a multiple of it is squared in CSTRfn ... 
#  
  fit = fitstruct$fit;

  basismat = datstruct$basismat;
# basisMat <- read.csv('CSTRbasismat.csv', header=FALSE)
#[N, nbasis] = size(datstruct.basismat);
  N <- dim(basismat)[1]
  nbasis <- dim(basismat)[2]
#           
  ind1    = 1:nbasis
  ind2    = ((nbasis+1):(2*nbasis))
# onesb is redefined later before it is used  
#  onesb <- rep(1, nbasis)
  zeromat <- array(0, dim=c(N, nbasis), 
            dimnames=dimnames(basismat) )
#
  Ccoef  = coef[ind1];
  Tcoef  = coef[ind2];
  
  Cwt  = as.vector(datstruct$Cwt)
  Twt  = as.vector(datstruct$Twt)
##
## 2.  Sres = matrix of residuals = 
##     datstruct$y - predicted 
##
# Sres = []
  yNames <- c("Conc", "Temp")
  fitNames <- yNames[as.logical(fit)]
  fit12 <- length(fitNames)

  yobs = datstruct$y;
  basisNames <- dimnames(basismat)
  Sres <- array(NA, dim=c(N, fit12), dimnames=
                list(basisNames[[1]], fitNames))

  if( fit[1]){
    resC = yobs[,1] - basismat%*%Ccoef;
    Sres[,1] <- resC/sqrt(Cwt) 
  }

  if( fit[2]){
#    resT = yobs[,fit12] - basismat%*%Tcoef;
    resT = yobs[,2] - basismat%*%Tcoef;
    Sres[,fit12] <- resT/sqrt(Twt)
  }
##
## 3.  Compute (Conc, Temp) and d/dt at quadrature points 
##
#  3.1.  Get basic model coefficients  
  kref   = fitstruct$kref;
  EoverR = fitstruct$EoverR;
  a      = fitstruct$a;
  b      = fitstruct$b;
    
#%  3.2.  basis function values at quadrature points

  quadmat  = datstruct$quadbasismat;
  Dquadmat = datstruct$Dquadbasismat;
# [nquad, nbasis] = size(quadmat);
  nquad <- dim(quadmat)[1]
  onesb = rep(1,nbasis)
  names(onesb) <- dimnames(quadmat)[[2]]
#  onesq = rep(nquad, 1);
  onesq <- rep(1, nquad)
  names(onesq) <- dimnames(quadmat)[[1]] 

#%  3.3.  set up the values of C and T at quad. pts.

  Chatquad = as.vector(quadmat%*%Ccoef)
  Thatquad = as.vector(quadmat%*%Tcoef)

  DC = Dquadmat%*%Ccoef;
  DT = Dquadmat%*%Tcoef;
##
## 4.  Right hand side of differential equation
##  

#% 4.1.  set up some constants that are required

  V      = fitstruct$V;
  rho    = fitstruct$rho;
  rhoc   = fitstruct$rhoc;
  delH   = fitstruct$delH;
  Cp     = fitstruct$Cp;
  Cpc    = fitstruct$Cpc;
  Tref   = fitstruct$Tref;
  
#%  these constants can vary.
#%  see function CSTR2in for other conditions

  Fc  = datstruct$Fc
  F.   = datstruct$F.
  CA0 = datstruct$CA0;
  T0  = datstruct$T0;
  Tc  = datstruct$Tcin;

#% 4.2.  compute multipliers of outputs

  Tdif    = 1./Thatquad - 1./Tref;
#
#  betaCC  = kref*exp(-1e4*EoverR*Tdif);
  log.betaCC <- (log(abs(kref))-1e4*EoverR*Tdif)
  oops <- (log.betaCC > max.log.betaCC)
  if(any(oops)){
    warning(sum(oops), " of ", length(log.betaCC),
            " values of log(abs(betaCC)) exceed the max = ",
            max.log.betaCC, ";  thresholding.")
    log.betaCC[oops] <- max.log.betaCC 
  }
  betaCC <- sign(kref)*exp(log.betaCC)
#                     
  TCfac   = -delH/(rho*Cp);
  betaTC  = TCfac*betaCC;
  aFc2b   = a*Fc^b;
  K1      = V*rho*Cp;
  K2      = 1./(2.*rhoc*Cpc);
  betaTT  = Fc*aFc2b/(K1*(Fc + K2*aFc2b));
  betaTT0 = F./V;

#% 4.3.  compute right sides of equations 

  DChat = (-(betaTT0 + betaCC)*Chatquad + betaTT0*CA0)
  DThat = (-(betaTT0 + betaTT)*Thatquad + betaTC*Chatquad
              + betaTT0*T0 + betaTT*Tc)

#  4.4.  Deviation between left and right hand sides   
  LC = DC - DChat
  LT = DT - DThat

#  4.5.  Quadrature weights   
  quadwts = datstruct$quadwts;
  rootwts = sqrt(quadwts);

#  lambdaC = lambda[1];
#  lambdaT = lambda[2];
  lambdaC.5 = sqrt(lambda[1]/Cwt) 
  lambdaT.5 = sqrt(lambda[2]/Twt)
##
## 5.  Lres = scaled deviations of Dy from predicted 
##    
#  Lres <- rbind(LC*rootwts*sqrt(lambdaC/Cwt),
#                LT*rootwts*sqrt(lambdaT/Twt) )
  Lres <- cbind(LConc = as.vector(LC)*rootwts*lambdaC.5,
                LTemp = as.vector(LT)*rootwts*lambdaT.5)

##
## 6.  Combine Sres and Lres
##  
  res  = list(Sres=Sres, Lres=Lres);
  out <- list(res=res) 
##
##% 7. compute gradient if required
##
  if(gradwrd){

#  7.1.  Derivatives of fit residuals
    
#    DSres = [];
#     if(fit[1])DSres = cbind( -basismat./sqrt(Cwt),zeromat)
#    if( fit[2])DSres = c(DSres, zeromat, -basismat./sqrt(Twt));
    DSres <- NULL
    if(fit[1])
      DSres <- cbind(-basismat/sqrt(Cwt), zeromat)
    if(fit[2])
      DSres <- rbind(DSres, cbind(zeromat, -basismat/sqrt(Twt)))
#  DSresMat <- read.csv('CSTR-DSres.csv', header=FALSE)
# d.DSres <- DSres - as.matrix(DSresMat)
# quantile(d.DSres)    
#           0%           25%           50%           75%          100% 
#-3.971928e-06  0.000000e+00  0.000000e+00  0.000000e+00  1.573097e-05
# sqrt(mean(DSres^2)) =0.218
# Reasonable accuracy.      
    
#  7.2.  Derivatives of weight functions
    
    DtbetaCC = (1e4*EoverR/Thatquad^2)*betaCC;
    DtbetaTC = TCfac*DtbetaCC;
    
#  7.3.  Derivatives of RHS of operators
    
    DcDChat  = -(betaCC + betaTT0);
    DtDChat  = -DtbetaCC*Chatquad;
    DcDThat  =  betaTC;
    DtDThat  = -(betaTT+betaTT0) + DtbetaTC*Chatquad;
    
#  7.4.  Operator derivatives
    
#    DcLC   = Dquadmat - (DcDChat*onesb).*quadmat;
#    DtLC   =          - (DtDChat*onesb).*quadmat;
#    DcLT   =          - (DcDThat*onesb).*quadmat;
#    DtLT   = Dquadmat - (DtDThat*onesb).*quadmat;
    DcLC   = Dquadmat - outer(DcDChat, onesb)*quadmat;
#    DcLC   = Dquadmat - tcrossprod(DcDChat, onesb)*quadmat;
    DtLC   =          - outer(DtDChat, onesb)*quadmat;
    DcLT   =          - outer(DcDThat, onesb)*quadmat;
    DtLT   = Dquadmat - outer(DtDThat, onesb)*quadmat;
    
#  7.5.  Multiply operator derivatives by root of
#    %  quadrature weights over root of SSE weights
    
#    wtmat   = rootwts*onesb;
#    DcLCmat = DcLC.*wtmat.*sqrt(lambdaC./Cwt);
#    DtLCmat = DtLC.*wtmat.*sqrt(lambdaC./Cwt);
#    DcLTmat = DcLT.*wtmat.*sqrt(lambdaT./Twt);
#    DtLTmat = DtLT.*wtmat.*sqrt(lambdaT./Twt);
    wtmat   = outer(rootwts, onesb) 
    DcLCmat = DcLC*wtmat*lambdaC.5
    DtLCmat = DtLC*wtmat*lambdaC.5
    DcLTmat = DcLT*wtmat*lambdaT.5
    DtLTmat = DtLT*wtmat*lambdaT.5

#  7.6.  Matrices of derivative of operator residuals
    
#    DLres  = [[DcLCmat, DtLCmat]; [DcLTmat, DtLTmat]]
    DLres <- rbind(cbind(DcLCmat, DtLCmat),
                   cbind(DcLTmat, DtLTmat))
#  DLresMat <- read.csv('CSTR-DLres.csv', header=FALSE)
# d.DLres <- DLres - as.matrix(DLresMat)
# quantile(d.DLres)    
#           0%           25%           50%           75%          100% 
#-0.0006900152  0.0000000000  0.0000000000  0.0000000000  0.0006814540 
# sqrt(mean(DLres^2))
# 0.812
# Not great but acceptable.      
    
#  7.7.  Combine with derivative of fit residuals
    
    out$Dres   = list(DSres=DSres, DLres=DLres)
    
  }
  
#  else
#     Dres = vector("list", 0)

  return(out)
}
