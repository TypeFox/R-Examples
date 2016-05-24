
clineLogit <- quote(4*(x-center)/width)
clineLogitRev <- quote(-4*(x-center)/width)

alphaExpF <-  function ( logitDelta , tailTau )
  substitute( exp(lD*(tT-1))/(1+exp(-lD)), list(lD=logitDelta, tT=tailTau));

pTailExpF <-  function( logitDelta , tailTau , logitU)
  substitute(alpha * exp( lU * tT) ,
             list(alpha= alphaExpF(logitDelta, tailTau),
                  tT=tailTau,
                  lU=logitU))

qTailExpF <-  function( logitDelta , tailTau , logitU)
  substitute(1 - alpha * exp( lU * tT) ,
             list(alpha= alphaExpF(logitDelta, tailTau),
                  tT=tailTau,
                  lU=logitU))

vTailExpF <-  function( logitDelta , tailTau , logitU)
  substitute( alpha*exp(lU*tT) - alpha^2*exp(2*lU*tT),
             list(alpha= alphaExpF(logitDelta, tailTau),
                  tT=tailTau,
                  lU=logitU))

sigmoidExpF <- function( logitU)
  substitute( 1 / ( 1 + exp(lU) ) , list ( lU=logitU))

vCenterExpF <- function( logitU)
  substitute(  1 / (2+exp(lU)+exp(-lU)), list ( lU=logitU))

muExpF <-  function( minMu, maxMu , clineExp)
  substitute( muA + (muB -muA)*cE , list(muA=minMu, muB=maxMu, cE=clineExp))


## 
## qV <- list()
## qV$MuL <-  quote(muL)
## qV$MuR <-  quote(muR)
## qV$DL <-  quote(deltaL)
## qV$DR <-  quote(deltaR)
## qV$LDL <-  quote(4*deltaL/width)
## qV$LDR <-  quote(4*deltaR/width)
## qV$TL <-  quote(tauL)
## qV$TR <-  quote(tauR)
## qV$DM <-  quote(deltaM)
## qV$LDM <-  quote(4*deltaM/width)
## qV$TM <-  quote(tauM)
## qV$VL <-  quote(varL)
## qV$VR <-  quote(varR)
## qV$VH <-  quote(varH)
## qV$KSQ <-  quote(kappa*kappa)
## leftTailMuExp <- muExpF( qV$MuL, qV$MuR, pTailExpF( qV$LDL, qV$TL,clineLogit))
## rightTailMuExp <- muExpF(qV$MuR,
##                          qV$MuL,
##                          pTailExpF( qV$LDR, qV$TR,clineLogitRev))
varExpF <-  function( leftVar , rightVar , kappaE, clineExp, vClineExp)
  substitute (vA + (vB - vA)*cE+4*kE*vE ,
              list(vA=leftVar,
                   vB=rightVar,
                   kE=kappaE,
                   cE=clineExp,
                   vE=vClineExp))
## leftTailVarExp <- varExpF(qV$VL,
##                           qV$VR,
##                           qV$VH,
##                           pTailExpF( qV$LDL, qV$TL,clineLogit),
##                           vTailExpF( qV$LDL, qV$TL,clineLogit))
                          
## pCenterMuExp <-  muExpF( qV$MuL, qV$MuR, sigmoidExpF( clineLogitRev))

## qV$LTCe <-  quote(x < center  - deltaL)

## qV$RTCe <-  quote(x > center  + deltaR)


## qV$MLTCe <-  quote(x < center  - deltaM)

## qV$MRTCe <-  quote(x > center  + deltaM)


guassianLLExpF <-  function( tValues, muExp, varExp )
  substitute( (-log(2*pi*vE) - (tV - mE)^2/ vE)/2,
             list(tV=tValues,mE=muExp,vE=varExp))
guassianLLSampleExpF <- function( sampleMean, sampleVar, nEff, muExp, varExp)
  substitute(- nE*(log(2*pi*vE)+(sV+(sM-mE)^2)/vE)/2,
             list(sM=sampleMean, sV=sampleVar,nE=nEff,mE=muExp,vE=varExp))
guassianThetaLLExpF <- function( sampleMean, sampleVar, nEff,distance, muExp, varExp)
  substitute(- sum(nE*(log(2*pi*vE)+(sV+(sM-mE)^2)/vE))/2,
             list(sM=sampleMean, sV=sampleVar,nE=nEff,
                  mE=eval(substitute(substitute(a,list(x=distance)),list(a=muExp))),
                  vE=eval(substitute(substitute(a,list(x=distance)),list(a=varExp)))))

