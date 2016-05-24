## So I have worked out the list of functions making up the various
## cline models.  I have center, width, direction, upper step delta
## (d2), upper tau (tau2), lower step delta (d1), lower tau (tau1),
## pMax, and pMin as parameters.  I have helper parameters gamma, u,
## intercept (A), and kappa.


## Frame level methods.

## ## Calling the long form of this function will compensate for the
## ## scaling caused by using pMin and pMax.
## eval.gamma <- function(width,pMin=0,pMax=1) {
##   return( 4/(width*(pMax-pMin))) }

## ## eval.gammaAlt <- 4/width

## ## For now, direction must either be 1 or -1.  This value must also
## ## always be fixed (at least for these cline models).
## eval.lambda <- function(gamma,direction) { return(gamma*direction) }

## ## Unnecessary helper function to document implict evaluation.
## eval.gammaStep <- function( gamma, step ) { return(gamma * step) }

eval.intercept <- function( gammaStep ) {
  return( 1/(1+exp(gammaStep))) }

eval.kappa <- function( tau, gammaStep ) {
  return( tau /(1+exp(-gammaStep)))}

## Model evaluation methods.

## ## intermediate functions.

## helper.u        <- function(x,center,lambda){
##   return((x- center) * lambda) }

## ## lowerU and upperU included for the sake of clarity.
## helper.lowerU   <- function(u,lowerStep){
##   return(u + lowerStep)}

## helper.upperU   <- function(u,upperStep){
##   return(u - upperStep)}


## cline model functions
cline.pCenter   <- function(u) {
  return(1/(1+ exp(-u))) }

cline.pLower    <- function(lowerU,lowerA,lowerK) {
  return( lowerA * exp( lowerK * lowerU))}

cline.pUpper <- function(upperU,upperA,upperK) {
  return( 1 - upperA * exp( -upperK * upperU)) }


step1VectorExpF <- function( conditionalExp, trueExp, falseExp )
  substitute(ifelse( cE, tE , fE ) ,
             list(cE=conditionalExp,tE = trueExp, fE=falseExp))

cline.exp.scale <- function(f) bquote(pMin+(pMax-pMin)*.(f))
cline.exp.sigmoid <- function(u) bquote(1/(1+ exp(-.(u))))

cline.exp.pLower    <- function(lowerU,lowerA,lowerK) bquote( .(lowerA) * exp( .(lowerK) * .(lowerU)))

cline.exp.pUpper <- function(upperU,upperA,upperK) bquote( 1 - .(upperA) * exp( -.(upperK) * .(upperU))) 
## Frame evaluation methods

## gamma     <- eval.gamma(width,pMin, pMax)

## lambda    <- eval.lambda(gamma,direction)
cline.exp.lower <- function(u,d1,tau1)
  cline.exp.pLower( bquote( .(u) + 4*.(d1)/width),
                   bquote( 1/(1+exp(4*.(d1)/width))),
                   bquote( .(tau1) /(1+exp(-4*.(d1)/width)) ))

cline.exp.upper <- function(u,d1,tau1)
  cline.exp.pUpper(bquote( .(u) - 4*.(d1)/width),
                   bquote( 1/(1+exp(4*.(d1)/width))),
                   bquote( .(tau1) /(1+exp(-4*.(d1)/width)) ))
##step1VectorExpF
cline.exp.stepBoth <- function(u,du,d.lo,d.up,lowerTail,upperTail)
  step1VectorExpF(bquote(.(d.lo) < - .(du)),
                  lowerTail,
                  step1VectorExpF(bquote( .(d.up) < .(du)),
                                  upperTail,
                                  cline.exp.sigmoid(u)))
cline.exp.stepLow <-function(u,du,d.lo,lowerTail)
  step1VectorExpF(bquote(.(d.lo) < - .(du)),
                  lowerTail,
                  cline.exp.sigmoid(u))
cline.exp.stepUp <-function(u,du,d.up,upperTail)
  step1VectorExpF(bquote( .(d.up) < .(du)),
                  upperTail,
                  cline.exp.sigmoid(u))

meta.tail.lower <- function(gamma,d1,tau1){
  lowerStep <- gamma*d1
  lowerA    <- eval.intercept(lowerStep)
  lowerK    <- eval.kappa(tau1,lowerStep)
  lower.tail <- function(u) {
    lowerU <- u + lowerStep
    return ( (lowerU<0)*( lowerA * exp( lowerK * lowerU)))}
  return(list(lowerStep=lowerStep,
              lowerA=lowerA,
              lowerK=lowerK,
              lower.tail=lower.tail))
}

meta.tail.upper <- function(gamma,d2,tau2){
  upperStep <- gamma*d2
  upperA    <- eval.intercept(upperStep)
  upperK    <- eval.kappa(tau2,upperStep)
  upper.tail <- function(u) {
    upperU <- u - upperStep
    return ( (upperU>0)*( 1 - upperA * exp( -upperK * upperU)))}
  return(list(upperStep=upperStep,
              upperA=upperA,
              upperK=upperK,
              upper.tail=upper.tail))
}

meta.tail.mirror <- function(gamma,d,tau){
  step <- gamma*d
  A    <- eval.intercept(step)
  K    <- eval.kappa(tau,step)
  lower.tail <- function(u) {
    lowerU <- u + step
    return ( (lowerU<0)*( A * exp( K * lowerU)))}
  upper.tail <- function(u) {
    upperU <- u - step
    return ( (upperU>0)*( 1 - A * exp( -K * upperU)))}
  return(list(lowerStep=step,
              lowerA=A,
              lowerK=K,
              lower.tail=lower.tail,
              upperStep=step,
              upperA=A,
              upperK=K,
              upper.tail=upper.tail))
}

meta.cline.func.stepBoth <- function(center, direction, gamma,
                                     lowerTail, upperTail){
  loStep=lowerTail$lowerStep
  upStep=upperTail$upperStep
  loTailFunc=lowerTail$lower.tail
  upTailFunc=upperTail$upper.tail
  cline.pCenter  <- function(u) {
    inCenter= (-loStep <= u) & (u <= upStep )
    return((inCenter)*(1/(1+ exp(-u)))) }
  cline.func <- function(x) {
    u <- (x - center) * gamma * direction
    return(cline.pCenter(u)+loTailFunc(u)+upTailFunc(u)) }
  return(cline.func)
}

meta.cline.func.noStep <- function(center, direction, gamma) {
    cline.func <- function(x) {
      u <- (x - center) * gamma * direction
      return((1/(1+ exp(-u)))) }
    return(cline.func)
}

meta.cline.func.lowStep <- function(center, direction, gamma, lowerTail){
  loStep=lowerTail$lowerStep
  loTailFunc=lowerTail$lower.tail
  cline.pCenter  <- function(u) {
    inCenter= -loStep <= u
    return((inCenter)*(1/(1+ exp(-u)))) }
  cline.func <- function(x) {
    u <- (x - center) * gamma * direction
    return(cline.pCenter(u)+loTailFunc(u)) }
  return(cline.func)
}

meta.cline.func.upStep <- function(center, direction, gamma, upperTail){
  upStep=upperTail$upperStep
  upTailFunc=upperTail$upper.tail
  cline.pCenter  <- function(u) {
    inCenter= (u <= upStep )
    return((inCenter)*(1/(1+ exp(-u)))) }
  cline.func <- function(x) {
    u <- (x - center) * gamma * direction
    return(cline.pCenter(u)+upTailFunc(u)) }
  return(cline.func)
}

meta.cline.func.pScale <- function(pMin,pMax,cline.func){
  new.cline.func<-function(x){
    return(pMin+(pMax-pMin)*cline.func(x))}
  return(new.cline.func)
}

## posJunk  <- function(junk, func) {
## junk.num <- as.numeric(junk)
## zeros <- numeric(length(junk.num))
## return(ifelse(junk>0,func(junk),zeros))
## }

# Edit on Sunday, January 2nd 2011: finished that last function, need
# to transfer it to a R-code file.
