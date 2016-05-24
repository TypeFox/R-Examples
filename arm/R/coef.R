
coef.sim <- function(object,...){
    ans <- object@coef
    return(ans)
}


coef.sim.polr <- function(object, slot=c("ALL", "coef", "zeta"),...){
  slot <- match.arg(slot)
  if(slot=="coef"){
    ans <- object@coef
  } else if(slot=="zeta"){
      ans <- object@zeta
  } else {
      ans <- cbind(object@zeta, object@coef)
  }
  return(ans)
}

coef.sim.merMod <- function(object,...){
  fef <- object@fixef
  ref <- object@ranef
  ans <- list("fixef" = fef, "ranef" = ref)
  return(ans)
}


fixef.sim.merMod <- function(object,...){
 ans <- object@fixef
 return(ans)
}

ranef.sim.merMod <- function(object,...){
 ans <- object@ranef
 return(ans)
}
