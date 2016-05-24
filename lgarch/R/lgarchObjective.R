lgarchObjective <-
function(pars, aux)
{
  #check parameters:
  if( any(is.na(pars)) || any(pars<=aux$lower) || any(pars>=aux$upper) ){
    chk.conds <- FALSE
  }else{
    chk.conds <- TRUE
  } #end if..else..

  if(chk.conds){
    #if mean-correction:
    if(aux$mean.correction){
      pars <- c(0,pars)
    }
    #recursion:
    uadj <- lgarchRecursion1(pars, aux)
    if(aux$yzeron > 0){
      uadj <- uadj[-aux$yzerowhere]
    }
    #compute objective value:
    if(aux$method=="ls"){
      objective.value <- sum(uadj^2)
    }
    if(aux$method=="ml"){
      sigma2u <- pars[aux$sigma2u.indx]
      objective.value <- -aux$ynonzeron*log(sigma2u)/2 - aux$ynonzeron*log(2*pi)/2 - sum(uadj^2/sigma2u)/2
    }
    if(aux$method=="cex2"){
      uadjElnz2 <- uadj + pars[aux$sigma2u.indx]
      objective.value <- -aux$ynonzeron*log(2*pi)/2 + sum(uadjElnz2 - exp(uadjElnz2))/2
    }
    #check objective value:
    if(is.na(objective.value) || abs(objective.value) == Inf){
      objective.value <- aux$objective.penalty
    }
  }else{
    objective.value <- aux$objective.penalty
  } #end if(chk.conds)
  return(objective.value)
}
