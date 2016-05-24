myoptim <- function(par, fn, gr,method,
                          lower, upper,
                          control, hessian,...)
{
  if(length(lower)==1) lower <- rep(lower,length(par))
  if(length(upper)==1) upper <- rep(upper,length(par))
  thenlm <- suppressWarnings(nlminb(par, fn, lower = lower, upper = upper))
  #if (thenlm$convergence!=0) print(paste("error ",thenlm$message))
  if (thenlm$convergence==0) theoptim <- list(par=thenlm$par,value=thenlm$objective,convergence=thenlm$convergence,message=thenlm$message)
  else {
    thenm <- Nelder_Mead(fn, thenlm$par, lower = lower, upper = upper)
    theoptim <- list(par=thenm$par,value=thenm$fval,convergence=ifelse(thenm$NM.result<0,thenm$NM.result,0),message=thenm$message)
  }
  return(theoptim)
}