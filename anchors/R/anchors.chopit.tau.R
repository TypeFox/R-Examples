#######################################################################
##
## Function: anchors.chopit.tau()
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created :  2008-04-20
##
## Extracted and refined from former chopit()  
#######################################################################
anchors.chopit.tau <- function(v0,gamma,v1,gamma1,toffset,
                               nobs, n.cat,n.tau.set,
                               nvars.gamma, tau.start.idx,
                               use.linear,debug,do.gr,verbose)
{
  
  if (debug > 1) cat("test chopit: use.linear?",use.linear,"\n")
  if (debug > 1) cat("test chopit:",
                               dim(v0),":",
                               dim(v1),":",
                               nobs, n.cat, n.tau.set, nvars.gamma,"\n")
  if (debug > 1) cat("test chopit:",tau.start.idx,"\n")

  gamma  <- matrix(gamma,nrow=nvars.gamma)

  if (debug > 1) print(gamma)
                                        #    if (debug > 1) print(v0)
  
  if (use.linear) {
    tau    <- cbind(v1 %*% gamma1 + toffset,  v0 %*% gamma)
  } else {
    tau    <- cbind(v1 %*% gamma1 + toffset,  as.matrix(exp(v0 %*% gamma)))

                                        #      vg1  <- v0 %*% gamma         ## vg: n x (n.cat-1)*n.self
                                        #      tau  <- as.matrix(exp(vg1) ) ## exponentiate all ## TT
                                        #      ## then sub back in non-exponentiated data
                                        #      tau[,tau.start.idx]  <- vg1[,tau.start.idx]
  }
  if (debug > 1) cat("test chopit: create taus...\n")

  ## verify that additions are positive, ELSE RETURN NULL
  if (any(tau[,-tau.start.idx] <= 0) || any(!is.finite(tau))) {
    if (debug > 1) {
      cat("test chopit: failure on tau...\n")
      print(tau)
    }
    return(NULL)
  }

  if (debug > 1) cat("test chopit: dim",dim(tau),"\n")

  
  ## create cumulative taus....
  taus <- tau
  taus <- Crowcumsum(taus,nobs,n.cat-1,n.tau.set);

  if (debug > 1) cat("test chopit: done creating taus...\n")
  
                                        #    print("TAUS XXX")
                                        #    cat(n,n.cat-1,tmp,"\n")

  
  return(taus)
}
