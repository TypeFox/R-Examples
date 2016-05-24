bmstep <- function(par, srchdirn, lower=NULL, upper=NULL, bdmsk=NULL, trace=0) {  
## Find maximum steplength from par along srchdirn given bounds and masks
# ?? do we want step length as multiple of srchdirn?, or actual max step length
# 20101031 -- issue of bounds not working correctly
#  - -Inf seems to upset bounds
#  - no upper bounds gives troubles (applies to Rcgmin too!)
#
# Input:
#  par = a vector containing the starting point
#  srchdirn = the direction of search (a vector)
#  lower = vector of lower bounds on parameters
#  upper = vector of upper bounds on parameters
#    Note: free parameters outside bounds will be adjusted to bounds.
#  bdmsk = control vector for bounds and masks. Parameters for which bdmsk are 1
#         are unconstrained or "free", those with bdmsk 0 are masked i.e., fixed.
#         For historical reasons, we use the same array as an indicator that a
#         parameter is at a lower bound (-3) or upper bound (-1)
#  trace = control of output: 0 for none (default), >0 for output
##
# Output:
#    A double giving the maximum steplength. Not bigger than maxstep.
#
########## length of vectors #########
n<-length(par)
############# bounds and masks ################
# check if there are bounds
  if(is.null(lower) || ! any(is.finite(lower))) nolower<-TRUE else nolower<-FALSE
  if(is.null(upper) || ! any(is.finite(upper))) noupper<-TRUE else noupper<-FALSE
# Next line NOT same as in bmchk(). Leave out bdmsk.
  if(nolower && noupper) bounds<-FALSE else bounds<-TRUE
  if (is.null(bdmsk)) bdmsk<-rep(1,n) # make sure we have values
  if (any(bdmsk==0)) {
     if (trace > 2) cat("Masks present -- adjusting search direction.\n")
     srchdirn[which(bdmsk==0)]<-0 # adjust search direction for masked elements
  }
  if(nolower) lower<-rep(-Inf,n)
  if(noupper) upper<-rep(Inf,n)
######## find maximum step (may be Inf) #############
# distance to bounds
  d2lo<-par-lower
  d2up<-upper-par
  if (trace>0) {
     cat("Distances to bounds, lower then upper\n")
     print(d2lo)
     print(d2up)
  }
  sslo<-rep(0,n)
  ssup<-sslo
  # Now want to get ssup -- stepsize to upper bound along directions where srchdirn>0
# Hard way, by loop
  for (i in 1:n) {
      if (bdmsk[i]==1) { # free parameter
          sdi<-srchdirn[i]
          if (sdi>0) ssup[i]<-d2up[i]/sdi
          if (sdi<0) sslo[i]<- -d2lo[i]/sdi
          # sdi==0, no changes
      }
  }
# another approach 20111022
#  suppressWarnings(ssup2<-(bdmsk*(d2up/srchdirn)))
#  suppressWarnings(sslo2<-((-1)*bdmsk*(d2lo/srchdirn)))
#  cat("sslo2 & ssup2\n")
#  print(sslo2)
#  print(ssup2)
#  cat("after adjustment\n")
#  ssup2[which(srchdirn<=0)]<-0
#  sslo2[which(srchdirn>=0)]<-0
#  print(sslo2)
#  print(ssup2)
#  ss<-c(sslo2, ssup2)
#  ss<-ss[which(ss>0)]
#  ms2<-min(ss)
#  cat("ms2=",ms2,"\n")

  if (trace>0) {
     cat("steplengths, lower then upper\n")
     print(sslo)
     print(ssup)
  }
  sslo<-sslo[which(sslo>0)]  
  ssup<-ssup[which(ssup>0)]  
  if (trace>0) {
     cat("steplengths, truncated, lower then upper\n")
     if (length(sslo)>0) print(sslo) else cat("sslo NULL\n")
     if (length(ssup)>0) print(ssup) else cat("ssup NULL\n")
  }
  if (is.null(sslo) && is.null(ssup)) {# Not needed, min will return Inf
       maxstep<-Inf
  } else {
       maxstep<-min(sslo,ssup)
  }
} ## end of bmstep.R
