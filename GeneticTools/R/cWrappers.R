# C Wrapper for the MDR call:
  mdr.C <- function(X, fold, status,t,cv,cvp,top,na,fix){
    .Call( "mdr", X, fold, status, t, cv, cvp,top,na,fix)
  } 

# C Wrapper for the distance calculations:
  calcDistances.C <- function(X){
    .Call("calcDistances",X)
  }