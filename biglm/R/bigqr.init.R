"bigqr.init" <-
function(p){

  rval<-list(D=numeric(p), rbar=numeric(choose(p,2)),
             thetab=numeric(p),
             ss=0, checked=FALSE,
             tol=numeric(p))
  class(rval)<-"bigqr"
  rval
}

