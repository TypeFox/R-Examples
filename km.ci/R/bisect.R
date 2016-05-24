"bisect" <-
function(fun,opar,lval,uval,tol=1e-7) {
 t1<-fun(lval,opar)
 t2<-fun(uval,opar)
 if(t1*t2 > 0 ) 
  stop("in bisect both function values have the same sign")
 if(t1>0) {
  t2<-uval
  uval<-lval
  lval<-t2
 }
 converged<-F
 while( !converged ) {
  t1<-(lval+uval)/2
  nf<-fun(t1,opar)
  if(abs(nf) < tol )
   return(t1)
  else if( nf<0)
   lval<-t1
  else
   uval<-t1
 }
}

