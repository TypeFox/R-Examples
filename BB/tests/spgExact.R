# Test for bad output when exact initial values are supplied (from John Nash)

require(BB)

maxfn<-function(x) {
	       n<-length(x)
       ss<-seq(1,n)
       f<-10-(crossprod(x-ss))^2
       f<-as.numeric(f)
       return(f)
}

x0<-c(1,2,3,4)
ans.spgm <- spg(x0, maxfn,control=list(maximize=TRUE,trace=TRUE), quiet=TRUE)

if ( 10.0 != ans.spgm$value)  stop(
   "spg is not returning correct value when given optimum starting point.")
   
if ( ! all( 0 == (ans.spgm$par - c(1,2,3,4))))  stop(
   "spg is not returning correct par when given optimum starting point.")
