# ACSPlan class definition and accessor function
setClass("ACSPlan", representation=representation(n="numeric",k="numeric"))

setGeneric("n",function(object) standardGeneric("n"))
setMethod("n","ACSPlan",function(object) object@n)

setGeneric("k",function(object) standardGeneric("k"))
setMethod("k","ACSPlan",function(object) object@k)


# plot method definition for objects of ACSPlan class
setMethod("plot", signature(x = "ACSPlan"),
      function(x, typeOC=c("exact", "napprox","ewmaSK","ewma2"),lam=1,xl=0.001, xu=0.1, xlabm="p",ylabm="L(p)",typem="l",...) {p=seq(from=xl,to=xu,length.out=1000);plot(x=p,y=OC(p,n(x),k(x), type=typeOC,lam),  xlab=xlabm,ylab=ylabm,type=typem, ...)     })

