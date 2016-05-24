### R code from vignette source 'brob.Rnw'

###################################################
### code chunk number 1: setClass
###################################################



###################################################
### code chunk number 2: brob.Rnw:120-129
###################################################
setClass("swift",
         representation = "VIRTUAL"
         )

setClass("brob",
         representation = representation(x="numeric",positive="logical"),
         prototype      = list(x=numeric(),positive=logical()),
         contains       = "swift"
         )


###################################################
### code chunk number 3: new
###################################################
new("brob",x=1:10,positive=rep(TRUE,10))


###################################################
### code chunk number 4: new_flaky_arguments
###################################################
new("brob",x=1:10,positive=c(TRUE,FALSE,FALSE))


###################################################
### code chunk number 5: validity_method
###################################################
.Brob.valid <- function(object){
  len <- length(object@positive)
  if(len != length(object@x)){
    return("length mismatch")
  } else {
    return(TRUE)
  }
}


###################################################
### code chunk number 6: call_setValidity
###################################################
setValidity("brob", .Brob.valid)


###################################################
### code chunk number 7: brob_definition
###################################################
"brob" <- function(x=double(),positive){
  if(missing(positive)){
    positive <- rep(TRUE,length(x))
  }
  if(length(positive)==1){
    positive <- rep(positive,length(x))
  }
  new("brob",x=as.numeric(x),positive=positive)
}


###################################################
### code chunk number 8: call_brob_recycling
###################################################
brob(1:10,FALSE)


###################################################
### code chunk number 9: use.function.is
###################################################
is(brob(1:5),"brob")


###################################################
### code chunk number 10: is.brob_definition
###################################################
is.brob <- function(x){is(x,"brob")}
is.glub <- function(x){is(x,"glub")}


###################################################
### code chunk number 11: as.brob_definition
###################################################
"as.brob" <- function(x){
  if(is.brob(x)){
    return(x)
  } else if(is.complex(x)) {
    warning("imaginary parts discarded")
    return(Recall(Re(x)))
  } else if(is.glub(x)){
    warning("imaginary parts discarded")
    return(Re(x))
  } else {
    return(brob(log(abs(x)), x>=0))
  }
}


###################################################
### code chunk number 12: as.brob_call
###################################################
as.brob(1:10)


###################################################
### code chunk number 13: setAs
###################################################



###################################################
### code chunk number 14: brob.Rnw:366-371
###################################################
setAs("brob", "numeric", function(from){
  out <- exp(from@x)
  out[!from@positive] <- -out[!from@positive]
  return(out)
} )


###################################################
### code chunk number 15: setMethodbrob
###################################################
setMethod("as.numeric",signature(x="brob"),function(x){as(x,"numeric")})


###################################################
### code chunk number 16: setAsbrobcomplex
###################################################
setAs("brob", "complex", function(from){
  return(as.numeric(from)+ 0i)
} )

setMethod("as.complex",signature(x="brob"),function(x){as(x,"complex")})


###################################################
### code chunk number 17: asCheck
###################################################
x <- as.brob(1:4)
x
as.numeric(x)


###################################################
### code chunk number 18: print_methods
###################################################
.Brob.print <- function(x, digits=5){
     noquote( paste(c("-","+")[1+x@positive],"exp(",signif(x@x,digits),")",sep=""))
   }


###################################################
### code chunk number 19: print.brob
###################################################
print.brob <- function(x, ...){
  jj <- .Brob.print(x, ...)
  print(jj)
  return(invisible(jj))
}


###################################################
### code chunk number 20: setmethodbrobshow
###################################################
setMethod("show", "brob", function(object){print.brob(object)})


###################################################
### code chunk number 21: as.brob14
###################################################
as.brob(1:4)


###################################################
### code chunk number 22: get.n.set
###################################################



###################################################
### code chunk number 23: brob.Rnw:465-469
###################################################
setGeneric("getX",function(x){standardGeneric("getX")})
setGeneric("getP",function(x){standardGeneric("getP")})
setMethod("getX","brob",function(x){x@x})
setMethod("getP","brob",function(x){x@positive})


###################################################
### code chunk number 24: setlength
###################################################



###################################################
### code chunk number 25: brob.Rnw:481-482
###################################################
setMethod("length","brob",function(x){length(x@x)})


###################################################
### code chunk number 26: setmethodSquareBrace
###################################################



###################################################
### code chunk number 27: brob.Rnw:492-499
###################################################
setMethod("[", "brob",
          function(x, i, j,  drop){
            if(!missing(j)){
              warning("second argument to extractor function ignored")
            }
            brob(x@x[i], x@positive[i])
          } )


###################################################
### code chunk number 28: setReplaceMethod
###################################################



###################################################
### code chunk number 29: brob.Rnw:512-529
###################################################
setReplaceMethod("[",signature(x="brob"),
                 function(x,i,j,value){
                   if(!missing(j)){
                     warning("second argument to extractor function ignored")
                   }
                   jj.x <- x@x
                   jj.pos <- x@positive
                   if(is.brob(value)){
                     jj.x[i] <- value@x
                     jj.pos[i] <- value@positive
                     return(brob(x=jj.x,positive=jj.pos))
                   } else {
                     x[i] <- as.brob(value)
                     return(x)
                   }
                 } )



###################################################
### code chunk number 30: .Brob.cPair
###################################################
.Brob.cPair <- function(x,y){
  x <- as.brob(x)
  y <- as.brob(y)
  brob(c(x@x,y@x),c(x@positive,y@positive))
}


###################################################
### code chunk number 31: setGeneric_cbrob
###################################################



###################################################
### code chunk number 32: brob.Rnw:568-569
###################################################
setGeneric(".cPair", function(x,y){standardGeneric(".cPair")})


###################################################
### code chunk number 33: setMethod.Cpair
###################################################



###################################################
### code chunk number 34: brob.Rnw:579-583
###################################################
setMethod(".cPair", c("brob", "brob"), function(x,y){.Brob.cPair(x,y)})
setMethod(".cPair", c("brob", "ANY"),  function(x,y){.Brob.cPair(x,as.brob(y))})
setMethod(".cPair", c("ANY", "brob"),  function(x,y){.Brob.cPair(as.brob(x),y)})
setMethod(".cPair", c("ANY", "ANY"),   function(x,y){c(x,y)})


###################################################
### code chunk number 35: cbrob
###################################################
"cbrob" <- function(x, ...) {
   if(nargs()<3)
      .cPair(x,...)
    else
      .cPair(x, Recall(...))
}


###################################################
### code chunk number 36: test.cbrob
###################################################
a <- 1:3
b <- as.brob(1e100)
cbrob(a,a,b,a)


###################################################
### code chunk number 37: sqrtmethod
###################################################



###################################################
### code chunk number 38: brob.Rnw:633-636
###################################################
setMethod("sqrt","brob", function(x){
 brob(ifelse(x@positive,x@x/2, NaN),TRUE)
} )


###################################################
### code chunk number 39: checklogsqrt
###################################################
sqrt(brob(4))


###################################################
### code chunk number 40: mathgeneric
###################################################



###################################################
### code chunk number 41: brob.Rnw:648-679
###################################################
setMethod("Math", "brob",
          function(x){
            switch(.Generic,
                   abs    = brob(x@x),
                   log    = {
                     out <- x@x
                     out[!x@positive] <- NaN
                     out
                   },
                   exp    = brob(x),
                   cosh   = {(brob(x) + brob(-x))/2},
                   sinh   = {(brob(x) - brob(-x))/2},
                   acos   =,
                   acosh  =,
                   asin   =,
                   asinh  =,
                   atan   =,
                   atanh  =,
                   cos    =,
                   sin    =,
                   tan    =,
                   tanh   =,
                   trunc  = callGeneric(as.numeric(x)),
                   lgamma =,
                   cumsum =,
                   gamma  =,
                   ceiling=,
                   floor  = as.brob(callGeneric(as.numeric(x))),
                   stop(paste(.Generic, "not allowed on Brobdingnagian numbers"))
                     )
          } )


###################################################
### code chunk number 42: checktrig
###################################################
sin(brob(4))


###################################################
### code chunk number 43: .brob.arithstuff
###################################################
.Brob.negative <- function(e1){
  brob(e1@x,!e1@positive)
}
.Brob.ds <- function(e1,e2){   # "ds" == "different signs"
  xor(e1@positive,e2@positive)
}

.Brob.add <- function(e1,e2){
  e1 <- as.brob(e1)
  e2 <- as.brob(e2)
  
  jj <- rbind(e1@x,e2@x)
  x1 <- jj[1,]
  x2 <- jj[2,]
  out.x <- double(length(x1))
  
  jj <- rbind(e1@positive,e2@positive)
  p1 <- jj[1,]
  p2 <- jj[2,]
  out.pos <- p1
  
  ds <- .Brob.ds(e1,e2)
  ss <- !ds                             #ss == "Same Sign"

  out.x[ss] <- pmax(x1[ss],x2[ss]) + log1p(+exp(-abs(x1[ss]-x2[ss])))
  out.x[ds] <- pmax(x1[ds],x2[ds]) + log1p(-exp(-abs(x1[ds]-x2[ds])))

  # Now special dispensation for 0+0:
  out.x[ (x1 == -Inf) & (x2 == -Inf)] <- -Inf
  out.pos <- p1
  out.pos[ds] <- xor((x1[ds] > x2[ds]) , (!p1[ds]) )
  return(brob(out.x,out.pos))
}

.Brob.mult <- function(e1,e2){
  e1 <- as.brob(e1)
  e2 <- as.brob(e2)
  return(brob(e1@x + e2@x, !.Brob.ds(e1,e2)))
}

.Brob.power <- function(e1,e2){
  stopifnot(is.brob(e1) | is.brob(e2))
  if(is.brob(e2)){ # e2 a brob => answer a brob (ignore signs)
    return(brob(log(e1) * brob(e2@x), TRUE))
  } else {  #e2 a non-brob (try to account for signs)
    s <- as.integer(2*e1@positive-1) #s = +/-1
    return(brob(e1@x*as.brob(e2),  (s^as.numeric(e2))>0))
  }
}

.Brob.inverse <- function(b){brob(-b@x,b@positive)}


###################################################
### code chunk number 44: setMethodArithUnary
###################################################



###################################################
### code chunk number 45: brob.Rnw:771-780
###################################################
setMethod("Arith",signature(e1 = "brob", e2="missing"),
          function(e1,e2){
            switch(.Generic,
                   "+" = e1,
                   "-" = .Brob.negative(e1),
                   stop(paste("Unary operator", .Generic,
                              "not allowed on Brobdingnagian numbers"))
                   )
          } )


###################################################
### code chunk number 46: check_minus_5
###################################################
-brob(5)


###################################################
### code chunk number 47: brob.arith
###################################################
.Brob.arith <- function(e1,e2){
  switch(.Generic,
         "+" = .Brob.add  (e1, e2),
         "-" = .Brob.add  (e1, .Brob.negative(as.brob(e2))),
         "*" = .Brob.mult (e1, e2),
         "/" = .Brob.mult (e1, .Brob.inverse(as.brob(e2))),
         "^" = .Brob.power(e1, e2),
         stop(paste("binary operator \"", .Generic, "\" not defined for Brobdingnagian numbers"))
         ) }


###################################################
### code chunk number 48: setMethodArith
###################################################
setMethod("Arith", signature(e1 = "brob", e2="ANY"), .Brob.arith)
setMethod("Arith", signature(e1 = "ANY", e2="brob"), .Brob.arith)
setMethod("Arith", signature(e1 = "brob", e2="brob"), .Brob.arith)


###################################################
### code chunk number 49: check_addition
###################################################
1e100 + as.brob(10)^100


###################################################
### code chunk number 50: brob.equalandgreater
###################################################
.Brob.equal <- function(e1,e2){
  (e1@x==e2@x) & (e1@positive==e2@positive)
}

.Brob.greater <- function(e1,e2){
  jj.x <- rbind(e1@x,e2@x)
  jj.p <- rbind(e1@positive,e2@positive)

  ds <- .Brob.ds(e1,e2)
  ss <- !ds                             #ss == "Same Sign"
  greater <- logical(length(ss))
  
  greater[ds] <- jj.p[1,ds]
  greater[ss] <- jj.p[1,ss] & (jj.x[1,ss] > jj.x[2,ss])
  return(greater)
}


###################################################
### code chunk number 51: brob.compare
###################################################
".Brob.compare" <- function(e1,e2){
  e1 <- as.brob(e1)
  e2 <- as.brob(e2)
  switch(.Generic,
         "==" =  .Brob.equal(e1,e2),
         "!=" = !.Brob.equal(e1,e2),
         ">"  =  .Brob.greater(e1,e2),
         "<"  = !.Brob.greater(e1,e2) & !.Brob.equal(e1,e2),
         ">=" =  .Brob.greater(e1,e2) |  .Brob.equal(e1,e2),
         "<=" = !.Brob.greater(e1,e2) |  .Brob.equal(e1,e2),
         stop(paste(.Generic, "not supported for Brobdingnagian numbers"))
         )
}


###################################################
### code chunk number 52: setMethodCompare
###################################################



###################################################
### code chunk number 53: brob.Rnw:876-879
###################################################
setMethod("Compare", signature(e1="brob", e2="ANY" ), .Brob.compare)
setMethod("Compare", signature(e1="ANY" , e2="brob"), .Brob.compare)
setMethod("Compare", signature(e1="brob", e2="brob"), .Brob.compare)


###################################################
### code chunk number 54: check.compare
###################################################
as.brob(10) < as.brob(11)
as.brob(10) <= as.brob(10)


###################################################
### code chunk number 55: brob.logic
###################################################
.Brob.logic <- function(e1,e2){
  stop("No logic currently implemented for Brobdingnagian numbers")
}


###################################################
### code chunk number 56: setmethodlogic
###################################################



###################################################
### code chunk number 57: brob.Rnw:910-913
###################################################
setMethod("Logic",signature(e1="swift",e2="ANY"), .Brob.logic)
setMethod("Logic",signature(e1="ANY",e2="swift"), .Brob.logic)
setMethod("Logic",signature(e1="swift",e2="swift"), .Brob.logic)


###################################################
### code chunk number 58: logchunk
###################################################
if(!isGeneric("log")){
  setGeneric("log",group="Math")
}


###################################################
### code chunk number 59: miscgenerics
###################################################



###################################################
### code chunk number 60: brob.Rnw:958-1009
###################################################
if(!isGeneric("sum")){
setGeneric("max", function(x, ..., na.rm = FALSE)
	{
		standardGeneric("max")
	},
	useAsDefault = function(x, ..., na.rm = FALSE)
	{
		base::max(x, ..., na.rm = na.rm)
	},
	group = "Summary")

setGeneric("min", function(x, ..., na.rm = FALSE)
	{
		standardGeneric("min")
	},
	useAsDefault = function(x, ..., na.rm = FALSE)
	{
		base::min(x, ..., na.rm = na.rm)
	},
	group = "Summary")

setGeneric("range", function(x, ..., na.rm = FALSE)
	{
		standardGeneric("range")
	},
	useAsDefault = function(x, ..., na.rm = FALSE)
	{
		base::range(x, ..., na.rm = na.rm)
	},
	group = "Summary")

setGeneric("prod", function(x, ..., na.rm = FALSE)
	{
		standardGeneric("prod")
	},
	useAsDefault = function(x, ..., na.rm = FALSE)
	{
		base::prod(x, ..., na.rm = na.rm)
	},
	group = "Summary")

setGeneric("sum", function(x, ..., na.rm = FALSE)
	{
		standardGeneric("sum")
	},
	useAsDefault = function(x, ..., na.rm = FALSE)
	{
		base::sum(x, ..., na.rm = na.rm)
	},
	group = "Summary")
}


###################################################
### code chunk number 61: brob.maxmin
###################################################
.Brob.max <- function(x, ..., na.rm=FALSE){
  p <- x@positive
  val <- x@x
  if(any(p)){
    return(brob(max(val[p])))
  } else {
    # all negative
    return(brob(min(val),FALSE))
  }
}

.Brob.prod <- function(x){
  p <- x@positive
  val <- x@x
  return(brob(sum(val),(sum(p)%%2)==0))
}

.Brob.sum <- function(x){
  .Brob.sum.allpositive( x[x>0]) -
  .Brob.sum.allpositive(-x[x<0]) 
}

.Brob.sum.allpositive <- function(x){
  if(length(x)<1){return(as.brob(0))}
  val <- x@x
  p <- x@positive
  mv <- max(val)
  return(brob(mv + log1p(sum(exp(val[-which.max(val)]-mv))),TRUE))
}


###################################################
### code chunk number 62: setmethodsummary
###################################################



###################################################
### code chunk number 63: brob.Rnw:1055-1067
###################################################
setMethod("Summary", "brob",
          function(x, ..., na.rm=FALSE){
            switch(.Generic,
                   max    =  .Brob.max( x, ..., na.rm=na.rm),
                   min    = -.Brob.max(-x, ..., na.rm=na.rm),
                   range  =   cbrob(min(x,na.rm=na.rm),max(x,na.rm=na.rm)),
                   prod   =  .Brob.prod(x),
                   sum    =  .Brob.sum(x),
                   stop(paste(.Generic, "not allowed on Brobdingnagian numbers"))
                   )
          }
          )


###################################################
### code chunk number 64: checksum
###################################################
sum(as.brob(1:100)) - 5050


###################################################
### code chunk number 65: factorial
###################################################
stirling <- function(x){sqrt(2*pi*x)*exp(-x)*x^x}


###################################################
### code chunk number 66: use.stirling
###################################################
stirling(100)
stirling(as.brob(100))


###################################################
### code chunk number 67: compare.two.stirlings
###################################################
as.numeric(stirling(100)/stirling(as.brob(100)))


###################################################
### code chunk number 68: stirling.of.1000
###################################################
stirling(1000)
stirling(as.brob(1000))


