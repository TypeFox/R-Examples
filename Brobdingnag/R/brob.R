setClass("swift",
         representation = "VIRTUAL"
         )

setClass("brob",
         representation = representation(x="numeric",positive="logical"),
         prototype      = list(x=numeric(),positive=logical()),
         contains       = "swift"
         )

".Brob.valid" <- function(object){
  len <- length(object@positive)
  if(len != length(object@x)){
    return("length mismatch")
  } else {
    return(TRUE)
  }
}
setValidity("brob", .Brob.valid)


"brob" <- function(x=double(),positive){
  if(missing(positive)){
    positive <- rep(TRUE,length(x))
  }
  if(length(positive)==1){
    positive <- rep(positive,length(x))
  }
  new("brob",x=as.numeric(x),positive=positive)
}

"is.brob" <- function(x){is(x,"brob")}
"is.glub" <- function(x){is(x,"glub")}

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

setAs("brob", "numeric", function(from){
  out <- exp(from@x)
  out[!from@positive] <- -out[!from@positive]
  return(out)
} )

setMethod("as.numeric",signature(x="brob"),function(x){as(x,"numeric")})

setAs("brob", "complex", function(from){
  return(as.numeric(from)+ 0i)
} )

setMethod("as.complex",signature(x="brob"),function(x){as(x,"complex")})

".Brob.print" <- function(x, digits=5){
     noquote( paste(c("-","+")[1+x@positive],"exp(",signif(x@x,digits),")",sep=""))
   }
    
"print.brob" <- function(x, ...){
  jj <- .Brob.print(x, ...)
  print(jj)
  return(invisible(jj))
}

setMethod("show", "brob", function(object){print.brob(object)})

setGeneric("getX",function(x){standardGeneric("getX")})
setGeneric("getP",function(x){standardGeneric("getP")})
setMethod("getX","brob",function(x){x@x})
setMethod("getP","brob",function(x){x@positive})
setMethod("length","brob",function(x){length(x@x)})

setGeneric("sign<-",function(x,value){standardGeneric("sign<-")})
setMethod("sign<-","brob",function(x,value){
  brob(x@x,value)
} )

setMethod("[", "brob",
          function(x, i, j,  drop){
            if(!missing(j)){
              warning("second argument to extractor function ignored")
            }
            brob(x@x[i], x@positive[i])
          } )

setReplaceMethod("[",signature(x="brob"),
                 function(x,i,j,value){
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

setGeneric(".cPair", function(x,y){standardGeneric(".cPair")})
setMethod(".cPair", c("brob", "brob"), function(x,y){.Brob.cPair(x,y)})
setMethod(".cPair", c("brob", "ANY"),  function(x,y){.Brob.cPair(x,as.brob(y))})
setMethod(".cPair", c("ANY", "brob"),  function(x,y){.Brob.cPair(as.brob(x),y)})
setMethod(".cPair", c("ANY", "ANY"),   function(x,y){c(x,y)})

"cbrob" <- function(x, ...) {
   if(nargs()<3)
      .cPair(x,...)
    else
      .cPair(x, Recall(...))
}

".Brob.cPair" <- function(x,y){
  x <- as.brob(x)
  y <- as.brob(y)
  brob(c(x@x,y@x),c(x@positive,y@positive))
}

setGeneric("log")

setMethod("sqrt","brob", function(x){
 brob(ifelse(x@positive,x@x/2, NaN),TRUE)
} )
          
setMethod("Math", "brob",
          function(x){
            switch(.Generic,
                   abs    = brob(x@x),
                   log    = {
                     out <- x@x
                     out[!x@positive] <- NaN
                     out
                   },
                   log10  = {
                     out <- x@x/log(10)
                     out[!x@positive] <- NaN
                     out
                   },
                   log2 = {
                     out <- x@x/log(2)
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

".Brob.negative" <- function(e1){
  brob(e1@x,!e1@positive)
}

".Brob.ds" <- function(e1,e2){   # "ds" == "different signs"
  xor(e1@positive,e2@positive)
}

".Brob.add" <- function(e1,e2){
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

".Brob.mult" <- function(e1,e2){
  e1 <- as.brob(e1)
  e2 <- as.brob(e2)
  return(brob(e1@x + e2@x, !.Brob.ds(e1,e2)))
}

".Brob.power"<- function(e1,e2){
  stopifnot(is.brob(e1) | is.brob(e2))
  if(is.brob(e2)){ # e2 a brob => answer a brob (ignore signs)
    return(brob(log(e1) * brob(e2@x), TRUE))
  } else {  #e2 a non-brob (try to account for signs)
    s <- as.integer(2*e1@positive-1) #s = +/-1
    return(brob(e1@x*as.brob(e2),  (s^as.numeric(e2))>0))
  }
}

".Brob.inverse" <- function(b){brob(-b@x,b@positive)}

setMethod("Arith",signature(e1 = "brob", e2="missing"),
          function(e1,e2){
            switch(.Generic,
                   "+" = e1,
                   "-" = .Brob.negative(e1),
                   stop(paste("Unary operator", .Generic,
                              "not allowed on Brobdingnagian numbers"))
                   )
          } )

".Brob.arith" <- function(e1,e2){
  switch(.Generic,
         "+" = .Brob.add  (e1, e2),
         "-" = .Brob.add  (e1, .Brob.negative(as.brob(e2))),
         "*" = .Brob.mult (e1, e2),
         "/" = .Brob.mult (e1, .Brob.inverse(as.brob(e2))),
         "^" = .Brob.power(e1, e2),
         stop(paste("binary operator \"", .Generic, "\" not defined for Brobdingnagian numbers"))
         ) }

setMethod("Arith", signature(e1 = "brob", e2="ANY"), .Brob.arith)
setMethod("Arith", signature(e1 = "ANY", e2="brob"), .Brob.arith)
setMethod("Arith", signature(e1 = "brob", e2="brob"), .Brob.arith)


".Brob.equal" <- function(e1,e2){
  (e1@x==e2@x) & (e1@positive==e2@positive)
}

".Brob.greater" <- function(e1,e2){
  jj.x <- rbind(e1@x,e2@x)
  jj.p <- rbind(e1@positive,e2@positive)

  ds <- .Brob.ds(e1,e2)
  ss <- !ds                             #ss == "Same Sign"
  greater <- logical(length(ss))
  
  greater[ds] <- jj.p[1,ds]
  greater[ss] <- jj.p[1,ss] & (jj.x[1,ss] > jj.x[2,ss])
  return(greater)
}

".Brob.compare" <- function(e1,e2){
   if( (length(e1) == 0) | (length(e2)==0)) {
       return(logical(0))
   }

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

setMethod("Compare", signature(e1="brob", e2="ANY" ), .Brob.compare)
setMethod("Compare", signature(e1="ANY" , e2="brob"), .Brob.compare)
setMethod("Compare", signature(e1="brob", e2="brob"), .Brob.compare)

".Brob.logic" <- function(e1,e2){
  stop("No logic currently implemented for Brobdingnagian numbers")
}

setMethod("Logic",signature(e1="swift",e2="ANY"), .Brob.logic)
setMethod("Logic",signature(e1="ANY",e2="swift"), .Brob.logic)
setMethod("Logic",signature(e1="swift",e2="swift"), .Brob.logic)

if(!isGeneric("max")){
setGeneric("max", function(x, ..., na.rm = FALSE)
	{
		standardGeneric("max")
	},
	useAsDefault = function(x, ..., na.rm = FALSE)
	{
		base::max(x, ..., na.rm = na.rm)
	},
	group = "Summary")
}

if(!isGeneric("min")){
setGeneric("min", function(x, ..., na.rm = FALSE)
	{
		standardGeneric("min")
	},
	useAsDefault = function(x, ..., na.rm = FALSE)
	{
		base::min(x, ..., na.rm = na.rm)
	},
	group = "Summary")
}

if(!isGeneric("range")){
setGeneric("range", function(x, ..., na.rm = FALSE)
	{
		standardGeneric("range")
	},
	useAsDefault = function(x, ..., na.rm = FALSE)
	{
		base::range(x, ..., na.rm = na.rm)
	},
	group = "Summary")
}

if(!isGeneric("prod")){
setGeneric("prod", function(x, ..., na.rm = FALSE)
	{
		standardGeneric("prod")
	},
	useAsDefault = function(x, ..., na.rm = FALSE)
	{
		base::prod(x, ..., na.rm = na.rm)
	},
	group = "Summary")
}

if(!isGeneric("sum")){
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

".Brob.max" <- function(x, ..., na.rm=FALSE){
  p <- x@positive
  val <- x@x
  if(any(p)){
    return(brob(max(val[p])))
  } else {
    # all negative
    return(brob(min(val),FALSE))
  }
}

".Brob.prod" <- function(x){
  p <- x@positive
  val <- x@x
  return(brob(sum(val),(sum(p)%%2)==0))
}

".Brob.sum" <- function(x){
  .Brob.sum.allpositive( x[x>0]) -
  .Brob.sum.allpositive(-x[x<0]) 
}

".Brob.sum.allpositive" <- function(x){
  if(length(x)<1){return(as.brob(0))}
  val <- x@x
  p <- x@positive
  mv <- max(val)
  return(brob(mv + log1p(sum(exp(val[-which.max(val)]-mv))),TRUE))
}

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


setMethod("plot",signature(x="brob",y="missing"),function(x, ...){plot.default(as.numeric(x), ...)})
setMethod("plot",signature(x="brob",y="ANY" ),function(x, y, ...){plot.default(as.numeric(x), as.numeric(y), ...)})
setMethod("plot",signature(x="ANY" ,y="brob"),function(x, y, ...){plot.default(as.numeric(x), as.numeric(y), ...)})


