setClass("glub",
         representation = representation(real="brob",imag="brob"),
         prototype      = list(real=new("brob"), imag=new("brob")),
         contains       = "swift"
         )

".Glub.valid" <- function(object){
  if(length(object@real) == length(object@imag)){
    return(TRUE)
  } else {
    return("length mismatch")
  }
}

setValidity("glub", .Glub.valid)

setAs("glub", "complex", function(from){
  complex(real=as.numeric(from@real), imaginary=as.numeric(from@imag))
} )

setMethod("as.complex",signature(x="glub"),function(x){as(x,"complex")})

setAs("glub", "numeric", function(from){
    warning("imaginary parts discarded in coercion; use as.complex() to retain them")
    as.numeric(Re(from))
} )
  
setMethod("as.numeric",signature(x="glub"),function(x){as(x,"numeric")})


"glub" <- function(real=double(), imag=double()){
  if(missing(imag)){
    imag <- 0
  }
  real <- as.brob(real)
  imag <- as.brob(imag)
  
  jj.x <- cbind(real@x,imag@x)
  jj.p <- cbind(real@positive,imag@positive)
  new("glub",
      real = brob(jj.x[,1],jj.p[,1]),
      imag = brob(jj.x[,2],jj.p[,2])
      )
}


setMethod("Re","glub",function(z){z@real})
setMethod("Im","glub",function(z){z@imag})

setMethod("length","glub",function(x){length(Re(x))})


setMethod("Mod", "glub", function(z){sqrt(Re(z)*Re(z) + Im(z)*Im(z))})

".Brob.arg" <- function(z){
  atan2(as.numeric(Im(z)),as.numeric(Re(z)))
}

".Glub.complex" <- function(z){
  switch(.Generic,
         Arg  = .Brob.arg(z),
         Conj = glub(Re(z),-Im(z)),
         stop(paste("Complex operator \"", .Generic, "\" not defined for Glub numbers"))
         )
}

setMethod("Complex","glub", .Glub.complex)

setGeneric("Re<-",function(z,value){standardGeneric("Re<-")})
setGeneric("Im<-",function(z,value){standardGeneric("Im<-")})

setMethod("Re<-","glub",function(z,value){
  return(glub(real=value, imag=Im(z)))
} )

setMethod("Im<-","glub",function(z,value){
  z <- as.glub(z)
  return(glub(real=z@real, imag=value))
} )

setMethod("Im<-","brob",function(z,value){
  return(glub(real=z, imag=value))
} )

"as.glub" <- function(x){
  if(is.glub(x)){
    return(x)
  } else if (is.brob(x)) {
    return(glub(real=as.brob(x),imag=as.brob(0)))
  } else {
    return(glub(real=as.brob(Re(x)),imag=as.brob(Im(x))))
  }
}

setMethod("[", "glub",
          function(x, i, j,  drop){
            if(!missing(j)){warning("second argument (j) ignored")}
            glub(x@real[i], x@imag[i])
          }
          )

setReplaceMethod("[",signature(x="glub"),
                 function(x,i,j,value){
                   if(!missing(j)){warning("second argument (j) ignored")}
                   value <- as.glub(value) 
                   x@real[i] <- Re(value)
                   x@imag[i] <- Im(value)
                   return(x)
                 }
                 )

setMethod(".cPair", c("glub", "glub"), function(x,y).Glub.cPair(x,y))
setMethod(".cPair", c("glub", "ANY"),  function(x,y).Glub.cPair(x,as.glub(y)))
setMethod(".cPair", c("ANY", "glub"),  function(x,y).Glub.cPair(as.glub(x),y))
setMethod(".cPair", c("complex", "brob"),  function(x,y).Glub.cPair(as.glub(x),y))
setMethod(".cPair", c("brob", "complex"),  function(x,y).Glub.cPair(as.glub(x),y))
setMethod(".cPair", c("glub", "brob"),  function(x,y).Glub.cPair(as.glub(x),y))
setMethod(".cPair", c("brob", "glub"),  function(x,y).Glub.cPair(as.glub(x),y))

".Glub.cPair" <- function(x,y){
  x <- as.glub(x)
  y <- as.glub(y)
  return(glub(.Brob.cPair(Re(x),Re(y)), .Brob.cPair(Im(x),Im(y))))
}

"print.glub" <- function(x,...){
  real <- .Brob.print(Re(x),...)
  imag <- .Brob.print(Im(x),...)
  jj <- noquote(paste(real,imag,"i  ",sep=""))
  print(jj)
}

setMethod("show", "glub", function(object){print.glub(object)})

setMethod("Math", "glub",
          function(x){
            switch(.Generic,
                   abs    = Mod(x),
                   log    = { glub(log(Mod(x)),Arg(x)) },
                   log10  = { glub(log10(Mod(x)),Arg(x)/log(10)) },
                   log2   = { glub(log2 (Mod(x)),Arg(x)/log( 2)) },
                   exp    = { exp(Re(x))*exp(1i*as.numeric(Im(x)))},
                   sqrt   = { exp(log(x)/2)},
                   cosh   = { (exp(x)+exp(-x))/2},
                   sinh   = { (exp(x)-exp(-x))/2},
                   tanh   = { (exp(x)-exp(-x))/(exp(x)+exp(-x))},
                   cos    = { (exp(1i*x)+exp(-1i*x))/(2 )},
                   sin    = { (exp(1i*x)-exp(-1i*x))/(2i)},
                   tan    = { (exp(1i*x)-exp(-1i*x))/(exp(1i*x)+exp(-1i*x))},
                   acos   = { -1i*log(   x + 1i*sqrt( 1-x*x)) },
                   acosh  = {     log(   x +    sqrt(-1+x*x)) },
                   asin   = { -1i*log(1i*x +    sqrt( 1-x*x)) },
                   asinh  = {     log(   x +    sqrt( 1+x*x)) },
                   atan   = { 0.5i*log((1i+x)/(1i-x)) },
                   atanh  = { 0.5 *log((1 +x)/(1 -x)) },
                   trunc  = callGeneric(as.complex(x)),
                   lgamma =,
                   cumsum =,
                   gamma  =,
                   ceiling=,
                   floor  = as.glub(callGeneric(as.complex(x))),
                   stop(paste(.Generic, "not allowed on Brobdingnagian numbers"))
                     )
          }
)

".Glub.negative" <- function(e1){
  glub(-Re(e1),-Im(e1))
}

".Glub.add" <- function(e1,e2){
  e1 <- as.glub(e1)
  e2 <- as.glub(e2)

  glub(Re(e1)+Re(e2),Im(e1)+Im(e2))
}

".Glub.mult" <- function(e1,e2){
  e1 <- as.glub(e1)
  e2 <- as.glub(e2)

glub(Re(e1)*Re(e2)-Im(e1)*Im(e2), Re(e1)*Im(e2)+Im(e1)*Re(e2))
}

".Glub.power" <- function(e1,e2){
  exp(e2*log(e1))
}

".Glub.inverse" <- function(e1){
  jj <- Re(e1)*Re(e1) + Im(e1)*Im(e1)
  glub(Re(e1)/jj, -Im(e1)/jj)
}


".Glub.arith" <- function(e1,e2){
  switch(.Generic,
         "+" = .Glub.add  (e1, e2),
         "-" = .Glub.add  (e1, .Glub.negative(e2)),
         "*" = .Glub.mult (e1, e2),
         "/" = .Glub.mult (e1, .Glub.inverse(e2)),
         "^" = .Glub.power(e1, e2),
         stop(paste("binary operator \"", .Generic, "\" not defined for Glub numbers"))
         )
  }


setMethod("Arith",signature(e1 = "glub", e2="missing"),
          function(e1,e2){
            switch(.Generic,
                   "+" = e1,
                   "-" = .Glub.negative(e1),
                   stop(paste("Unary operator", .Generic,
                              "not allowed on Brobdingnagian numbers"))
                   )
          }
          )

setMethod("Arith",signature(e1 = "glub", e2="glub"), .Glub.arith)
setMethod("Arith",signature(e1 = "glub", e2="ANY" ), .Glub.arith)
setMethod("Arith",signature(e1 = "ANY" , e2="glub"), .Glub.arith)
setMethod("Arith",signature(e1= "brob"   , e2="complex"), .Glub.arith)
setMethod("Arith",signature(e1= "complex", e2="brob"   ), .Glub.arith)
setMethod("Arith",signature(e1= "glub"   , e2="complex"), .Glub.arith)
setMethod("Arith",signature(e1= "complex", e2="glub"   ), .Glub.arith)
setMethod("Arith",signature(e1= "glub", e2="brob"), .Glub.arith)
setMethod("Arith",signature(e1= "brob", e2="glub"), .Glub.arith)



".Glub.equal" <- function(e1,e2){
  (Re(e1) == Re(e2)) & ( Im(e1) == Im(e2))
}

".Glub.compare" <- function(e1,e2){
  e1 <- as.glub(e1)
  e2 <- as.glub(e2)
  switch(.Generic,
         "==" =  .Glub.equal(e1,e2),
         "!=" = !.Glub.equal(e1,e2),
         stop(paste(.Generic, "not supported for Glub numbers"))
         )
}

setMethod("Compare", signature(e1="glub",e2="glub"), .Glub.compare)
setMethod("Compare", signature(e1="glub",e2="ANY" ), .Glub.compare)
setMethod("Compare", signature(e1="ANY", e2="glub"), .Glub.compare)

setMethod("Compare", signature(e1="brob", e2="glub"), .Glub.compare)
setMethod("Compare", signature(e1="glub", e2="brob"), .Glub.compare)


".Glub.prod" <- function(z){
  out <- as.glub(1)
  for(i in 1:length(z)){
    out <- out * z[i]
  }
  return(out)
}

".Glub.sum" <- function(x){
  glub(sum(Re(x)),sum(Im(x)))
}

setMethod("Summary", "glub",
          function(x, ..., na.rm=FALSE){
            switch(.Generic,
                   prod   =  .Glub.prod(x),
                   sum    =  .Glub.sum(x),
                   stop(paste('\"', .Generic, '()\" not allowed on Glubbdubdribbian numbers',sep=""))
                   )
          }
          )
setMethod("plot",signature(x="glub",y="missing"),function(x, ...){plot.default(as.complex(x), ...)})
setMethod("plot",signature(x="glub",y="ANY" ),function(x, y, ...){plot.default(as.complex(x), as.complex(y), ...)})
setMethod("plot",signature(x="ANY" ,y="glub"),function(x, y, ...){plot.default(as.complex(x), as.complex(y), ...)})
