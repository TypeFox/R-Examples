default_bands_MS<- function(object,  length.out=101){

  kn <- getKnots(object)
  bands <- c(kn[1])
  for( i in 2:length(kn)){
    bands <- c(bands, seq(kn[i-1], kn[i], length.out = length.out)[-1])
  }
  bands
}


default_bands_BS<- function(object,  length.out=101){

  kn <- getKnots(object)[getOrder(object):(length(getKnots(object))-getDegree(object))]
  bands <- c(kn[1])
  for( i in 2:length(kn)){
    bands <- c(bands, seq(kn[i-1], kn[i], length.out = length.out)[-1])
  }
  bands
}


default_bands_TPS<- function(object,  length.out=101){
#assume min(t)=0, max(t)=
  kn <- c(getMin(object), getKnots(object), getMax(object))
  bands <- c(0)
  for( i in 2:length(kn)){
    bands <- c(bands, seq(kn[i-1], kn[i], length.out = length.out)[-1])
  }
  bands
}


setGeneric("default_bands",function(object, length.out, ...)standardGeneric("default_bands"))
setMethod("default_bands",
          signature("MSplineBasis","numeric"),
          function(object, length.out) default_bands_MS(object=object, length.out=length.out))

setMethod("default_bands",
          signature("BSplineBasis","numeric"),
          function(object, length.out) default_bands_BS(object=object, length.out=length.out))

setMethod("default_bands",
          signature("TPSplineBasis","numeric"),
          function(object, length.out) default_bands_TPS(object=object, length.out=length.out))




