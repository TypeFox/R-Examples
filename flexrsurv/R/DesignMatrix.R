# Methods for classes DesignMatrix, DesignMatrixNPH, DesignMatrixNPHNLL defined in AllClass.R 

DesignMatrixNPHNLL <- function(Z, listsplinebasis, timesplinebasis){
#  Z : matrix of raw variables
#  listsplinebasis : list of spline parameters for each Zi.
#  each splineparameter is S4 class AnySplineBasis (SplinBasis or descendant)
  
  nZ <- dim(Z)[2]

#  cat("DesignMatrixNPHNLL \n")
#  cat(nZ, dim(Z))
#  cat(dimnames(Z)[[2]])
#  cat("\nDesignMatrixNPHNLL2 \n")
  if(length(listsplinebasis)!= nZ){
    stop("Incorrect number od Spline Parameters")
  }
  index <- matrix(1, nrow=nZ, ncol=2)

  if(is.null(dimnames(Z)[[2]])) {
    namesZ <- paste("VNPHNLL", 1:nZ, sep="")
  } else {
    namesZ <- dimnames(Z)[[2]] 
  }

  
  ZZ <- NULL
  for (i in 1:nZ){
    ZZi <- evaluate(listsplinebasis[[i]], Z[,i], intercept=FALSE, xname=namesZ[i])
#    cat(i, dim(ZZi), namesZ[i], dimnames(ZZi)[[2]])
    if(i>1) {
      index[i,1] <- index[i-1, 2]+1
      index[i,2] <- index[i-1, 2]+ dim(ZZi)[2] 
    }
    else {
      index[i,2] <- dim(ZZi)[2] 
    }
    ZZ <- cbind(ZZ, ZZi)
  }
  nparam <- dim(ZZ)[2]

#  print(index)
  
  signature <- matrix(0, ncol=nZ, nrow=nparam)
  for (i in 1:nZ){
    signature[index[i,1]:index[i,2],i] <- 1
  }

     
  
  new("DesignMatrixNPHNLL",
      DM=ZZ,
      nObs=dim(Z)[1],
      nZ=nZ,
      nparam=nparam,
      signature=signature,
      index=index,
      listSplineBasis=listsplinebasis,
      names=namesZ,
      TSplineBasis=timesplinebasis,
      nTbasis=getNBases(timesplinebasis)

      )
}

# extractores, Getteurs 

setGeneric("getDesignMatrix",function(object)standardGeneric("getDesignMatrix"))
setMethod("getDesignMatrix",signature("DesignMatrix"),function(object)object@DM)
#setMethod("getDesignMatrix",signature("DesignMatrixNPHNLL"),function(object)object@DM))

setGeneric("getNparam",function(object)standardGeneric("getNparam"))
setMethod("getNparam",signature("DesignMatrix"),function(object)dim(object@DM)[2])
setMethod("getNparam",signature("DesignMatrixNPHNLL"),function(object) object@nparam)

setGeneric("getNvar",function(object)standardGeneric("getNvar"))
setMethod("getNvar",signature("DesignMatrix"),function(object)dim(object@DM)[2])
setMethod("getNvar",signature("DesignMatrixNPHNLL"),function(object)object@nZ)

setGeneric("getNobs",function(object)standardGeneric("getNobs"))
setMethod("getNobs",signature("DesignMatrix"),function(object)dim(object@DM)[1])
setMethod("getNobs",signature("DesignMatrixNPHNLL"),function(object)object@nObs)


setGeneric("getNames",function(object)standardGeneric("getNames"))
setMethod("getNames",signature("DesignMatrix"),function(object)dim(names(object@DM)[[2]]))
setMethod("getNames",signature("DesignMatrixNPHNLL"),function(object) object@names)


setGeneric("getSignature",function(object)standardGeneric("getSignature"))
setMethod("getSignature",signature("DesignMatrixNPHNLL"),function(object) object@signature)

setGeneric("getIndex",function(object)standardGeneric("getIndex"))
setMethod("getIndex",signature("DesignMatrixNPHNLL"),function(object) object@index)

