# TODO: Add comment
#
# Author: iovleff
###############################################################################

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterDiagGaussian-method
setMethod(
  "missingValues",
  c("ClusterHeterogeneous"),
  function(x)
  {
    nbData <- length(x@ldata)
    res <- vector("list", nbData)
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        res[[l]]  <- cbind(x@ldata[[l]]@missing, (x@ldata[[l]]@data)[x@ldata[[l]]@missing]);
        colnames(res[[l]])[3] <- "value";
      }
    }
    return(res)
  }
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterDiagGaussianComponent-method
setMethod(
  "missingValues",
  c("ClusterDiagGaussianComponent"),
  function(x)
  { res = cbind(x@missing, x@data[x@missing]);
    colnames(res)[3] <- "value";
    return(res)
  }
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterDiagGaussian-method
setMethod(
  "missingValues",
  c("ClusterDiagGaussian"),
  function(x) { return(missingValues(x@component));}
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterGammaComponent-method
setMethod(
    "missingValues",
    c("ClusterGammaComponent"),
    function(x)
    { res = cbind(x@missing, x@data[x@missing]);
      colnames(res)[3] <- "value";
      return(res)
    }
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterGamma-method
setMethod(
  "missingValues",
  c("ClusterGamma"),
  function(x) { return(missingValues(x@component));}
)


#' @rdname missingValues-methods
#' @aliases missingValues,ClusterCategoricalComponent-method
setMethod(
  f="missingValues",
  signature=c("ClusterCategoricalComponent"),
  function(x)
  { res = cbind(x@missing, x@data[x@missing]);
    colnames(res)[3] <- "value";
    return(res)
  }
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterCategorical-method
setMethod(
    f="missingValues",
    signature=c("ClusterCategorical"),
    function(x) { return(missingValues(x@component));}
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterPoissonComponent-method
setMethod(
  "missingValues",
  c("ClusterPoissonComponent"),
  function(x)
  { res = cbind(x@missing, x@data[x@missing]);
    colnames(res)[3] <- "value";
    return(res)
  }
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterPoisson-method
setMethod(
  "missingValues",
  c("ClusterPoisson"),
  function(x) { return(missingValues(x@component));}
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterKernelComponent-method
setMethod(
    "missingValues",
    c("ClusterKernelComponent"),
    function(x)
    {
      return(NULL)
    }
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterKernel-method
setMethod(
    "missingValues",
    c("ClusterKernel"),
    function(x) { return(NULL);}
)
