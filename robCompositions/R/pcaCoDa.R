#' Robust principal component analysis for compositional data
#' 
#' This function applies robust principal component analysis for compositional
#' data.
#' 
#' The compositional data set is transformed using the ilr tranformation.
#' Afterwards, robust principal component analysis is performed.  Resulting
#' loadings and scores are back-transformed to the clr space where the
#' compositional biplot can be shown.
#' 
#' \code{mult_comp} is used when there are more than one group of compositional
#' parts in the data. To give an illustrative example, lets assume that one
#' variable group measures angles of the inner ear-bones of animals which sum
#' up to 100 and another one having percentages of a whole on the thickness of
#' the inner ear-bones included. Then two groups of variables exists which are
#' both compositional parts. The ilr-transformation is then internally applied
#' to each group independently whenever the \code{mult_comp} is set correctly.
#' 
#' @aliases pcaCoDa print.pcaCoDa
#' @param x compositional data
#' @param method either \dQuote{robust} (default) or \dQuote{standard}
#' @param mult_comp a list of numeric vectors holding the indices of linked
#' compositions
#' @return \item{scores }{scores in clr space} \item{loadings }{loadings in clr
#' space} \item{eigenvalues }{eigenvalues of the clr covariance matrix}
#' \item{method }{method} \item{princompOutputClr }{output of \code{princomp}
#' needed in \code{plot.pcaCoDa}}
#' @author K. Hron, P. Filzmoser, M. Templ
#' @importFrom stats princomp
#' @seealso \code{\link{print.pcaCoDa}}, \code{\link{plot.pcaCoDa}}
#' @references Filzmoser, P., Hron, K., Reimann, C. (2009) Principal Component
#' Analysis for Compositional Data with Outliers. \emph{Environmetrics},
#' \bold{20}, 621-632.
#' @keywords multivariate
#' @export
#' @importFrom MASS cov.mve
#' @examples
#' 
#' data(expenditures)
#' p1 <- pcaCoDa(expenditures)
#' p1
#' plot(p1)
#' 
#' ## just for illustration how to set the mult_comp argument
#' p1 <- pcaCoDa(expenditures, mult_comp=list(c(1,2,3),c(4,5)))
#' p1
pcaCoDa <- function(x, method="robust",mult_comp=NULL){
  
  # Closure problem with ilr transformation
  ilrV <- function(x){
    # ilr transformation
    x.ilr=matrix(NA,nrow=nrow(x),ncol=ncol(x)-1)
    for (i in 1:ncol(x.ilr)){
      x.ilr[,i]=sqrt((i)/(i+1))*log(((apply(as.matrix(x[,1:i]), 1, prod))^(1/i))/(x[,i+1]))
    }
    return(x.ilr)
  }
  if(is.null(mult_comp)){
    xilr <- ilrV(x)
  }else{
    xilr <- do.call("cbind",lapply(mult_comp,function(xx)ilrV(x[,xx])))
  }		
  if( method == "robust"){
    cv <- robustbase::covMcd(xilr, cor=FALSE)
    pcaIlr <- suppressWarnings(princomp(xilr, covmat=cv, cor=FALSE))
    eigenvalues <- eigen(cv$cov)$values
  } else if (method =="mve"){
    cv <- MASS::cov.mve(xilr)
    pcaIlr <- suppressWarnings(princomp(xilr, covmat=cv, cor=FALSE))
    eigenvalues <- eigen(cv$cov)$values
  } else {
    pcaIlr <- princomp(xilr, cor=FALSE)
    eigenvalues <- eigen(cov(xilr))$values
  }
  # construct orthonormal basis
  if(is.null(mult_comp)){
    V <- matrix(0, nrow=ncol(x), ncol=ncol(x)-1)
    for( i in 1:ncol(V) ){
      V[1:i,i] <- 1/i
      V[i+1,i] <- (-1)
      V[,i] <- V[,i]*sqrt(i/(i+1))
    }
  }else{
    V <- matrix(0, nrow=length(unlist(mult_comp)), ncol=length(unlist(mult_comp))-length(mult_comp))
    l <- sapply(mult_comp,length)
    start <- c(1,cumsum(l[-length(l)]))
    cumsum(l[-length(l)])
    end <- cumsum(l-1)
    start2 <- c(1,cumsum(l[-length(l)])+1)
    end2 <- cumsum(l)
    for(j in 1:length(mult_comp)){
      ind <- start[j]:end[j]
      ind2 <- start2[j]:end2[j]
      for( i in 1:length(ind)){
        V[ind2[1:i],ind[i]] <- 1/i
        V[ind2[i]+1,ind[i]] <- (-1)
        V[,ind[i]] <- V[,ind[i]]*sqrt(i/(i+1))
      }
    }
  }
  
 
  
  
  # robust ilr result - back-transformed to clr-space
  
  loadings <- V %*% pcaIlr$loadings	
  if(is.null(mult_comp)){
    if(!is.null(names(x))) dimnames(loadings)[[1]] <- names(x)
  }else{
    if(!is.null(names(x))) dimnames(loadings)[[1]] <- colnames(x)[unlist(mult_comp)]
  }
  pcaClr <- pcaIlr
#	pcaClr$scores <- pcaIlr$scores %*% t(V)
  pcaClr$scores <- pcaIlr$scores 
  pcaClr$loadings <- loadings
  
  res <- list(scores = pcaClr$scores,
      loadings = loadings,
      eigenvalues = eigenvalues,
      method = method,
      princompOutputClr = pcaClr,
      mult_comp = mult_comp
  )
  class(res) <- "pcaCoDa"
  invisible(res)
}

#' @rdname pcaCoDa
#' @export
#' @method print pcaCoDa
#' @param ... additional parameters for print method passed through
print.pcaCoDa <- function(x, ...){
  ## percentage of explained variability for clr transformed data
  eV <- x$eigenvalues / sum(x$eigenvalues)
  eVcum <- cumsum(x$eigenvalues) / sum(x$eigenvalues)
  cat("\n-------------------")
  cat("\n Percentages of explained variability for compositional data after clr transformation \n")
  print(eVcum)
  cat("\n-------------------\n\n")	
}
