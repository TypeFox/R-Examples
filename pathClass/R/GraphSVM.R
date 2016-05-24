#' Implementation of a supervised classification framework introduced by Franck Rapaport et al., 2007.
#'
#' \code{mapping} must be a data.frame with at least two columns. The column names have to be \code{c('probesetID','graphID')}.
#' Where 'probesetID' is the probeset ID present in the expression matrix (i.e. \code{colnames(x)}) and 'graphID' is any ID that
#' represents the nodes in the diffusionKernel (i.e. \code{colnames(diffusionKernel)} or \code{rownames(diffusionKernel)}). The purpose of the this mapping is that
#' a gene or protein in the network might be represented by more than one probe set on the chip. Therefore, the algorithm must know
#' which genes/protein in the network belongs to which probeset on the chip.
#'
#' @param x a p x n matrix of expression measurements with p samples and n genes.
#' @param y a factor of length p comprising the class labels.
#' @param DEBUG should debugging information be plotted.
#' @param scale a character vector defining if the data should be centered and/or scaled.
#' Possible values are \emph{center} and/or \emph{scale}. Defaults to \code{c('center', 'scale')}.
#' @param Cs soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
#' @param stepsize amount of features that are discarded in each step of the feature elimination. Defaults to 10\%.
#' @param mapping a mapping that defines how probe sets are summarized to genes.
#' @param diffusionKernel the diffusion kernel which was pre-computed by using the function \code{\link{calc.diffusionKernel}}
#' @param useOrigMethod use the method originally decribed in the paper by Franck Rapaport et al. 2007
#' @return a graphSVM object \item{features}{the selected features} \item{error.bound}{the span bound of the model} \item{fit}{the fitted SVM model}
#' @references Rapaport F. et al. (2007). Classification of microarray data using gene networks. \emph{BMC Bioinformatics}
#' @export
#' @callGraphPrimitives
#' @note We combined the original method with a Recursive Feature Elimination in order to allow a feature selection.
#' The optimal number of features is found by using the span estimate. See Chapelle, O., Vapnik, V., Bousquet, O., and Mukherjee, S. (2002). Choosing multiple parameters for support vector machines. \emph{Machine Learning}, 46(1), 131-159.
#'
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(Biobase)
#' data(sample.ExpressionSet)
#' x <- t(exprs(sample.ExpressionSet))
#' y <- factor(pData(sample.ExpressionSet)$sex)
#' # create the mapping
#' library('hgu95av2.db')
#' mapped.probes <- mappedkeys(hgu95av2REFSEQ)
#' refseq <- as.list(hgu95av2REFSEQ[mapped.probes])
#' times <- sapply(refseq, length)
#' mapping <- data.frame(probesetID=rep(names(refseq), times=times), graphID=unlist(refseq), row.names=NULL, stringsAsFactors=FALSE)
#' mapping <- unique(mapping)
#' library(pathClass)
#' data(adjacency.matrix)
#' matched <- matchMatrices(x=x, adjacency=adjacency.matrix, mapping=mapping)
#' dk <- calc.diffusionKernel(L=matched$adjacency, is.adjacency=TRUE, beta=0) # beta should be tuned
#' res.gSVM <- crossval(matched$x, y, theta.fit=fit.graph.svm, folds=5, repeats=2, DEBUG=TRUE, parallel=FALSE, Cs=10^(-3:3), mapping=matched$mapping, diffusionKernel=dk)
#' }
fit.graph.svm <- function(x, y, DEBUG=FALSE, scale=c('center', 'scale'), Cs=10^c(-3:3), stepsize=0.1,  mapping, diffusionKernel, useOrigMethod=FALSE){
  best.bound = Inf
  feat = colnames(x)

  if(missing(diffusionKernel)) stop('Need the diffusion kernel for classification')
  if(missing(mapping)) stop('You must provide a mapping.')
  result <- NULL
  
  if(!useOrigMethod){

    while(NCOL(x) > 1){
      if(DEBUG) cat(NCOL(x),' Features left.\n')

      fit = calc.graph.svm(x=x, y=y, Cs=Cs, R=diffusionKernel, scale=scale, DEBUG=DEBUG)
      
      if(fit$error.bound <= best.bound){
        best.bound = fit$error.bound
        feat  = colnames(x)			
        best.fit   = fit
        if(DEBUG) cat('Model Updated. Spanbound=',best.bound,', C=',best.fit$C,', ', length(feat),'features.\n')
      }
      
      ord     <-  order(fit$w)
      remove  <- colnames(x)[ord[1:round(NCOL(x)*stepsize)]]
      keep    <- setdiff(colnames(x), remove)
      x       <- x[ , keep, drop=FALSE]
      diffusionKernel <- diffusionKernel[keep, keep, drop=FALSE]
    }
    if(DEBUG) cat('Best Model is: Spanbound=',best.fit$error.bound,', C=',best.fit$C,',', length(best.fit$features),'features.')
    
    result <- list(features=feat, error.bound = best.bound, fit=best.fit)
  } else{
    fit = calc.graph.svm(x=x, y=y, Cs=Cs, R=diffusionKernel, scale=scale, DEBUG=DEBUG, calcW=FALSE)
    result <- list(features=colnames(x), error.bound = fit$error.bound, fit=fit)
  }
  
  class(result) <- 'graphSVM'
  return(result)
}


# x: data
# y: vector (length n)
# Cs: regularization parameters
# R: diffusion kernel matrix between network node pairs
# scale: scale/center data?
# calcW:calc feature influences?
# DEBUG:print debug messages?
# orig : use original implementation
calc.graph.svm = function(x, y, Cs, R, scale, calcW=TRUE, DEBUG=FALSE){
  
  if(missing(x))     stop('No epxression matrix provided.')
  if(missing(y))     stop('No class-lables provided.')
  if(missing(Cs))    stop('No tuning parameter \'C\' provided.')
  if(missing(R))     stop('No kernel \'R\' provided.')
  if(missing(scale)) stop('Parameter \'scale\' must be in \'scale\', \'center\' or NULL.')
  if(length(levels(factor(y))) != 2) stop('y must have 2 levels.')

  scale.mean <- scale.std <- NULL

  if(!is.null(scale)){
    scale <- tolower(scale)

    if("center" %in% scale){
      x = scale(x,center=T)
      ## save centering coefficient
      scale.mean = attr(x,"scaled:center")
      names(scale.mean) = colnames(x)
    }
    
    if("scale" %in% scale){
      x = scale(x,center=F,scale=T)
      ## save scaling coefficient    
      scale.std = attr(x,"scaled:scale")
      names(scale.std) = colnames(x)
    }
  }

  ## unknowns = setdiff(colnames(x), PPIgenes)
  K = matrix(0, ncol=NROW(x), nrow=NROW(x))	
  ## if(length(PPIgenes) > 0){				
  ##   KPPI = x[,PPIgenes,drop=FALSE]%*%R%*%t(x[,PPIgenes,drop=FALSE])
  ##   K = K + KPPI
  ## }

  KPPI = x %*% R %*% t(x)
  K = K + KPPI

  ## if(length(unknowns) > 0) # experimentell!!
  ##   K = K + min(R)*tcrossprod(x[,unknowns,drop=FALSE])
  best.bound = Inf		
  for(C in Cs){
    if(DEBUG) cat('Trying C=',C,'\n')
    K2 = as.kernelMatrix(K + 1/C*diag(NROW(x)))
    fit.tmp = ksvm(K2, y, C=Inf, type="C-svc", shrinking=FALSE)			
    cat("finished!\n")
    ##bound = spanbound(fit.tmp, K2, as.numeric(as.matrix(y)))
    bound = spanbound(fit.tmp, K2, sign(as.numeric(y) - 1.5))
    if(bound < best.bound){
      fit = fit.tmp
      best.bound = bound
      Cbest = C			
    }
  }
  if(DEBUG) cat('Best C=',Cbest,'\n')
  
  svs = unlist(alphaindex(fit))
  alpha = as.matrix(unlist(coef(fit)))
  if(calcW){
    if(DEBUG) cat("Calculating feature influences...")
    w = double(ncol(x))
    names(w) = colnames(x)
    for(i in colnames(x)){			
      s = setdiff(colnames(x), i)		
      dKi = 2*x[svs,s,drop=FALSE]%*%R[s,i,drop=FALSE]%*%t(x[svs,i,drop=FALSE]) + R[i,i]*tcrossprod(x[svs,i,drop=FALSE])			
      w[i] = t(alpha)%*%dKi%*%alpha
    }	
    #w[unknowns] = (t(alpha)%*%x[svs,unknowns,drop=FALSE])^2
    if(DEBUG) cat("finished!\n")		
  }
  else{
    w = NULL		
  }	
  list(fit.svm=fit, C=Cbest, R=R, w=w, error.bound=best.bound, xsvs=x[svs,,drop=FALSE], scale.mean=scale.mean, scale.std=scale.std, features=colnames(x))
}

#' Predict Method for Graph-SVM Fits
#'
#' Obtains predictions from a fitted graphSVM object.
#'
#' @param object a fitted object of class inheriting from 'graphSVM'
#' @param newdata a matrix with variables to predict
#' @param type \code{response} gives the predictions \code{class} gives the predicted classes.
#' @param ... currently ignored.
#' @return the predictions.
#' @export
#' @callGraphPrimitives
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(pathClass)
#' data(example_data)
#' matched <- matchMatrices(x=x, adjacency=adjacency.matrix, mapping=mapping)
#' dk <- calc.diffusionKernel(L=matched$adjacency, is.adjacency=TRUE, beta=0) # beta should be tuned
#' fit <- fit.graph.svm(matched$x[1:5,], y[1:5], DEBUG=TRUE, mapping=matched$mapping, diffusionKernel=dk)
#' predict(fit, newdata=matched$x[6:10,])
#' }
predict.graphSVM <- function(object, newdata, type="response", ...){

  fit <- object$fit
  
  ## do the prediction only with those genes
  ## that were use for training
  newdata <- newdata[,fit$features]

  if(!is.null(fit$scale.mean))
    newdata <- scale(newdata, center=fit$scale.mean[fit$features], scale=FALSE)
  
  if(!is.null(fit$scale.std))
    newdata <- scale(newdata, center=FALSE, scale=fit$scale.std[fit$features])

  ## Ktst = newdata[,fit$PPIgenes,drop=FALSE] %*% fit$R %*% t(fit$xsvs[,fit$PPIgenes,drop=FALSE])
  Ktst = newdata %*% fit$R %*% t(fit$xsvs)
  
  alpha = as.matrix(unlist(coef(fit$fit.svm)))	
  yhat = Ktst%*%alpha - b(fit$fit.svm)			

  if(type == "class" | is.null(fit$w))
    yhat = sign(yhat)			

  return(yhat)
}

#' Calculation of diffusion kernel matrix
#'
#' Calculation of diffusion kernel matrix
#'
#' @param L Laplacian or transition probability matrix
#' @param is.adjacency is L a laplace or adjacency matrix
#' @param beta beta parameter of the diffusion kernel. beta controls the extent of diffusion.
#' @return the diffusion kernel
#' @references Schoelkopf B, Tsuda K, Vert JP: Kernel Methods in Computational Biology, MIT Press; 2004.
#' @export
#' @callGraphPrimitives
#'
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
calc.diffusionKernel = function(L, is.adjacency=FALSE, beta=0){
  if(missing(L)) stop('You have to provide either a laplace or a adjacency matrix to calculate the diffusion kernel.')
  method="thresh"
  if(is.adjacency){
    dnames <- dimnames(L)
    L = graph.adjacency(L, mode='undirected', diag=FALSE)
    L = graph.laplacian(L, normalized=TRUE)
    dimnames(L) = dnames
  }
  
  if(method == "thresh"){
    eig = eigen(L)
    ## Let's discard 80% of the eigenvector with largest eigenvalues
    ncomp = round(0.8*ncol(L))+1
    V = eig$vectors[,ncol(L):ncomp]
    ## diff kernel
    R = V %*% diag(exp(-beta*eig$values[ncol(L):ncomp])) %*% t(V)
  }
## erstmal nur mehtod = 'thres zulassen'
##  else if(method == "pseudoinv")
##    R = pseudoinverse(L)
##  else if(method == "LLE"){	
##    if(!is.null(G)){
##      G = ugraph(G)
##      W = as(G, "matrix")
##      deg = rowSums(W)
##      Dinv = diag(1/deg)
##      Dinv[Dinv == Inf] = 1	
##      P = Dinv%*%W
##    }		
##    else{
##      P = L
##    }					
##    I = diag(ncol(P)) 
##    T = I - P
##    M = t(T)%*%(T)	
##    lam = eigen(M, symmetric=TRUE, only.values=TRUE)
##    K = lam$values[1]*I - M
##    E = I - matrix(1,ncol=ncol(I), nrow=nrow(I))/ncol(I)
##    R = E%*%K%*%E		
##  }
  dimnames(R) = dimnames(L)
  ## Normalize Kernel
  KernelDiag <- sqrt(diag(R) + 1e-10)
  R.norm <- R/(KernelDiag %*% t(KernelDiag))
  R.norm
}
