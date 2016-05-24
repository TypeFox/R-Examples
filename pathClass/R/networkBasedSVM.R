#' Implementation of the network-based Support Vector Machine introduced by Yanni Zhu et al., 2009.
#'
#' \code{mapping} must be a data.frame with at least two columns. The column names have to be \code{c('probesetID','graphID')}.
#' Where 'probesetID' is the probeset ID present in the expression matrix (i.e. \code{colnames(x)}) and 'graphID' is any ID that
#' represents the nodes in the diffusionKernel (i.e. \code{colnames(diffusionKernel)} or \code{rownames(diffusionKernel)}). The purpose of the this mapping is that
#' a gene or protein in the network might be represented by more than one probe set on the chip. Therefore, the algorithm must know
#' which genes/protein in the network belongs to which probeset on the chip.
#'
#' @param exps a p x n matrix of expression measurements with p samples and n genes.
#' @param y a factor of length p comprising the class labels.
#' @param DEBUG should debugging information be plotted.
#' @param n.inner number of fold for the inner cross-validation.
#' @param scale a character vector defining if the data should be centered and/or scaled.
#' Possible values are \emph{center} and/or \emph{scale}. Defaults to \code{c('center', 'scale')}.
#' @param sd.cutoff a cutoff on the standard deviation (sd) of genes. Only genes with sd > sd.cutoff stay in the analysis.
#' @param lambdas a set of values for lambda  regularization parameter of the L\eqn{_\infty}-Norm. Which, if properly chosen,
#'        eliminates factors that are completely irrelevant to the response, what in turn
#'        leads to a factor-wise (subnetwork-wise) feature selection. The 'best' lambda is
#'        found by an inner-cross validation.
#' @param adjacencyList a adjacency list representing the network structure. The list can be generated from a
#'        adjacency matrix by using the function \code{\link{as.adjacencyList}}
#' @return a networkBasedSVM object containing \item{features}{the selected features} \item{lambda.performance}{overview how different values of lambda performed in the inner cross validation} \item{fit}{the fitted network based SVM model}
#' @references Zhu Y. et al. (2009). Network-based support vector machine for classification of microarray samples. \emph{BMC Bioinformatics}
#' @export
#' @callGraphPrimitives
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
#' ad.list <- as.adjacencyList(matched$adjacency)
#' res.nBSVM <- crossval(matched$x, y, theta.fit=fit.networkBasedSVM, folds=3, repeats=1, DEBUG=TRUE, parallel=FALSE, adjacencyList=ad.list, lambdas=10^(-1:2), sd.cutoff=50)
#' }
fit.networkBasedSVM <- function (exps, y, DEBUG=FALSE, n.inner=3, scale=c('center', 'scale'), sd.cutoff=1, lambdas=10^(-2:4), adjacencyList){

  if(is.factor(y)) y <- sign(as.numeric(y) - 1.5)
  
  if(missing(adjacencyList)) stop('You need to provide an adjacency matrix.')

  ## discard low variance genes
  sd.exps <- apply(exps,2,sd)
  keep.sd <- colnames(exps)[which(sd.exps>=sd.cutoff)]
  if(DEBUG) cat("After fitering based on sd, ",length(keep.sd), " genes remained",sep="")

  keep.net  <- adjacencyList[(adjacencyList[,1] %in% keep.sd)&(adjacencyList[,2] %in% keep.sd),]
  keep.nnb  <- table(as.vector(keep.net))
  keep.nb   <- keep.net[,2]
  keep.gene <- names(keep.nnb)
  keep.exps <- exps[,keep.gene]
#  keep.name <- genename[keep.gene,]
#  keep.name$keepid<-c(1:dim(keep.exps)[2])
  keep.nnb  <-data.frame(keep.nnb)
  colnames(keep.nnb)<-c("id","nnb")
  duplicate <- which(keep.net[,1]>keep.net[,2])
  keep.net  <- keep.net[-duplicate,]

  ## things for inner-cv
  n.lambdas   <- length(lambdas)
  n.trn       <- length(y)
  inner.perm  <- sample(1:n.trn)
  lambda.perf <- matrix(0, ncol=n.lambdas, nrow=n.inner)
  dimnames(lambda.perf) <- list(paste('inner.fold', 1:n.inner, sep=''),lambdas)

  for(i in 1:n.inner){
    inner.tst = inner.perm[seq(i, n.trn, by=n.inner)]
    inner.trn = setdiff(1:n.trn, inner.tst)

    for(z in 1:n.lambdas){
      if(DEBUG) cat("\ninner CV fold ", i,". Testing Lambda ",lambdas[z],sep="")
      inner.fit  <- calc.networkBasedSVM(keep.exps[inner.trn,], y[inner.trn], lambda=lambdas[z], nnb=keep.nnb, unique.net=keep.net, keep.gene=keep.gene, scale=scale, DEBUG=DEBUG)
      inner.pred <- predict.networkBasedSVM(inner.fit, keep.exps[inner.tst,])
      l.perf    <- calc.auc(inner.pred, y[inner.tst])
      if(l.perf < 0.5) l.perf <- 1 - l.perf
      lambda.perf[i,z] <- l.perf
    }
  }

  ## select best lambda
  average.lambda.perf <- colMeans(lambda.perf)
  best.lambda         <- lambdas[which.max(average.lambda.perf)]
  best.lambda.auc     <- average.lambda.perf[which.max(average.lambda.perf)]
  if(DEBUG) cat("\nLambda choosen as: ",best.lambda,sep="")

  ## train final model
  complete.model <- calc.networkBasedSVM(keep.exps, y, lambda=best.lambda, nnb=keep.nnb, unique.net=keep.net, keep.gene=keep.gene, scale=scale, DEBUG=DEBUG)
  result = list(beta=complete.model$beta.hat, features=complete.model$features, lambda=complete.model$lambda, lambda.performance=lambda.perf, scale.mean=complete.model$scale.mean, scale.std=complete.model$scale.std)
  class(result) <- 'networkBasedSVM'

  return(result)
}


calc.networkBasedSVM <- function(x, tr.y, lambda, nnb, unique.net, keep.gene, scale, DEBUG=FALSE){

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

  ## anz samples
  n      <- nrow(x)
  n.gene <- ncol(x)
  p      <- sum(nnb[,2])/2
  
  ## Matrix of numeric constraint coefficients, one row per constraint, one column per variable
  f.con <- create.constraints.matrix(x, tr.y, nnb, unique.net, keep.gene)
  ## Vector of character strings giving the direction of the constraint
  f.dir <- c(rep(">=",dim(f.con)[1]))
  ## Vector of numeric values for the right-hand sides of the constraints
  b <- c(rep(1,nrow(x)),rep(0,2*p+2*n.gene+nrow(x)))
  ## Numeric vector of coefficients of objective function
  f.obj <- cbind(c(rep(0,2+2*n.gene),rep(1,nrow(x)),rep(lambda,p)))
  gc()
  
  ## beta = coefficient vector
  if(DEBUG) cat('\toptimizing using LpSolve...')
  beta.tr<-lp(direction="min",objective.in=f.obj,const.mat=f.con,const.dir=f.dir,const.rhs=b)
  if(DEBUG) cat('Done.')
  beta.pos<-c(beta.tr$solution[1],beta.tr$solution[3:(2+n.gene)])
  beta.neg<-c(beta.tr$solution[2],beta.tr$solution[(2+n.gene+1):(2+2*n.gene)])
  beta.est<-as.matrix(beta.pos-beta.neg)

  list(beta.hat=beta.est, lambda=lambda, scale.mean=scale.mean, scale.std=scale.std, features=colnames(x))
}

#' Predict Method for Network-based SVM Fits
#'
#' Obtains predictions from a fitted networkBasedSVM object.
#'
#' @param object a fitted object of class inheriting from 'networkBasedSVM'
#' @param newdata a matrix with variables to predict
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
#' ad.list <- as.adjacencyList(matched$adjacency)
#' fit = fit.networkBasedSVM(matched$x[1:5,], y[1:5], DEBUG=TRUE,  adjacencyList=ad.list, lambdas=10^(-1:2), sd.cutoff=50)
#' predict(fit, newdata=matched$x[6:10,])
#' }
predict.networkBasedSVM <- function(object, newdata, ...){

  ## do the prediction only with those genes
  ## that were use for training
  newdata <- newdata[,object$features]

  if(!is.null(object$scale.mean))
    newdata <- scale(newdata, center=object$scale.mean[object$features], scale=FALSE)
  
  if(!is.null(object$scale.std))
    newdata <- scale(newdata, center=FALSE, scale=object$scale.std[object$features])

  cbind(rep(1,nrow(newdata)), newdata) %*% object$beta
}

create.constraints.matrix <- function(tr.exps, tr.y, nnb, unique.net, keep.gene){
  ##########################
  ## generate Aalt.matrix ##
  ##########################

  ## anz samples
  n<-dim(tr.exps)[1]
  n.gene<-dim(tr.exps)[2]
  p<-sum(nnb[,2])/2

  ## anz samples + 2*Summe Verbindungen/2 + 2*Anz gene + anz samples = 276770
  a.nrow<-n+2*p+2*n.gene+n

  ## 2 + 2*Anz gene + Anz samples + Summe Verbindungen/2 = 151519
  a.ncol<-2+2*n.gene+n+p

  a<-matrix(0,nrow=a.nrow,ncol=a.ncol)

  ## erste spalte = y
  a[1:n,1]<-tr.y;

  ## zweite spalte = -y
  a[1:n,2]<--tr.y

  ## expressions werte mit y und -y multiplizieren
  a[1:n,3:(n.gene+2)]<-tr.exps*tr.y; a[1:n,(n.gene+3):(2*n.gene+2)]<-tr.exps*(-tr.y)

  ## einheits matrix samples x samples
  a[1:n,(2*n.gene+3):(a.ncol-p)]<-diag(1,n,n)

  ## einheits matrix (2 anz gene + anz samples) x (2 anz gene + anz samples)
  a[(n+2*p+1):a.nrow,3:(a.ncol-p)]<-diag(1,2*n.gene+n,2*n.gene+n)

  temp.size<-size<-0
  K<-length(table(unique.net[,1]))
  for(k in 1:K){
    temp.size<-table(unique.net[,1])[k];
    size<-size+temp.size
    subnet<-unique.net[(size-temp.size+1):size,2]
    ## node<-as.numeric(names(temp.size)); ## Does not work if we have real node-names like NP_1234
    node<-names(temp.size);
    loc.node<-which(keep.gene==node)
    for(l in 1:temp.size){
      loc.end<-which(keep.gene==subnet[l])
      a[n+2*(size-temp.size+l)-1,2+loc.node]<-a[n+2*(size-temp.size+l)-1,2+n.gene+loc.node]<--1/sqrt(nnb$nnb[nnb$id==node])
      a[n+2*(size-temp.size+l),2+loc.end]<-a[n+2*(size-temp.size+l),2+n.gene+loc.end]<--1/sqrt(nnb$nnb[nnb$id==subnet[l]])
    }
  }

  a[(n+1):(n+2*p),(2+2*n.gene+n+1):a.ncol]<-bdiag(p)

  a
}

## create block diagonal matrix with block (1,1)
## n=number of blocks
bdiag<-function(n){
  a<-matrix(0,nrow=2*n,ncol=n)
  for(i in 1:n){
    a[(2*i-1),i]<-1
    a[(2*i),i]<-1
  }
  return(a) 
}

#' Uses a adjacency matrix to create a adjacency list
#'
#' Uses a adjacency matrix to create a adjacency list as needed for \code{\link{fit.networkBasedSVM}}.
#'
#' @param adjacency.matrix a adjacency matrix.
#' @param skip.redundant.nodes if \code{TRUE} and the graph is undirected only the upper triangular matrix (including the diagonal) is used to create the adjacency list.
#' @param is.directed determines wether or not the graph is directed.
#' @return an adjacency list.
#' @export
#' @callGraphPrimitives
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(pathClass)
#' data(adjacency.matrix)
#' ad.list <- as.adjacencyList(adjacency.matrix)
#' }
as.adjacencyList <- function (adjacency.matrix, skip.redundant.nodes = TRUE, is.directed = FALSE){
  if (ncol(adjacency.matrix) != nrow(adjacency.matrix)) {
    stop("Matrix is non square")
  }
  if (is.null(rownames(adjacency.matrix)) & is.null(colnames(adjacency.matrix))) 
    rownames(adjacency.matrix) <- colnames(adjacency.matrix) <- 1:ncol(adjacency.matrix)
  else if (is.null(rownames(adjacency.matrix)) & !is.null(colnames(adjacency.matrix))) 
    rownames(adjacency.matrix) <- colnames(adjacency.matrix)
  else if (!is.null(rownames(adjacency.matrix)) & is.null(colnames(adjacency.matrix))) 
    colnames(adjacency.matrix) <- rownames(adjacency.matrix)
  
  adjacency.list = cbind(rep(rownames(adjacency.matrix), times = ncol(adjacency.matrix)), rep(rownames(adjacency.matrix), each = nrow(adjacency.matrix)))

  if (skip.redundant.nodes & !is.directed) {
    print("Skipping lower triangular")
    adjacency.matrix[lower.tri(adjacency.matrix, diag = F)] = 0
  }
  
  adjacency.list = adjacency.list[as.logical(adjacency.matrix == 1), ]
  adjacency.list
}
