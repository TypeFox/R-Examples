#' Reweighted Recursive Feature Elimination (RRFE)
#'
#' Implementation of the Reweighted Recursive Feature Elimination (RRFE) algorithm.
#' \code{mapping} must be a data.frame with at least two columns. The column names have to be \code{c('probesetID','graphID')}.
#' Where 'probesetID' is the probeset ID present in the expression matrix (i.e. \code{colnames(x)}) and 'graphID' is any ID that
#' represents the nodes in the graph (i.e. \code{colnames(Gsub)} or \code{rownames(Gsub)}). The purpose of the this mapping is that
#' a gene or protein in the network might be represented by more than one probe set on the chip. Therefore, the algorithm must know
#' which genes/protein in the network belongs to which probeset on the chip. However, the method is able to use all feature when one
#' sets the parameter \code{useAllFeatures} to \code{TRUE}. When doing so, RRFE assigns the minimal wheight returned by GeneRank
#' to those genes which are not present in \code{Gsub}.
#'
#' @param x a p x n matrix of expression measurements with p samples and n genes.
#' @param y a factor of length p comprising the class labels.
#' @param DEBUG should debugging information be plotted.
#' @param scale a character vector defining if the data should be centered and/or scaled.
#' Possible values are \emph{center} and/or \emph{scale}. Defaults to \code{c('center', 'scale')}.
#' @param Cs soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
#' @param stepsize amount of features that are discarded in each step of the feature elimination. Defaults to 10\%.
#' @param useAllFeatures should all features be used for classification (see also \code{Details}).
#' @param mapping a mapping that defines how probe sets are summarized to genes.
#' @param Gsub an adjacency matrix that represents the underlying biological network.
#' @param d the damping factor which controls the influence of the network data and the fold change on the ranking of the genes.
#' Defaults to 0.5 
#' @return a RRFE fit object.
#' \item{features}{the selected features}
#' \item{error.bound}{the span bound of the model}
#' \item{fit}{the fitted SVM model}
#' @references Johannes M, et al. (2010). Integration Of Pathway Knowledge Into A Reweighted Recursive Feature Elimination Approach For Risk Stratification Of Cancer Patients. \emph{Bioinformatics}
#' @export
#' @callGraphPrimitives
#' @note The optimal number of features is found by using the span estimate. See Chapelle, O., Vapnik, V., Bousquet, O., and Mukherjee, S. (2002). Choosing multiple parameters for support vector machines. \emph{Machine Learning}, 46(1), 131-159.
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
#' res.rrfe <- crossval(x, y, DEBUG=TRUE, theta.fit=fit.rrfe, folds=3, repeats=1, parallel=TRUE, Cs=10^(-3:3), mapping=mapping, Gsub=adjacency.matrix, d=1/2)
#' }
fit.rrfe = function(x, y, DEBUG=FALSE, scale=c('center', 'scale'), Cs=10^c(-3:3), stepsize=0.1, useAllFeatures=F,  mapping, Gsub, d=0.5){

  best.bound = Inf
  feat = colnames(x)

  if(missing(mapping)) stop('You must provide a mapping.')
  if(missing(Gsub))    stop('RRFE needs the underlying grap structure.')

  ## discard rows from the mapping which are not needed
  ## because they are not present in x
  ## which might happen if the used does not use the full chip
  mapping <- mapping[mapping[,'probesetID'] %in% colnames(x),]
  if(nrow(mapping) == 0) stop('Probeset IDs of the mapping do not fit to the column names of x!')

  int <- intersect(rownames(Gsub), mapping[,"graphID"])
  Gsub <- Gsub[int, int]
  mapping <- mapping[mapping[,'graphID'] %in% int,]
  ## overlap <- rownames(Gsub) %in% mapping[,"graphID"]
  ## Gsub    <- Gsub[overlap, overlap]
  ## if(nrow(Gsub) == 0 || ncol(Gsub) == 0) stop('The row/column names of Gsub do not fit to the graph ID of the mapping!')

  ## 3. den chip nun dem mapping anpassen
  x <- x[,mapping[,'probesetID']]
  
  if(DEBUG) cat('Calculating GeneRanks...')
  ## calculate the fold change for a gene by using all probes
  ## targeting that gene
  exprs.sum <- summarizeProbes(exprs=t(x), mapping=mapping, method="foldChange", groups=y, adjacency=Gsub)
  ranks     <- getRanking(Gsub, exprs.sum, d=d)
  ## distribute the rank of the gene back to its probes
  ranks     <- desummarize.ranks(ranks, mapping)
  if(DEBUG) cat('done\n')

  if(useAllFeatures){
    ## create vector with minimum google weigth
    ranks.complete        <- rep(min(ranks),NCOL(x))
    names(ranks.complete) <- colnames(x)
    ## save the 'real' weigth for those genes pathway
    ## knowledge was available
    ranks.complete[names(ranks)] <- ranks
    ranks <- ranks.complete
  }
  ## use only probes we have ranks for
  x <- x[,names(ranks)]

  while(NCOL(x) > 1){		
    if(DEBUG) cat(NCOL(x),' Features left.\n')

    fit = svm.fit(x=x, y=y, Cs=Cs, scale=scale, DEBUG=DEBUG)
    
    if(fit$error.bound <= best.bound){
      best.bound = fit$error.bound
      feat  = colnames(x)			
      best.fit   = fit
      if(DEBUG) cat('Model Updated. Spanbound=',best.bound,', C=',best.fit$C,', ', length(feat),'features.\n')
    }

    fit$w = as.vector(fit$w) * ranks[colnames(fit$w)]
    
    ord      = order(fit$w)
    remove   = colnames(x)[ord[1:round(NCOL(x)*stepsize)]]
    x        = x[, setdiff(colnames(x), remove), drop=FALSE]
  }
  if(DEBUG) cat('Best Model is: Spanbound=',best.fit$error.bound,', C=',best.fit$C,',', length(best.fit$features),'features.')
  
  result <- list(features=feat, error.bound = best.bound, fit=best.fit)
  class(result) <- 'rrfe'
  return(result)
}

#' Predict Method for RRFE Fits
#'
#' Obtains predictions from a fitted RRFE object.
#'
#' @param object a fitted object of class inheriting from 'rrfe'
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
#' fit = fit.rrfe(x[1:5,], y[1:5], DEBUG=T, mapping=mapping, Gsub=adjacency.matrix)
#' predict(fit, newdata=x[6:10,])
#' }
predict.rrfe = function(object, newdata, type='response', ...){
    svm.predict(object$fit, newdata, type)
}

## Transforms the GeneRank result as
## as needed for RRFE
getRanking <- function(graph, exprs, d=0.5){

  adjac.matrix <- NULL
  
  if(class(graph) == "graphNEL"){
    igraph <- igraph.from.graphNEL(graph)
    adjac.matrix <- get.adjacency(igraph)
  }
  else if(class(graph) == "dgCMatrix" || class(graph) == "matrix")
    adjac.matrix = graph

  rm(graph)

  if(!all(colnames(adjac.matrix) %in% names(exprs)))
    stop("Couldn't put expr Values into right order, because the names didn't fit...\n")

  ## put the exprs values into right order
  ## as in the adjac matrix
  exprs <- exprs[colnames(adjac.matrix)]

  ranks        <- geneRank(ex=exprs,W=adjac.matrix,d=d)
  rank.names   <- colnames(adjac.matrix)

  ## take the reciprocal of the rank, so that the biggest
  ## google-rank is ranked on place one after applying rank().
  ## And then take the reciprocal again to tranlate it into
  ## rank-based google weights
  ranks        <- 1/rank(1/ranks)
  names(ranks) <- rank.names

  ranks  
}


#' Calculate GeneRanks as used by RRFE
#'
#' Uses the GeneRank  to calculate the ranks for genes. Afterwards the ranks 
#' are transformed as needed for the RRFE algorithm.
#'
#' @param x a p x n matrix of expression measurements with p samples and n genes.
#' @param y a factor of length p comprising the class labels.
#' @param mapping a mapping that defines how probe sets are summarized to genes.
#' @param Gsub an adjacency matrix that represents the underlying biological network.
#' @param method see help of \link{summarizeProbes}
#' @param d the damping factor which controls the influence of the network data and the fold change on the ranking of the genes.
#' Defaults to 0.5 
#' @return a ranking of the genes for which pathway knowledge was available.
#' @references Johannes M, et al. (2010). Integration Of Pathway Knowledge Into A Reweighted Recursive Feature Elimination Approach For Risk Stratification Of Cancer Patients. \emph{Bioinformatics}
#' @export
#' @callGraphPrimitives
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(pathClass)
#' data(example_data)
#' ranks = getGeneRanks(x, y, mapping=mapping, Gsub=adjacency.matrix)
#' }
getGeneRanks <- function(x, y, mapping, Gsub, method="foldChange", d=0.5){

  if(missing(mapping)) stop('You must provide a mapping.')
  if(missing(Gsub))    stop('RRFE needs the underlying grap structure.')

  ## discard rows from the mapping which are not needed
  ## because they are not present in x
  ## which might happen if the user does not use the full chip
  mapping <- mapping[mapping[,'probesetID'] %in% colnames(x),]
  if(nrow(mapping) == 0) stop('Probeset IDs of the mapping do not fit to the column names of x!')

  int <- intersect(rownames(Gsub), mapping[,"graphID"])
  Gsub <- Gsub[int, int]
  mapping <- mapping[mapping[,'graphID'] %in% int,]

  ## 3. fit chip to mapping
  x <- x[,mapping[,'probesetID']]
  
  ## calculate the fold change for a gene by using all probes
  ## targeting that gene
  exprs.sum <- summarizeProbes(exprs=t(x), mapping=mapping, method=method, groups=y, adjacency=Gsub)
  ranks     <- getRanking(Gsub, exprs.sum, d=d)
  ## distribute the rank of the gene back to its probes
  ranks     <- desummarize.ranks(ranks, mapping)

  return(ranks)
}
