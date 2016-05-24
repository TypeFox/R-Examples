# Various utility functions
#
# Date: 
#   Revised: February 19, 2015
#
# Author: Jiarui Ding <jiaruid@cs.ubc.ca>
#   Department of Computer Science, UBC
#   Department of Molecular Oncology, BC Cancer Agency 
#

# Various functions for data pre-processing
#
# Date: 
#   Revised: February 18, 2015
#
# Author: Jiarui Ding <jiaruid@cs.ubc.ca>
#   Department of Computer Science, UBC
#   Department of Molecular Oncology, BC Cancer Agency 
#


## Remove the colnames with too many NAs
RemoveLowCoverage = function(X, ratio=0.5) {
  count = colSums(is.na(X))
  id = which(count >= ratio * nrow(X))
  
  return(id)
}


#==============================================================================
#' A mixture modelling approach to estiamte whether a gene is expressed in a study
#'  given RNA-seq gene expression data
#' 
#' @export
#' @param expr A matrix of RNA-seq gene expression values where 
#' each row corresponds to a patient and each column is a gene.
#' Typically the expression of each gene is the log2 transformed RSEM value. 
#' @param show.plot Logical, specifying whether to plot results
#' @param loglik Logical, whether plot the log-likelihoods
#' @param xlab xlab of the plot
#' @param ylab ylab of the plot
#' @param ... Arguments for plotting
#' @return A weight vector representing whether individual genes 
#' are expressed in the study
#' @importFrom stats setNames quantile sd 
#' @examples
#' data(expr)
#' weight = EstimateExpression(expr)
EstimateExpression = function(expr, show.plot=FALSE, loglik=TRUE, 
                      xlab="Expression", ylab="Density", ...) {
  expr.quantile = apply(expr, 2, quantile, 0.90, na.rm=TRUE)
  id = is.na(expr.quantile)
  gene.no.expr = names(which(id))
  
  expr.quantile = expr.quantile[!id]
  out = DetectBoxplotOutlier(expr.quantile)
  expr.quantile.out = expr.quantile[out$outlier.id]
  expr.quantile = expr.quantile[out$keep.id]
  
  th = quantile(expr.quantile, 0.30)
  mu = c(mean(expr.quantile[expr.quantile<=th]), 
         mean(expr.quantile[expr.quantile>th]))
  sigma  = c(sd(expr.quantile[expr.quantile<=th]), 
             sd(expr.quantile[expr.quantile>th]))
  lambda = c(0.30, 0.70)
  
  model = MixGaussFitEM(expr.quantile, lambda=lambda, 
                        mu=mu, sigma=sigma, K=2)
  if(show.plot == TRUE) {
    MixModelPlot(model,  xlab2=xlab, ylab2=ylab, 
                 breaks=40, loglik=loglik, ...)
  }
  
  weight = model$posterior[, 2]
  names(weight) = names(expr.quantile)
  
  weight.no.expr = setNames(rep(0, length(gene.no.expr)), nm=gene.no.expr)
  id = names(which(expr.quantile.out > model$mu[2]))
  weight.out.up = setNames(rep(1, length(id)), id)
  
  id = names(which(expr.quantile.out < model$mu[1]))
  weight.out.down = setNames(rep(0, length(id)), id)
  
  weight.out = c(weight.out.down, weight.out.up)
  weight = c(weight, weight.no.expr, weight.out)
  
  return(weight)
}

#==============================================================================
#' Impute missing values (NAs) using K-nearest neighbour averaging
#' 
#' @export
#' @param X A matrix of real values where each row corresponds to a patient 
#' and each column is a gene.
#' @param ratio The rows (columns) with more than (default 50\%) of 
#' missing values are removed
#' @param ... Arguments to be passed to impute.knn
#' @return A matrix without NAs
#' @examples
#' data(expr)
#' expr.norm = ImputeKnn(expr)
ImputeKnn = function(X, ratio=0.5, ...) {  
  id.col = RemoveLowCoverage(X, ratio)
  id.row = RemoveLowCoverage(t(X), ratio)
  
  if (length(id.col) > 0)
    X = X[, -id.col]
  if (length(id.row) > 0)
    X = X[-id.row, ]
  
  X = t(X)
  X = impute::impute.knn(data=X, ...)
  X = t(X$data)
  
  return(X)
}


#==============================================================================
#' Quantile normalize a matrix
#' 
#' @export
#' @param X A matrix of real values where each row corresponds to a patient 
#' and each column is a gene
#' @return The normalized matrix of \code{X}
#' @examples
#' data(expr)
#' expr.quantile = QuantileNorm(expr)
QuantileNorm = function(X) {
  patient = rownames(X)
  gene = colnames(X)
  
  X = preprocessCore::normalize.quantiles(t(X))
  X = t(X)
  
  rownames(X) = patient
  colnames(X) = gene
  
  return(X)
}


#==============================================================================
#' @importFrom graphics par lines
#' @importFrom stats predict
NormExprSub = function(cna.logr, expr, sigma=0.2, 
                       type="svm", show.plot=FALSE, show.norm=TRUE) {
  id = order(cna.logr)
  expr = expr[id]
  cna.logr = cna.logr[id]
  
  x = cna.logr
  if(type == "gp") {
    options = gptk::gpOptions()
    options$kern$comp = list("rbf", "white")
    model = gptk::gpCreate(1, 1, as.matrix(cna.logr), as.matrix(expr), options)
    
    y = gptk::gpPosteriorMeanVar(model, as.matrix(cna.logr), varsigma.return=TRUE)
    expr.pred = y$mu
  } else {
    model = e1071::svm(x, expr)
    expr.pred = predict(model, x)
  }
  expr.corr = expr - expr.pred
  
  if (show.plot) {
    if (show.norm == TRUE) {
      op = par(mfrow=c(1,2))
    }
    NeatPlot(cna.logr, expr, xlab="CNA", ylab="Expression", 
             col=AddTrans("dodgerblue", 0.2))
    lines(cna.logr, expr.pred)
    
    if(type == "gp") {
      lines(cna.logr, y$mu - sqrt(y$varsigma)*2, col="red")
      lines(cna.logr, y$mu + sqrt(y$varsigma)*2, col="red")
    }    
    if (show.norm == TRUE) {
      NeatPlot(cna.logr[order(id)], expr.corr[order(id)], xlab="CNA", 
               ylab="Expression", col=AddTrans("dodgerblue", 0.2))
      par(op)
    } 
  }
  
  return(expr.corr[order(id)])
}


#==============================================================================
#' Remove the cis-effects of copy number alterations on gene expression
#' 
#' @export
#' @param cna.logr A matrix of copy number alterations log2 ratio 
#' where each row corresponds to a patient and each column is a gene.
#' @param expr A matrix of gene expression values where 
#' each row corresponds to a patient and each column is a gene.
#' @param gene A character vector of gene HGNC symbols, 
#' default for all genes with both gene expression and 
#' copy number log2 ratio data
#' @param type A character, either Gaussian process regression  ("gp") or 
#' suppor vector machine regression ("svm")
#' @param debug Logical, specifying whether debug information should be printed
#' @param show.plot Logical, specifying whether to plot the orginal expression 
#' and the normalized expression for a gene
#' @param show.norm Logical, specifying whether to plot the express of 
#' a gene after normalization, only used when \code{show.plot = TRUE}.
#' @return The normalized expression matrix
#' 
#' @importFrom graphics hist
#' @importFrom stats median
#' 
#' @examples
#' data(cna.logr, expr)
#' expr.norm = NormExpr(cna.logr, expr, gene="PTEN")
NormExpr = function(expr, cna.logr, gene, type="gp", 
                    debug=FALSE, show.plot=FALSE, show.norm=TRUE) {
  #
  gene.common = intersect(colnames(expr), colnames(cna.logr))
  if(!missing(gene)) {
    gene = intersect(gene, gene.common)
  } else {
    gene = gene.common
  }
  
  sample = intersect(rownames(expr), rownames(cna.logr))
  expr.corr = expr[sample, gene, drop=FALSE]
  for (gene.i in gene) {
    if(debug) {
      cat(gene.i, "\n")
    }
    
    sigma = median(abs(cna.logr[sample, gene.i])) / 1
    if(sigma != 0) {
      y.corr = NormExprSub(cna.logr[sample, gene.i], expr[sample, gene.i], 
                           show.plot=show.plot, show.norm=show.norm, 
                           sigma=sigma, type=type)
    } else {
      y.corr = expr[sample, gene.i] - median(expr[sample, gene.i])
      if(show.plot == TRUE) {
        hist(y.corr, main="", xlab="Exprssion")
      }
    }
    expr.corr[, gene.i] = y.corr
  }
  return(expr.corr)
}

#' Filter network 
#' 
#' @export
#' @param net List, a gene interaction network
#' @param weight The weights of genes, could from the function \code{EstimateExpression}
#' @param min.weight Filter the connected genes with weights less than \code{min.weight}
#' @param min.conn.strength The minimum gene connection strength
#' @param remove.self.connection Logical, whether removing self-connections or not
#' @param max.num.conn Only keep the top max.conn genes
#' @param min.num.conn The minimum number of connections required for a gene to be considered for trans-analysis
#' @param debug Logical, specifying whether debug information should be printed
#' 
#' @return The filtered network
#' @examples
#' data(net)
#' net.filt = FilterNetwork(net)
#' 
#==============================================================================
FilterNetwork = function(net, weight, min.weight=0.8, min.conn.strength=0.4,  
                         min.num.conn=5, max.num.conn=50, remove.self.connection=TRUE, 
                         debug=FALSE) { 
  
  weight.status = !missing(weight)
  
  net = sapply(net, function(z) {
    z = sort(z[z >= min.conn.strength], decreasing=TRUE)
    
    if (weight.status == TRUE) {
      gene = names(z)
      weight.gene = weight[gene]
      
      id = weight.gene >= min.weight
      id = id & !is.na(id)
      z = z[id]
    }
    
    if (length(z) > max.num.conn) {
      z = z[1:max.num.conn]
    }
    if (length(z) < min.num.conn) {
      z = NULL
    }
    return (z)
  })
  
  net = Filter(length, net)

  return(net)
}

#==============================================================================
ConvertMutType = function(mut) {
  id = mut[, "variant_type"] %in% c("HOMD", "HLAMP")
  mut.somatic = mut[!id, , drop=FALSE]
  mut.cna     = mut[id, , drop=FALSE]
  
  sample.gene = paste(mut.somatic[, "sample"], 
                      mut.somatic[, "hgnc_symbol"], 
                      sep="_")
  id = duplicated(sample.gene) | duplicated(sample.gene, fromLast=TRUE)
  
  if(sum(id) > 0) {
    mut.somatic[id, "variant_type"] = "COMPLEX" 
  }
  mut = rbind(mut.somatic, mut.cna)
  
  return(mut)
}


#==============================================================================
Trim = function (str) {
  # returns a string w/o leading or trailing whitespaces
  
  gsub("^\\s+|\\s+$", "",  str)
}


#==============================================================================
## From gtools package, to prevent confict with e1071
invalid = function(x) {
  if (missing(x) || is.null(x) || length(x)==0) {
    return(TRUE)
  }
  
  if (is.list(x)) {
    return (all(sapply(x, invalid)))
  } else if (is.vector(x)) {
    return (all(is.na(x)))
  } else {
    return (FALSE)
  }
}


