#' Upper triangular loadings restriction
#'
#' Returns the restriction matrix corresponding to Geweke & Zhou (1996)
#' upper-triangular positive diagonal condition
#' '
#' @param p number of rows (variables)
#' @param k number of columns (factors)
#' @return p by k matrix
#' @export

utri_restrict = function(p, k) {
  out = matrix(1, nrow=p, ncol=k)
  out[upper.tri(out)] = 0
  diag(out) = 2
  return(out)
}

#' Initialize a bfa model
#'
#' This function accepts a data matrix \code{D} and specified options,
#' returning an S3 object of class bfa.
#'
#' @param x A formula or matrix
#' @param data The data if x is a formula
#' @param num.factor number of factors
#' @param restrict A matrix or list giving restrictions on factor loadings. A matrix should be the 
#' same size as the loadings matrix. Acceptable values are 0 (identically 0), 1 (unrestricted), 
#' or 2 (strictly positive). List elements should be character vectors of the form c('variable',1, '>0')
#' where 'variable' is the manifest variable, 1 is the factor, and '>0' is the restriction. Acceptable
#' restrictions are '>0' or '0'.
#' @param normal.dist A character vector specifying which variables should be treated as observed Gaussian
#' @param center.data Center each margin
#' @param scale.data Scale each margin to unit mean/variance
#' @param init Initialize the factor loadings 
#' @param ... ignored
#' @return An S3 object of class \code{bfa}. 

bfa_model <- function(x, data=NULL, num.factor=1, restrict=NA, 
                      normal.dist=NA, center.data=TRUE, scale.data=FALSE, 
                      init=TRUE, ...) {
  more_args = list(...)
  if (is.matrix(x)) D = x
  if(inherits(x, "formula")) {
    fr = model.frame(x, data=data, na.action=NULL)
    d = dim(fr)
    if(any(is.na(normal.dist))) {
      normal.dist = rep(0, d[2])
      fillnd = TRUE
    } else {
      fillnd = FALSE
    }
    D = matrix(NA, nrow=d[2], ncol=d[1])
    for (i in 1:d[2]) {
      if(class(fr[,i])[1] %in% c('ordered', 'numeric', 'integer')) {
        D[i,] = as.numeric(fr[,i])
        if('numeric' %in% class(fr[,i]) && fillnd) {
          normal.dist[i] = 1
        } else {
          if(fillnd) normal.dist[i] = 0
        }
      } else {
        stop("Data contain variables which are not ordinal. Each variable must be an integer/numeric vector or ordered factor")
      }
    }
    rownames(D) = colnames(fr); colnames(D) = rownames(fr)
  }
  
  if(is.numeric(normal.dist))   normal.var = matrix(normal.dist)
  if(any(is.na(normal.dist)))   normal.var = matrix(rep(0.0, nrow(D)))
  if(is.character(normal.dist)) normal.var = matrix(as.numeric(rownames(D) %in% normal.dist))
  
  D = as.matrix(D)
	D[is.infinite(D)] = NA
	D = t(scale(t(D), center=center.data, scale=scale.data))
  
  mInd = which(is.na(D[as.logical(normal.var),]), arr.ind=TRUE)
  
	p = dim(D)[1]
	n = dim(D)[2]
  k = num.factor

  U = D
  argsorts = matrix(0, nrow=p, ncol=n)
  Ra =  matrix(0, nrow=p, ncol=n)
  for (i in 1:p) {
    levs = sort(unique(D[i,]))
    asrow = which(is.na(D[i,]))
    rarow = rep(0, length(asrow))
    
    for (j in 1:length(levs)) {
      lev = levs[j]
      inds = which(D[i,]==lev)
      asrow = c(asrow, inds)
      Ra[i,inds] = rep(j, length(inds))
    }
    argsorts[i,] = asrow
    
    df = ecdf(D[i,])
    nobs = n - sum(is.na(D[i,]))
    U[i,] = df(D[i,])*nobs/(nobs+1)
  }
    
  argsorts = argsorts - 1 
  Ra = Ra - 1
  
	Z = qnorm(U)*sqrt(k)
  for(i in 1:p) {
    if(normal.var[i]>0) Z[i,] = D[i,]
  }

	maxes = matrix(rep(0, p*n), nrow=p)
	L = 0
	for (i in 1:p) {
		m = sort(unique(Z[i,]))#[-length(m)], 99999.0)
    m[length(m)] = 99999.0
		maxes[i,1:length(m)] = m
		L = max(L, length(m))
	}
	
	maxes = maxes[,1:L]

  Z[is.na(Z)] = rnorm(length(Z[is.na(Z)]))
  
  scores = matrix(rnorm(k*n), nrow=k)
  loadings = matrix(rnorm(k*p), nrow=p)
  
  if (any(is.na(restrict))) {
    loadings.restrict = matrix(1, nrow=p, ncol=k)
  } else if(any(restrict == "upper.tri")) {
    loadings.restrict = utri_restrict(p,k)
  } else if (is.list(restrict)) {
    loadings.restrict = matrix(1, nrow=p, ncol=k)
    for (r in restrict) {
      rowi = try(which(rownames(D)==r[1]))
      if(class(rowi) == "try-error") stop("Variable ",r[1]," not found in data")
      coli = as.numeric(r[2])
      #print(rowi); print(coli)
      if(r[3]=='>0') {
        loadings.restrict[rowi, coli] = 2
      } else if(as.numeric(r[3])==0) {
        loadings.restrict[rowi, coli] = 0
      }
    }
  } else if (is.matrix(restrict)) {
    loadings.restrict = restrict
  }
  
  
  if(init) {   
    initialized=TRUE
    sm = apply(Z, 1, mean)/apply(Z, 1, sd)
    sm = ifelse(abs(sm>0.5), sign(sm)*0.5, sm)
    for (i in 1:p) {
      for (j in 1:k) {
        if (loadings.restrict[i,j]==1) {
          loadings[i,j] = sm[i]
        } else {
          loadings[i,j] = loadings.restrict[i,j]
        }
      }
    }
    .updateScores(Z, loadings, scores)
  }
  
  if(any(is.null(more_args$init.fa))) more_args$init.fa=FALSE
  if(more_args$init.fa) {
    mf = factanal(t(Z), factors=num.factor, scores='regression', lower=0.05)
    loadings = (1/mf$uniqueness)*matrix(mf$loadings, ncol=num.factor)
    scores = t(mf$scores)
  }
  
  if(is.null(colnames(D))) colnames(D) = paste('Obs',1:ncol(D))
  if(is.null(rownames(D))) rownames(D) = paste('Var',1:nrow(D))
  
	out = list(data=D, ldata=Z, ranks=Ra, maxes=maxes, 	
				argsorts=argsorts, P = dim(D)[1], N = dim(D)[2], K=num.factor,
				obslabel=colnames(D), varlabel=rownames(D),
				nsim=0, nburn=0, thin=1, 
				
				original.data=data,
        mInd = mInd,
        
        loadings = loadings,
        post.loadings.mean = matrix(rep(0,k*p), ncol=k), 
        post.loadings.var  = matrix(rep(0,k*p), ncol=k),
        
        scores = scores,
        post.scores.mean = matrix(rep(0,k*n), nrow=k),
        post.scores.var  = matrix(rep(0,k*n), nrow=k),
        
        tauinv = rgamma(p, 1), 
        rho = rbeta(k, 10, 10),
        
        error_var_i = normal.var,
        loadings.restrict = loadings.restrict
        )
  
	class(out) = "bfa"
	attr(out, 'init')=initialized
  return(out)
}

#' Print method for bfa object
#' @param x A bfa object
#' @param ... Ignored
#' @method print bfa
#' @export
print.bfa <- function(x, ...) {
	cat("bfa model object\n")
	}

#' Extract posterior means from bfa object
#' @param x A bfa object
#' @param ... Ignored
#' @return A list with elements loadings and scores containing MCMC sample means
#' @method mean bfa
#' @export
mean.bfa <- function(x, ...) {
  return(list(loadings=x$post.loadings.mean, scores=x$post.scores.mean))
}

#' Extract posterior means from bfa object
#' @param x A bfa object
#' @param ... Ignored
#' @return A list with elements loadings and scores containing MCMC sample variances
#' @export
var.bfa <- function(x, ...) {
  return(list(loadings=x$post.loadings.var, scores=x$post.scores.var))
}

#' Extract samples of implied regression coefficients from a bfa object
#' @param object A bfa object
#' @param responses A character vector containing one or more response variables
#' @param scale Whether to compute regression coefficients from factor loadings on
#' the correlation scale; recommended if object is a copula or mixed factor model.
#' @param ... Ignored
#' @return An array of dimension length(index) x p-length(index) x (no. of mcmc 
#' samples) with posterior samples of regression coefficients
#' @method coef bfa
#' @export
coef.bfa <- function(object, responses, scale = attr(object, "type")!="gauss", ...) {
  
  index = match(responses, colnames(object$original.data))
  pl = object$post.loadings
  ns = dim(pl)[3]
  out = array(NA, dim = c(length(index), object$P - length(index), 
                          ns))
  dimnames(out) = list(object$varlabel[index], object$varlabel[-index], 
                       NULL)
  for (i in 1:ns) {
    if (attr(object, "type") == "gauss") {
      u = object$post.sigma2[i, ]
    }
    else {
      u = 1/(1 + rowSums(pl[, , i]^2))
      pl[, , i] = pl[, , i] * sqrt(u)
    }
    mat = t(pl[-index, , i]) %*% woodbury(pl[-index, , i], u[-index])
    out[, , i] = pl[index, , i] %*% mat
  }
  return(out)
}

#' HPD intervals from a bfa object
#' @param obj A bfa object
#' @param prob Target probability content (see ?HPDinterval)
#' @param loadings Compute interval for the loadings 
#' @param scores Compute interval for the scores
#' @param ... Ignored
#' @return A list with elements loadings.lower, loadings.upper, scores.lower, scores.upper which 
#' are matrices of dimension p x k or n x k, or NA's if loadings or scores is FALSE
#' @method HPDinterval bfa
#' @export

HPDinterval.bfa = function(obj, prob=0.95, loadings=TRUE, scores=FALSE, ...) {
  load.lo = load.hi = score.lo = score.hi = NA
  p = obj$P
  k = obj$K
  n = obj$N
  out=list()
  if(loadings) {
    co = get_coda(obj, loadings, scores=FALSE, scale = attr(obj, "type")!="gauss")
    inter = HPDinterval(co)
    load.lo = matrix(inter[,1], nrow=p)
    load.hi = matrix(inter[,2], nrow=p)
    rownames(load.lo) = rownames(load.hi) = obj$varlabel
    colnames(load.lo) = colnames(load.hi) = paste("Factor",1:k, sep='')
    out$loadings.lower = data.frame(load.lo)
    rownames(out$loadings.lower)=obj$varlabel
    out$loadings.upper = data.frame(load.hi)
    rownames(out$loadings.upper)=obj$varlabel
  }
  if(scores) {
    co = get_coda(obj, loadings=FALSE, scores, scale = attr(obj, "type")!="gauss")
    inter = HPDinterval(co)
    score.lo = matrix(inter[,1], nrow=n)
    score.hi = matrix(inter[,2], nrow=n)
    rownames(score.lo) = rownames(score.hi) = obj$obslabel
    colnames(score.lo) = colnames(score.hi) = paste("Factor",1:k, sep='')
    out$scores.lower = data.frame(score.lo)
    rownames(out$scores.lower)=obj$obslabel
    out$scores.upper = data.frame(score.hi)
    rownames(out$scores.upper)=obj$obslabel
  }
  for(m in out) {
    colnames(m)=paste("Factor",1:obj$K, sep='')
  }
  return(out)
}
