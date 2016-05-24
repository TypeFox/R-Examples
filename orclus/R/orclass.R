################################
### definition of orclass method:
################################

orclass <- function (x, ...) 
    UseMethod("orclass")

###################################
### main orclass.default - function:
###################################

orclass.default <- function(x, grouping, k, l, k0, a = 0.5, prior = NULL, inner.loops = 1, predict.train = "nearest", verbose = TRUE, ...){
  # x:            data input matrix
  # grouping      vector of class labels
  # k:            number of clusters
  # l:            final subspace dimension
  # k0:           number of initial clusters
  # a:            reduction rate for subspace dimension over the iterations
  # prior         class prior probabilities
  # inner.loops:  number of iterations for each pair of 'number of clusters' and 'subspace diension' during the iteration process
  # predict.train: "nearest" / "fuzzywts" (type) if prediction of training data is desired
  # verbose:      logical specifying whether iteration process should be printed

  # compute (relative) class frequencies and priors
  relfreqs  <- table(grouping)
  relpriors <- relfreqs / sum(relfreqs)
  if (is.null(prior)) prior <- relpriors
  if(length(prior) != length(relfreqs)) stop("Length of prior does not match number of classes in grouping.")
  if(any(prior == 0)) stop("Any priors of 0 or empty classes!")
  
  # compute orclus model
  orclus.res <- orclus(x=x, k=k, l=l, k0=k0, a = a, inner.loops=inner.loops, verbose=verbose)# ,...)

  ### compute cluster.posteriors
  xtab <- table(orclus.res$cluster, grouping)
  
  for(i in 1: length(prior)) xtab[,i] <- xtab[,i] * prior[i] / relpriors[i]
  cluster.posteriors <- xtab / rowSums(xtab) 
  
  # compute cluster frequencies
  cluster.freqs <- table(orclus.res$cluster)
  rel.clusfreqs <- cluster.freqs/sum(cluster.freqs) 
  # prior correction of cluster frequencies
  prior.corrected.cluster.freqs <- rowSums(xtab)
  rel.pricorr.clusfreqs <- prior.corrected.cluster.freqs/sum(prior.corrected.cluster.freqs)
  
  ### compute criteria of discriminatory power
  # compute train.error
  train.error <- sum(apply(cluster.posteriors, 1, function(z) return(1-max(z))) *  rel.pricorr.clusfreqs)
  # compute entropy 
  cl.entros <- apply(cluster.posteriors, 1, function(z) return( -sum(log2(z[which(z !=0)]) * z[which(z !=0)]) ))
  entro <- sum(cl.entros * rel.pricorr.clusfreqs)
  # compute gini
  gini <- sum((1-rowSums(cluster.posteriors^2)) *  rel.pricorr.clusfreqs)
  # compute chisq
  chisq <- chisq.test(xtab)$p.value # note: xtab is already prior corrected
  
  # summarize criteria in one vector
  sumcrit <- c(train.error = train.error, entropy = entro, gini = gini, chisquare.p = chisq)
  
  # call
  cl <- match.call()
  cl[[1]] <- as.name("orclass")
  
  # create result.list
  orclass.res <- list(orclus.res = orclus.res, cluster.posteriors = cluster.posteriors 
                      , cluster.priors = rel.pricorr.clusfreqs, purity = sumcrit 
                      , prior = prior, predict.train = NULL, orclass.call = cl
                      )
  
  # create result-object of class orclus
  class(orclass.res) <- "orclass"
  
  # training prediction, if specified
  if (predict.train == "nearest") orclass.res$predict.train <- predict(orclass.res, x, type = "nearest")
  if (predict.train == "fuzzywts") orclass.res$predict.train <- predict(orclass.res, x, type = "fuzzywts")
  
  return(orclass.res)   
  }


##################################
### orclass formula interface
##################################

orclass.formula <- function(formula, data = NULL, ...)
{
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval.parent(m$data))) 
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  grouping <- model.response(m)
  x <- model.matrix(Terms, m)
  #x <- model.frame(Terms, m)
  xvars <- as.character(attr(Terms, "variables"))[-1]
  if ((yvar <- attr(Terms, "response")) > 0) 
    xvars <- xvars[-yvar]
  xlev <- if (length(xvars) > 0) {
    xlev <- lapply(m[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  xint <- match("(Intercept)", colnames(x), nomatch = 0)
  if (xint > 0) 
    x <- x[, -xint, drop = FALSE]
  res <- orclass(x, grouping, ...)
  res$terms <- Terms
  res$contrasts <- attr(x, "contrasts")
  res$xlevels <- xlev
  attr(res, "na.message") <- attr(m, "na.message")
  if (!is.null(attr(m, "na.action"))) 
    res$na.action <- attr(m, "na.action")
  res
}



########################
### predict function
########################

predict.orclass <- function(object, newdata, type = "nearest", ...){
  # several error checks 
  if (class(object) != "orclass") stop("Object must be of class orclass!")
  if(is.data.frame(newdata)) newdata <- as.matrix(newdata)
  if(!is.matrix(newdata)) stop("Newdata must be of type matrix or data frame!")
  if(dim(object$orclus.res$centers)[2] != dim(newdata)[2]) stop("Newdata must be of same dimension as the explanatory input variables data set!")
  
  predicted.clusters <- predict(object$orclus.res, newdata)
  
  # assign posteriors
  if (type == "nearest"){ 
    posteriors <- object$cluster.posteriors[predicted.clusters$cluster,]
    }
  
  if (type == "fuzzywts"){
     # compute cluster membership values (for m=2) [coloumn wise]
     funx <- function(x) return(1/sum((x[1] * x[2:length(x)])^2))
     i <- 1
     dummy <- cbind(predicted.clusters$distances[,i], 1/predicted.clusters$distances)
     memberships.i <- apply(dummy, 1, funx)
     memberships <- memberships.i
     for(i in 2:ncol(predicted.clusters$distances)){
       dummy <- cbind(predicted.clusters$distances[,i], 1/predicted.clusters$distances)
       memberships.i <- apply(dummy, 1, funx)
       memberships <- cbind(memberships, memberships.i)       
      }
    # weight clusterwise posteriors with cluster memberships
    posteriors <- memberships %*% object$cluster.posteriors
    } 
  
  # assign classes
  rownames(posteriors) <- rownames(newdata)
  # without randomization (use first in case of multiple maxima): classes <- colnames(posteriors)[apply(posteriors, 1, which.max)]
  # randomization in case of several maxima as which.max returns first maximum
  class.obs <- function(z){res <- which(z == max(z)); if (length(res)>1) res <- sample(res,1); return(res)} 
  classes <- colnames(posteriors)[apply(posteriors, 1, class.obs)]
  
  # create results object
  res <- list()
  res$class       <- classes
  res$posterior   <- posteriors
  res$cluster     <- predicted.clusters$cluster
  res$distances   <- predicted.clusters$distances
  return(res)
  }



########################
### print function
########################

print.orclass <- function(x, ...){
  cat("\nOrclass results\n\n")
  cat("Call:\n")
  print(x$orclass.call)

  cat("\nDimension of subspaces:", x$orclus.res$subspace.dimension, "\n")
  cat("\nNumber of subgroups   :", nrow(x$orclus.res$center), "\n")
  
  cat("\nCluster priors        :", x$cluster.priors, "\n")
  cat("\nClass priors          :", x$prior, "\n")
  cat("\nCluster posteriors    :\n")
  print(x$cluster.posteriors)
  
  cat("\nCluster sparsity coefficient: ", x$orclus.res$sparsity.coefficient, "\n")
  cat("\nCluster purity :\n")
  print(x$purity)
  cat("\n\n")
  
  invisible(x)
  }

