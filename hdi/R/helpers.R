lasso.cv <- function(x, y, nfolds = 10, grouped = nrow(x) > 3 * nfolds,...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 25 Mar 2013, 17:08
  
  fit.cv <- cv.glmnet(x, y, nfolds = nfolds, grouped = grouped, ...)
  ## Use default value of "lambda1.se" in cv.glmnet optimal lambda sel.
  sel <- predict(fit.cv, type = "nonzero") ## Intercept??? Exceptions???
  sel[[1]] ## ugly...
}

lasso.firstq <- function(x, y, q, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  2 Apr 2013, 13:42

  ## Use glmnet (dfmax = q+1 because of Intercept)
  fit <- glmnet(x, y, dfmax = q, ...)   ## only need partial path
  m   <- predict(fit, type = "nonzero") ## determine non-zero coefs

  ## determine largest model that is <= q
  delta <- q - unlist(lapply(m, length)) ## deviation from desired model size
  delta[delta < 0] <- Inf ## overshooting not allowed

  take <- which.min(delta) ## takes first occurrence
  m[[take]]
}

lm.pval <- function(x, y, exact = TRUE, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  2 Apr 2013, 11:34

  fit.lm <- lm(y ~ x, ...) ## Intercept??? Exceptions???
  fit.summary <- summary(fit.lm)

  tstat <- coef(fit.summary)[-1, "t value"] ## Intercept??? Exceptions???

  ## return p-values
  setNames(2 * (if(exact) ## Use appropriate t-dist
                  pt(abs(tstat), df = fit.lm$df.residual, lower.tail = FALSE)
                else  ## p-values based on *normal* distribution
                  pnorm(abs(tstat), lower.tail = FALSE)),
           colnames(x))
}

lm.ci <- function(x, y, level = 0.95, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 20 Feb 2014, 13:43

  fit.lm <- lm(y ~ x, ...) ## Intercept??? Exceptions???

  confint(fit.lm, level = level)[-1,, drop = FALSE]
}

glm.pval <- function(x, y, family = "binomial", verbose = FALSE, ...)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Ruben Dezeure based on lm.pval, Date:  30 Sept 2013, 18:04

    fit.glm <- glm(y ~ x, family = family, ...) ## Intercept??? Exceptions???
    fit.summary <- summary(fit.glm)

    if(!fit.glm$converged & verbose){ ## should be consistent with lm.pval?
      print(fit.summary)
    }

    pval.sel <- coef(fit.summary)[-1,4] ## dangerous with [,4]???

##-     if(family %in% c("poisson", "binomial")){
##-       zstat <- fit.summary$coefficients[-1, "z value"]
##-         ## Intercept??? Exceptions???
##-
##-       ## p-values based on *normal* distribution
##-       pval.sel <- 2 * pnorm(abs(zstat), lower.tail = FALSE)
##-     }else{
##-       tstat <- fit.summary$coefficients[-1, "t value"]
##-       ## Intercept??? Exceptions???
##-       pval.sel <- 2 * pt(abs(tstat), df = fit.lm$df.residual,
##-                          lower.tail = FALSE)
##-     }
    names(pval.sel) <- colnames(x)
    pval.sel
}

fdr.adjust <- function(p)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 17 Jul 2013, 16:42

  p.fin <- p
  use <- (p < 1)
  if(any(use)){
    p.use <- p[use]

    lp <- length(p.use) ## as in p.adjust
    i <- lp:1L
    o <- order(p.use, decreasing = TRUE)
    ro <- order(o)
    p.fin[use] <- pmin(1, cummin(p.use[o] / i))[ro]
  }
  p.fin
}

calc.ci <- function(bj, se, level = 0.95)
{
  ## Purpose:
  ## calculating confidence intervals with the given coefficients, standard
  ## errors and significance level.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Feb 2014, 14:27

  quant <- qnorm(1 - (1 - level) / 2)
  list(lci = bj - se * quant,
       rci = bj + se * quant)
}


p.adjust.wy <- function(cov, pval, N = 10000)
{
  ## Purpose:
  ## multiple testing correction with a Westfall young-like procedure as
  ## in ridge projection method, http://arxiv.org/abs/1202.1377 P.Buehlmann
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## cov: covariance matrix of your estimator
  ## pval: the single testing pvalues
  ## N: the number of samples to take for the empirical distribution
  ##    which is used to correct the pvalues
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Feb 2014, 14:27

  ## Simulate distribution
  zz  <- mvrnorm(N, rep(0, ncol(cov)), cov)
  zz2 <- scale(zz, center = FALSE, scale = sqrt(diag(cov)))
  Gz  <- apply(2 * pnorm(abs(zz2),lower.tail = FALSE), 1, min)

  ## Corrected p-values  pcorr
  ecdf(Gz)(pval)
}

preprocess.group.testing <- function(N, cov, conservative)
{
  if(conservative){
    ## No preprocessing to perform
    NULL
  } else {
    ## Simulate distribution
    zz  <- mvrnorm(N, rep(0, ncol(cov)), cov)
    scale(zz, center = FALSE, scale = sqrt(diag(cov)))
  }
}

calculate.pvalue.for.group <- function(brescaled, group, individual,
                                       Delta = NULL, conservative = TRUE, zz2)
{
  ## Purpose:
  ## calculation of p-values for groups
  ## using the maximum as test statistic
  ## http://arxiv.org/abs/1202.1377 P.Buehlmann
  ## Author: Ruben Dezeure 2 May 2014
  if(is.list(group)){
    pvalue <- lapply(group,
                     calculate.pvalue.for.group,
                     brescaled = brescaled,
                     individual = individual,
                     Delta = Delta,
                     conservative = conservative,
                     zz2 = zz2)
    pvalue <- unlist(pvalue)
  }else{
    p <- length(brescaled)

    if(!is.logical(group)){
      stopifnot(all(group <= p) & all(group >= 1))
      tmp <- logical(length(brescaled))
      tmp[group] <- TRUE
      group <- tmp
    }

    stopifnot(is.logical(group))
    stopifnot(length(group) == length(brescaled))

    if(conservative){ ## (very) conservative alternative proposed by Nicolai
      pvalue <- min(individual[group])
      ## if(correct){ ## only for alternative method
      pvalue <- min(1, p * pvalue) ## Bonferroni correction
      ##}
    }else{
      ## Max test statistics according to http://arxiv.org/abs/1202.1377
      ## P.Buehlmann
      if(is.null(zz2))
        stop("you need to preprocess zz2 by calling preprocess.group.testing")
      if(sum(group) > 1){
        group.coefficients <- abs(brescaled[group])
        max.coefficient    <- max(group.coefficients)
        group.zz           <- abs(zz2[,group])

        if(!is.null(Delta)){ ## special case Ridge method
          group.Delta <- Delta[group]
          group.zz    <- sweep(group.zz, 2, group.Delta, "+")
        }

        Gz <- apply(group.zz, 1, max)

        ## Determine simulated p-value
        pvalue <- 1 - ecdf(Gz)(max.coefficient)
      }else{ ## if group consists of a single variable
        pvalue <- individual[group]
      }
    }
  }
  return(pvalue)
}

calculate.pvalue.for.cluster <- function(hh, p, pvalfunction,
                                         clusterextractfunction,
                                         alpha, verbose = FALSE) {
  ## Remark: not the cleanest code yet
  ## Purpose:
  ## calculation of p-values hierarchically for a cluster as input
  ## Author: Ruben Dezeure 10 April 2014

  ## test the clusters
  upper.signif     <- rep(TRUE, p) ## the significant clusters one level up
  upper.clust      <- list()       ## the upper significant cluster
  upper.clust[[1]] <- 1:p
  upper.clust.ind  <- list()

  ord <- hh$order## ordering from the cluster

  clusters   <- list()
  pvalue     <- list()
  leftChild  <- list()
  rightChild <- list()
  iter       <- 1
  old.cluster.lengths <- 0
  start <- TRUE
  for(nclust in 1:p){ ## test the subclusters of the significant current clusters
    if(verbose)
      cat("cutting tree into", nclust, "groups\n")

    cut.level <- clusterextractfunction(nclust)## TODO clusterextractfunction
    ## cutree(hh,k=nclust)

    all.signif <- all(upper.signif)
    if(!all.signif){
      cut.level[!upper.signif] <- 0
      ## ignore clusters that are not significant above
    }
    clusters.this.level <- sapply(setdiff(unique(cut.level), 0),
                                  FUN = "==", cut.level)
    cluster.lengths     <- sort(apply(clusters.this.level, 2, sum))

    if(!(identical(cluster.lengths,old.cluster.lengths))){
      ## we went down the tree one level for those clusters we are interested in
      old.cluster.lengths <- cluster.lengths
      if(verbose)
        cat("total number of clusters at this level =",
            ncol(clusters.this.level), "\n")

      stopifnot(ncol(clusters.this.level) > 0)

      current.clusts <- list()
      for(j in 1:ncol(clusters.this.level)){
        current.clusts[[j]] <- which(clusters.this.level[,j])
      }
      ## one of the upper clusters has children
      ## find out which are the children
      are.children <- rep(FALSE,ncol(clusters.this.level))
      has.children <- rep(FALSE,length(upper.clust))
      for(j in 1:length(current.clusts)){## check if current cluster is a child
        for(i in 1:length(upper.clust)){ ## check if it has a child
          if(all(is.element(current.clusts[[j]],upper.clust[[i]])) &&
             (start || !(length(current.clusts[[j]]) ==
                         length(upper.clust[[i]]))))
            { ## current cluster is a child of above
              if(verbose){
                cat("the subcluster of length\n")
                cat(length(current.clusts[[j]]), "\n")
                cat("fits in the cluster of length\n")
                cat(length(upper.clust[[i]]), "\n")
                cat("saving the child for testing\n")
              }
              are.children[j] <- TRUE
              has.children[i] <- TRUE
            }
        }
      }
      if(verbose)
        cat("there are", sum(are.children), "children that we will test.\n")

      if(sum(are.children) > 0){ ## else: we have not encountered children in
                               ##       this level of tree

        ## if there are children it is only possible that there are two
        ## children clusters of
        ## one upper cluster that has children, this is due to the iteration
        ## over cutree
        ## TODO
        if(!start){
          ## add leftchild and rightchild
          ## we assume always that nchildren =2 at this point!
          if(mean(((1:length(ord)))
                  [ord %in%current.clusts[[which(are.children)[2]]]] / p) <
             mean(((1:length(ord)))
                  [ord %in%current.clusts[[which(are.children)[1]]]] / p)){
            left  <- iter + 1
            right <- iter
          }else{
            left  <- iter
            right <- iter + 1
          }
          parent.ind <- upper.clust.ind[[which(has.children)]]
          ## Need to know the index of this in clusters!
          leftChild[[parent.ind]]  <- left
          rightChild[[parent.ind]] <- right
        }
        ## test them
        clusters.to.test <- clusters.this.level[, are.children, drop = FALSE]
        pvals <- pvalfunction(clusters.to.test = clusters.to.test,
                              nclust = nclust)## TODO pvalfunction
        ## mapply(group.testing.function,
        ## group=split(clusters.to.test,col(clusters.to.test)))

        ## save the children and result if it was significant or not
        for(i in 1:sum(are.children)){
          clusters[[iter]]   <- current.clusts[[which(are.children)[i]]]
          leftChild[[iter]]  <- -1  ## initialise on a leaf
          rightChild[[iter]] <- -1 ## initialise on a leaf
          pvalue[[iter]]     <- pvals[i]
          iter <- iter + 1
        }
        ## update the list of clusters that are significant to
        ## keep branching out,
        upper.clust <- upper.clust[!has.children]
        upper.clust.ind <- upper.clust.ind[!has.children]
        tmp.iter <- length(upper.clust)+1
        for(i in which(pvals <= alpha)){
          upper.clust[[tmp.iter]] <- which(clusters.to.test[,i])
          upper.clust.ind[[tmp.iter]] <- iter - sum(are.children) + (i-1)
          ## index in clusters[[]]
          tmp.iter <- tmp.iter + 1
        }
        upper.signif <- 1:p %in% unique(unlist(upper.clust))
        ## done
        if(length(upper.clust) == 0){
          ## there are no more upper significant clusters
          if(verbose){
            cat("reached an end to the cluster tree\n")
            cat("the lowest clusters have sizes:\n")
            cat(apply(clusters.to.test,2,sum), "\n")
          }
          break
        }
      }
    }
    if(start)
      start <- FALSE
  }
  ## return output in the style of lowerbound method
  list(
    clusters   = clusters,
    pval       = unlist(pvalue),
    leftChild  = unlist(leftChild),
    rightChild = unlist(rightChild),
    alpha      = alpha,
    hh         = hh)
}

get.clusterGroupTest.function <- function(group.testing.function, x)
{
  ## Purpose:
  ## facilitate the creation of clusterGrouptest based on the
  ## group.testing.function
  ## Author: Ruben Dezeure 5th of August 2014
  stopifnot(is.function(group.testing.function))

  ## Return  'clusterGroupTest' function :
  function(hcloutput,
           dist = as.dist(1 - abs(cor(x))),
           alpha = 0.05,
           method = "average",
           conservative = TRUE) {
    ## optional argument: hcloutput = the result from a hclust call
    hh <- if(missing(hcloutput))
            hclust(dist, method = method)
          else
            hcloutput

    clusterextractfunction <- function(nclust){
      ## function to extract the correct level of the tree where the number
      ## of clusters = nclust
      cutree(hh, k = nclust)
    }

    pvalfunction <- function(clusters.to.test, nclust){
      ## function to calculate the p-value for certain clusters
      mapply(group.testing.function,
             conservative = conservative,
             group = split(clusters.to.test, col(clusters.to.test)))
    }
    structure(c(calculate.pvalue.for.cluster(hh = hh,
                                           p = ncol(x),
                                           pvalfunction = pvalfunction,
                                           alpha = alpha,
                                           clusterextractfunction =
                                             clusterextractfunction),
                method = "clusterGroupTest"),
              class = c("clusterGroupTest", "hdi"))
  }
}

switch.family <- function(x, y, family)
{
  switch(family,
         "binomial" = {
           fitnet        <- cv.glmnet(x, y, family = "binomial",
                                      standardize = FALSE)
           glmnetfit     <- fitnet$glmnet.fit
           netlambda.min <- fitnet$lambda.min
           netpred <- predict(glmnetfit, x, s = netlambda.min,
                              type = "response")
           betahat <- predict(glmnetfit, x, s = netlambda.min,
                              type = "coefficients")
           betahat <- as.vector(betahat)
           pihat   <- netpred[,1]

           diagW <- pihat * (1 - pihat)
           W     <- diag(diagW)
           xl    <- cbind(rep(1, nrow(x)), x)

           ## Adjusted design matrix
           xw <- sqrt(diagW) * x

           ## Adjusted response
           yw <- sqrt(diagW) * (xl %*% betahat + solve(W, y - pihat))
         },
         {
           stop("The provided family is not supported (yet). Currently supported are gaussian and binomial.")
         })
  return(list(x = xw, y = yw))
}


sandwich.var.est.stderr <- function(x,y,betainit,Z){
  ## Purpose:
  ## an implementation of the calculation of the robust standard error
  ## based on the sandwich variance estimator from
  ## http://arxiv.org/abs/1503.06426
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x: the design matrix
  ## y: the response vector
  ## betainit: the initial estimate
  ## Z:       the residuals of the nodewis regressions
  ## ----------------------------------------------------------------------
  ## Return values:
  ## se.robust: the robust estimate for the standard error of the corresponding de-sparsified lasso fit
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 18 Mai 2015 (initial version)

  n <- nrow(x)
  p <- ncol(x)

  ## Check if normalization is fulfilled
  if(!isTRUE(all.equal(rep(1, p), colSums(Z * x) / n, tolerance = 10^-8))){
    ## no need to print stuff to the user, this is only an internal detail
    rescale.out <- score.rescale(Z = Z, x = x)
    Z <- rescale.out$Z
    ## scaleZ <- rescale.out$scaleZ
  }

  if(length(betainit) > ncol(x)){
    x <- cbind(rep(1,nrow(x)),x)
  }else{
    if(all(x[,1]==1)){
      ## ok, we have included the intercept
    }else{
      ## hmm, should we substract the mean from y?
    }
  }
  eps.tmp <- as.vector(y - x%*%betainit)

  ## should these esp.tmp be forced to have mean 0?

  sigmahatZ.direct <- sqrt(colSums(sweep(eps.tmp*Z, MARGIN = 2,
                                         STATS = crossprod(eps.tmp,Z)/n,
                                         FUN = `-`)^2))
  ## rm(eps.tmp)
  ## return se.robust :
  sigmahatZ.direct/n ## this is the s.e. of bproj from pval.score,
  ## if we multiply bproj with 1/this, we get on the N(0,1) scale
}

do.initial.fit <- function(x,y,
                           initial.lasso.method = c("scaled lasso","cv lasso"),
                           lambda,
                           verbose = FALSE)
{
  ## Purpose:
  ## This function performs the initial fit of the high dimensional linear model
  ## used in ridge.proj and lasso.proj
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x: the design matrix
  ## y: the response vector
  ## initial.lasso.method: the method to use to tune the lasso
  ## lambda: OPTIONAL, the tuning parameter for the lasso. If this is provided,
  ##         it overrides initial.lasso.method
  ## verbose: let's you know if the user is overriding initial.lasso.method with a value for lambda
  ## ----------------------------------------------------------------------
  ## Return values:
  ## betalasso: the lasso coefficients, excluding the intercept
  ## sigmahat: estimate for the noise standard  deviation
  ## intercept: intercept, in case the lasso was fitted with it
  ## lambda: the tuning parameter for the lasso that was used, not available for the scaled lasso
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 16 Oct 2015 (initial modified version for the package)

  no.lambda.given <- missing(lambda) || is.null(lambda)
  if(!no.lambda.given && verbose) {
    cat("A value for lambda was provided to the do.initial.fit function,\n")
    cat("the initial.lasso.method option was therefore ignored.\n")
    cat("We now do a lasso with the tuning parameter instead of a self-tuning procedure.\n") 
  }

  if(no.lambda.given){
    ## tune the lasso
    switch(initial.lasso.method,
           "scaled lasso"={
             scaledlassofit <- scalreg(X=x,y=y)
             lambda <- NULL## no way to extract lambda for this ?
           },
           "cv lasso"={
             glmnetfit <- cv.glmnet(x=x,y=y)
             lambda <- glmnetfit$lambda.1se
           },
           {
             stop("Not sure what lasso.method you want me to use for the initial fit. The only options for the moment are: 1)scaled lasso 2)cvlasso")
           })
  }else{
    ## fit for a range of lambda
    glmnetfit <- glmnet(x=x,y=y)
  }

  if(no.lambda.given && identical(initial.lasso.method,"scaled lasso")){
    intercept <- 0
    betalasso <- scaledlassofit$coefficients
    sigmahat  <- scaledlassofit$hsigma
    residual.vector <- y-x%*%betalasso
  }else{
    if((nrow(x) - sum(as.vector(coef(glmnetfit, s = lambda)) != 0)) <= 0){
      ## Problem: when the fixed lambda you chose sets n==p, you've used all
      ## your degrees of freedom 
      ##- -> refit
      ## This only occurs if the user provides a lambda to use
      message("Refitting using cross validation: your lambda used all degrees of freedom")
      glmnetfit <- cv.glmnet(x=x,y=y)
      ## setting a lambda here would interpolate solutions, :/
      lambda <- glmnetfit$lambda.1se
    }

    intercept <- coef(glmnetfit,s=lambda)[1]
    betalasso <- as.vector(coef(glmnetfit,s=lambda))[-1]## leaving out the intercept
    residual.vector <- y-predict(glmnetfit,newx=x,s=lambda)
    sigmahat <- sqrt(sum((residual.vector)^2)/
                     (nrow(x)-sum(as.vector(coef(glmnetfit,s=lambda))!=0)))
  }
  ## return
  list(betalasso=betalasso,
       sigmahat=sigmahat,
       intercept=intercept,
       lambda=lambda)
}

initial.estimator <- function(betainit,x,y,sigma)
{
  ## check if betainit is correctly provided
  if(!((is.numeric(betainit) && length(betainit) == ncol(x)) ||
       (betainit %in% c("scaled lasso","cv lasso"))))
  {
    stop("The betainit argument needs to be either a vector of length ncol(x) or one of 'scaled lasso' or 'cv lasso'")
  }

  warning.sigma.message <- function() {
    warning("Overriding the error variance estimate with your own value. The initial estimate implies an error variance estimate and if they don't correspond the testing might not be correct anymore.") 
  }

  if(is.numeric(betainit))
  {
    beta.lasso <- betainit

    if(is.null(sigma))
    {
      stop("Not sure what variance estimate to use here")
      ## if betainit comes from the lasso, the below should be good
      residual.vector <- y- x%*%betainit
      sigmahat <- sqrt(sum((residual.vector)^2)/
                       (nrow(x)-sum(betainit!=0)))
    }else{
      warning.sigma.message()
      sigmahat <- sigma
    }
  }else{
    initial.fit <- do.initial.fit(x=x,y=y,
                                  initial.lasso.method=betainit)
    beta.lasso <- initial.fit$betalasso

    if(is.null(sigma))
    {
      sigmahat <- initial.fit$sigmahat
    }else{
      warning.sigma.message()
      sigmahat <- sigma
    }
  }
  ## return
  list(beta.lasso = beta.lasso,
       sigmahat = sigmahat)
}
