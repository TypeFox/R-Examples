rvalues <- function(data, family = gaussian, hypers = "estimate",
                    prior = "conjugate", alpha.grid = NULL, ngrid = NULL,
                    smooth = "none", control = list())  {
  
  mod <- match.call(expand.dots=FALSE)
  if(is.null(alpha.grid)) {
    ### initialize alpha_grid if not entered by user
    alpha.grid <- MakeGrid(nunits = nrow(data), type = "log", ngrid = ngrid)
  }
  if(!is.null(alpha.grid)) {
    ngrid <- length(alpha.grid)
    
    # check that alpha.grid values are strictly between 0 and 1  
    alpha.check <- all((alpha.grid > 0) & (alpha.grid < 1))
    if(!alpha.check) {
        stop("All the alpha.grid values must be strictly between zero and one")
    }
  }
  
  #### switch based on family argument
  if(class(family)=="newfam") {
    family <- family
  }
  else {
    if (is.character(family))   {
      family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family))  { 
      family <- family()
    }
    if (is.null(family$family)) {
      stop("family argument is not valid")
    }
  }
 
  ###  setup control parameters
  con <- list(tol = 1e-04, maxiter = 500, smooth = FALSE)
  con[(namc <- names(control))] <- control
  
  estimate <- data[,1]
  nuisance <- data[,2]
  switch(family$family,
         gaussian={
             ### Here, nuisance = std_err
             
           if(prior=="conjugate") {
             ###  Use Normal prior
             estimate <- data[,1]
             nuisance <- data[,2]
             res <- rvalue.agrid.nn(estimate,nuisance,hypers,alpha.grid,smooth)
           }
           else {
             ###  Use nonparametric prior
             npfit <- npmle(data,family=gaussian,maxiter=con$maxiter,
                            tol=con$tol,smooth=con$smooth)
             ###  posterior mean
             PM <- npfit$post.mean
             
             lb <- min(npfit$support) - .01
             theta.alpha <- ThetaQuantiles(npfit$Fhat, alpha.grid, lbd = lb, 
                                          ubd = max(npfit$support))
             
             theta.probs <- -diff(c(1,npfit$Fhat(theta.alpha)))
             tmp <- NPagrid(estimate = data[,1], nuisance = data[,2], theta.alpha, 
                            theta.probs, alpha.grid, smooth, family = "gaussian")
             
             bar <- data.frame( RValue=tmp$rvalues, RV.rank=rank(tmp$rvalues),MLE.rank=rank(-estimate),
                               PM.rank=rank(-PM),MLE=estimate,SE=nuisance,PostMean=PM,
                               PVal = pnorm( estimate/nuisance, lower.tail=FALSE ),
                               PVal.rank = rank( -estimate/nuisance ) )  
  
             ord <- order(tmp$rvalues, -PM )
             sorted.dframe <- bar[ord,]
             
             res <- list()
             res$main <- sorted.dframe
             res$aux <- list(mix.prop=npfit$mix.prop,support=npfit$support,mixcdf=npfit$Fhat,
                             alpha.grid=tmp$alpha,Vmarginals=tmp$lam,Vmarginals.smooth=tmp$lamsmooth,
                             V=tmp$V,unsorted=bar,prior="nonparametric",family="gaussian",
                             smooth=smooth, conv=npfit$conv)
              res$rvalues <- tmp$rvalues
           }
         },
         poisson={
           if(prior=="conjugate") {
             estimate <- data[,1]
             nuisance <- data[,2]
             res <- rvalue.agrid.pg(estimate,nuisance,hypers,alpha.grid,smooth)
           }
           else {
             npfit <- npmle(data,family=poisson,maxiter=con$maxiter,
                            tol=con$tol,smooth=con$smooth)
             ### should this be on the log scale?
             PM <- npfit$post.mean
             
             lb <- exp(log(min(npfit$support)) - 2)
             theta.alpha <- ThetaQuantiles(npfit$Fhat, alpha.grid, lbd = lb, 
                                          ubd = max(npfit$support))
             theta.probs <- -diff(c(1,npfit$Fhat(theta.alpha)))
             tmp <- NPagrid(estimate = data[,1], nuisance = data[,2], theta.alpha, 
                            theta.probs, alpha.grid, smooth, family = "poisson")
             bar <- data.frame( RValue=tmp$rvalues, RV.rank=rank(tmp$rvalues),
                                MLE.rank=rank(-estimate/nuisance), PM.rank=rank(-PM),
                                xx = estimate, eta = nuisance, PostMean = PM)  
  
             ord <- order(tmp$rvalues, -PM )
             sorted.dframe <- bar[ord,]
             
             res <- list()
             res$main <- sorted.dframe
             res$aux <- list(mix.prop=npfit$mix.prop,support=npfit$support,mixcdf=npfit$Fhat,
                             alpha.grid=tmp$alpha,Vmarginals=tmp$lam,Vmarginals.smooth=tmp$lamsmooth,
                             V=tmp$V,unsorted=bar,prior="nonparametric",family="poisson",
                             smooth=smooth, conv=npfit$conv)
             res$rvalues <- tmp$rvalues
           }
         },
         binomial={
           if(prior=="conjugate") {
             estimate <- data[,1]
             nuisance <- data[,2]
             res <- rvalue.agrid.bb(estimate,nuisance,hypers,alpha.grid,smooth)
           }
           else if(prior=="nonparametric") {
             npfit <- npmle(data,family=binomial,maxiter=con$maxiter,
                            tol=con$tol,smooth=con$smooth)
             #PM <- npmixapply(npfit,function(x) { x })
             PM <- npfit$post.mean
             
             min.sup <- log(min(npfit$support)) - log(1 - min(npfit$support))
             lb <- exp(min.sup - 2)/( 1 + exp(min.sup - 2))
             
             theta.alpha <- ThetaQuantiles(npfit$Fhat, alpha.grid, lbd = lb, 
                                          ubd = max(npfit$support))
             
             theta.probs <- -diff(c(1,npfit$Fhat(theta.alpha)))
             tmp <- NPagrid(estimate = data[,1], nuisance = data[,2], theta.alpha, 
                            theta.probs, alpha.grid, smooth, family = "binomial")
                            
             mle <- rep(0,length(estimate))
             mle[nuisance > 0] <- estimate[nuisance > 0]/nuisance[nuisance > 0]
             mle[nuisance == 0] <- 0 
             bar <- data.frame(RValue=tmp$rvalues, RV.rank=rank(tmp$rvalues),
                               MLE.rank=rank(-mle), PM.rank=rank(-PM),xx=estimate,
                               nn=nuisance,PostMean=PM )  
  
             ord <- order(tmp$rvalues, -PM )
             sorted.dframe <- bar[ord,]
             
             res <- list()
             res$main <- sorted.dframe
             res$aux <- list(mix.prop=npfit$mix.prop,support=npfit$support,mixcdf=npfit$Fhat,
                             alpha.grid=tmp$alpha,Vmarginals=tmp$lam,Vmarginals.smooth=tmp$lamsmooth,
                             V=tmp$V,unsorted=bar,prior="nonparametric",family="binomial",
                             smooth=smooth, conv=npfit$conv)
             res$rvalues <- tmp$rvalues 
           }
         },
         Gamma={
             ### assumes that X|theta,shapes ~ Gamma(shapes, theta) (scale-form)
             ### and theta ~ InvGamma(a, b)
             ### This implies that theta|X,shapes ~ InvGamma(shapes + a, X + b)
             estimate <- data[,1]
             nuisance <- data[,2]
             res <- rvalue.agrid.gg(estimate, nuisance, hypers, alpha.grid, smooth)
         },
         tdist={
            if(prior=="conjugate") {
                stop("Must use nonparametric prior with the t family")
            }
            else if(prior=="nonparametric") {
                npfit <- npmle(data,family=family,maxiter=con$maxiter,
                            tol=con$tol,smooth=con$smooth)
                #PM <- npmixapply(npfit,function(x) { x })
                PM <- npfit$post.mean
                
                lb <- min(npfit$support) - .01
                theta.alpha <- ThetaQuantiles(npfit$Fhat, alpha.grid, lbd = lb, 
                                               ubd = max(npfit$support))
             
                theta.probs <- -diff(c(1,npfit$Fhat(theta.alpha)))
                tmp <- NPagrid(estimate = data[,1], nuisance = data[,2], theta.alpha, 
                               theta.probs, alpha.grid, smooth, family = "tdist", df=family$df)
                                
                bar <- data.frame( RValue=tmp$rvalues, RV.rank=rank(tmp$rvalues),MLE.rank=rank(-estimate),
                               PM.rank=rank(-PM),MLE=estimate,SE=nuisance,PostMean=PM,
                               PVal = pt( estimate/nuisance, df=family$df, lower.tail=FALSE),
                               PVal.rank = rank( -estimate/nuisance ) )  
  
                 ord <- order(tmp$rvalues, -PM )
                 sorted.dframe <- bar[ord,]
             
                 res <- list()
                 res$main <- sorted.dframe
                 res$aux <- list(mix.prop=npfit$mix.prop,support=npfit$support,
                             mixcdf=npfit$Fhat,alpha.grid=tmp$alpha,
                             Vmarginals=tmp$lam,Vmarginals.smooth=tmp$lamsmooth,
                             V=tmp$V,unsorted=bar,prior="nonparametric",
                             family="tdist", df=family$df, smooth = smooth, conv=npfit$conv)
                 res$rvalues <- tmp$rvalues
            }
         },
  )
  ans <- list()
  class(ans) <- "rvals"
  ans$main <- res$main
  ans$aux <- res$aux
  ans$rvalues <- res$rvalues
  ans$call <- mod
  return(ans)
}

################################################################################
#############  Print function  ########################################
################################################################################

print.rvals <- function(x, ...) {
  ### print method for r-value objects.
  
  nnr <- min(nrow(x$main),10)
  print(x$main[1:nnr,])
}


################################################################################
######  Functions Calling C  ###################################################
################################################################################


VVcut <- function(Vmat, lam_fun, nunits, ngrid, alpha.grid) {
   res <- .Call("Vcut",Vmat,lam_fun,as.integer(nunits),as.integer(ngrid),alpha.grid,PACKAGE="rvalues");
   return(res) 
}

