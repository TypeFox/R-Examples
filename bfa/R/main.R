
.bfa <- function(x, data=NULL, num.factor=1, restrict=NA, normal.dist=NA, 
                center.data=TRUE, scale.data=TRUE, nsim=0, nburn=0, thin=1,
                print.status=500, keep.scores=FALSE, keep.loadings=TRUE,
                loading.prior="pointmass", factor.scales=FALSE, px=TRUE,
                coda="loadings", coda.scale=TRUE, imh=FALSE, imh.iter=500,
                imh.burn=500, ...) {
  if (class(x) != 'bfa') {
    model = bfa_model(x, data, num.factor, restrict, normal.dist, center.data,
                    scale.data, ...)
  } else {
    model=x
  }
  
  if(imh) {
    cat("Building proposal...")
    model = fit_bfa(model, imh.iter, imh.burn, 1, print.status,
                   keep.scores, keep.loadings, loading.prior,
                   factor.scales, px, coda, coda.scale, save_max=TRUE, quiet=TRUE, ...)

    P = model$P
    means <- vector("list", P)
    covs  <- vector("list", P)
    precs <- vector("list", P)
    
    levs = rep(0, P)
    
    for (ind in 1:P) {
      ind.m = which(model$maxes[ind,]==99999.0) - 1
      levs[ind] = ind.m
      if (model$error_var_i[ind]>0) {
        prop.cov = matrix(1)
        prop.prec = matrix(1)
        prop.mean = matrix(1)
      } else if (ind.m > 1) {
        ps = t(model$post.cutpoints[ind,1:ind.m,])
        sample.cov = 2.4*2.4*cov(as.matrix(ps))/ind.m
        prop.cov = chol(sample.cov)
        prop.prec = solve(sample.cov)
        prop.mean = apply(as.matrix(ps), 2, mean)
      } else {
        ps = t(model$post.cutpoints[ind,1:ind.m,])
        
        prop.cov = sqrt(var(c(ps)))
        prop.prec = 1/prop.cov
        prop.mean = mean(ps)
      }
      covs[[ind]] = prop.cov
      precs[[ind]] = prop.prec
      means[[ind]] = prop.mean
    }
    cat("Done.\n")
    model = fit_bfa(model, nsim, nburn, thin, print.status,
                   keep.scores, keep.loadings, loading.prior,
                   factor.scales, px, coda, coda.scale, imh=1, 
                   prop.cov=covs, prop.prec=precs, prop.mean=means, 
                   nlevels=matrix(levs), df=100,  ...)
  } else {
    model = fit_bfa(model, nsim, nburn, thin, print.status,
                   keep.scores, keep.loadings, loading.prior,
                   factor.scales, px, coda, coda.scale,  ...)
  }
  return(model)

}


# Perform MCMC model fitting for a bfa model
#
# This function performs a specified number of MCMC iterations and
# returns an object containing summary statistics from the MCMC samples
# as well as the actual samples if keep.scores or keep.loadings are \code{TRUE}.
# Default behavior is to save only the loadings. 
# 
# Prior parameters:
# loadings.var: Factor loading prior variance
# tau.a, tau.b: Gamma hyperparameters (scale=1/b) for factor precisions (if factor.scales=T)
# rho.a, rho.b: Beta hyperparameters for point mass prior
# sigma2.a, sigma2.b: Gamma hyperparameters for error precisions
# gdp.alpha, gdp.beta: GDP prior parameters
# 
# @param model an object of type bfa, as returned by bfa(data)
# @param nsim number of iterations past burn-in
# @param nburn number of initial (burn-in) iterations to discard
# @param thin keep every thin'th MCMC sample (i.e. save nsim/thin samples)
# @param print.status how often to print status messages to console
# @param keep.scores save samples of factor scores
# @param keep.loadings save samples of factor loadings
# @param loading.prior Specify point mass ("pointmass", default) or normal priors ("normal") 
# @param factor.scales Include a separate scale parameter for each factor
# @param px Use parameter expansion for ordinal/copula/mixed factor models (recommended)
# @param coda create \code{mcmc} objects to allow use of functions from the 
# \code{coda} package: "all" for loadings and scores, "loadings" or "scores" for one or the
# other, or "none" for neither
# @param coda.scale put the loadings on the correlation scale when creating \code{mcmc} objects
# @param ... Prior parameters and other (experimental) arguments (see details)
# @return The S3 \code{bfa} object \code{model}, now with posterior samples/summaries.

fit_bfa <- function(model, nsim, nburn, thin=1, print.status=500,
                   keep.scores=FALSE, keep.loadings=TRUE, loading.prior="pointmass",
                   factor.scales=FALSE, px=TRUE, coda="loadings", coda.scale=TRUE,  ...) {
  
  more_args = list(...)
  if (is.null(more_args$init)) more_args$init=TRUE
  more_args$px = as.numeric(px)
  
  rhoa  = ifelse(is.null(more_args$rho.a),     1.0, more_args$rho.a)
  rhob  = ifelse(is.null(more_args$rho.b),     1.0, more_args$rho.b)
  s     = ifelse(is.null(more_args$s),         1.0, more_args$s)
  taua  = ifelse(is.null(more_args$tau.a),     1.0, more_args$tau.a)
  taub  = ifelse(is.null(more_args$tau.b),     1.0, more_args$tau.b)
  s2a   = ifelse(is.null(more_args$sigma2.a),  0.5, more_args$sigma2.a)
  s2b   = ifelse(is.null(more_args$sigma2.b),  1.0, more_args$sigma2.b)
  alpha = ifelse(is.null(more_args$gdp.alpha), 3.0, more_args$gdp.alpha)
  beta  = ifelse(is.null(more_args$gdp.beta),  1.0, more_args$gdp.beta)
  bp    = ifelse(is.null(more_args$bp),       10.0, more_args$gdp.alpha)
  bq    = ifelse(is.null(more_args$bq),        5.0, more_args$gdp.beta)
  

  model$priors = list(rhoa=rhoa, rhob=rhob, s=s, taua=taua, taub=taub, 
                      s2a=s2a, s2b=s2b, alpha=alpha, beta=beta,
                      bp=bp, bq=bq)
  
  if (model$nsim == 0 && more_args$init==TRUE && more_args$px>0) {
    model = .fit(model, max(1000, nburn/10), 0, thin=1, factor.scales=factor.scales,
                 print.status=max(1000, nburn/2)+1,
                 keep.scores=FALSE, keep.loadings=FALSE, 
                 loading.prior=loading.prior, coda="none", 
                 save_max=0, px=0, more_args$error_var_i, quiet=TRUE )
  }
  
  model = .fit(model, nsim, nburn, thin=thin, print.status=print.status,
                   keep.scores=keep.scores, keep.loadings=keep.loadings, 
                   loading.prior=loading.prior, coda=coda, coda.scale=coda.scale,
                   px=more_args$px, factor.scales=factor.scales, ...)
  
}

########################################################
# Actual model fitting fcn
########################################################

.fit <- function(model, nsim, nburn, thin=1, print.status=500,
                   keep.scores=FALSE, keep.loadings=TRUE, loading.prior="pointmass",
                   coda="loadings", coda.scale=TRUE,  ...) {
      
  K = model$K; P = model$P; N = model$N
  
  more_args = list(...)
  
  more_args$loadings_var  = ifelse(is.null(more_args$loadings.var), 1.0, more_args$loadings.var)
  more_args$factor_scales = ifelse(is.null(more_args$factor.scales), FALSE, more_args$factor.scales)
  
  if (loading.prior=="normal")     more_args$method = 0
  if (loading.prior=="pointmass")  more_args$method = 1
  if (loading.prior=="gdp")        more_args$method = 2
  
  if(is.null(more_args$imh)) more_args$imh = 0
  if(is.null(more_args$px)) more_args$px = 1
  if(is.null(more_args$df)) more_args$df = 100
  if(is.null(more_args$save_max)) more_args$save_max = 0
  if(is.null(more_args$positive)) more_args$positive = FALSE
  
  more_args$save_max = as.numeric(more_args$save_max)
  more_args$error_var_i = model$error_var_i

	model$nsim  = nsim
	model$nburn = nburn
	model$thin  = thin
  
  more_args$mInd = model$mInd

	sim = .MCMCstep(model$ldata, model$loadings, model$scores, model$tauinv, 
		              model$rho, model$ranks, model$maxes, model$argsorts, model$loadings.restrict,
		              model$priors, nsim, nburn, thin, print.status, keep.scores, keep.loadings, more_args)
  
  if (keep.scores)   dim(sim$Fp)=c(K,N,nsim/thin)
  if (keep.loadings) dim(sim$Ap)=c(P,K,nsim/thin)
  if (!is.null(more_args$save_max) && more_args$save_max>0) dim(sim$maxp) = c(dim(model$maxes), nsim/thin)	
 
  dim(sim$Fp.mean) = c(K,N)
  dim(sim$Fp.var)  = c(K,N)
  dim(sim$Ap.mean) = c(P,K)
  dim(sim$Ap.var)  = c(P,K)
  dim(sim$pnz)     = c(P,K)

  model$post.loadings.mean = sim$Ap.mean
  model$post.loadings.var  = sim$Ap.var
  model$post.loadings      = sim$Ap
	
  model$post.scores.mean   = sim$Fp.mean
  model$post.scores.var    = sim$Fp.var
  model$post.scores        = sim$Fp 
  
  model$post.sigma2        = 1/sim$sigma2inv
  model$post.loadings.prob = 1.0-sim$pnz
  
 if (!is.null(more_args$save_max) && more_args$save_max>0) {
    model$post.cutpoints = sim$maxp
 }
 
  if (coda!="none") {
    if (coda=="all") {
      model$loadings.mcmc = get_coda(model, loadings=TRUE, scores=FALSE, scale=coda.scale)
      model$scores.mcmc   = get_coda(model, loadings=FALSE, scores=TRUE)
    } else if (coda=="loadings") {
      model$loadings.mcmc = get_coda(model, loadings=TRUE, scores=FALSE, scale=coda.scale)
    } else if (coda=="scores") {
      model$scores.mcmc   = get_coda(model, loadings=FALSE, scores=TRUE)
    }
  }

  return(model)
}


# Imputation on the original scale of the data
#
# This function imputes data on the original scale given a current value of the latent
# continuous data, using the empirical cdf.
#
# @param D Observed data
# @param Z Latent data (standardized)
# @return A vector of length \code{length(D[is.na(D)])} containing the
# imputed values
# @export


#impute <- function(D, Z) {
#  out = NA
#  for (i in 1:dim(D)[1]) {
#    im = quantile(D[i,], pnorm(Z[i,is.na(Z[i,]), type=1, na.rm=TRUE)
#    out = c(out, im)
#  }
#}

unload <- function() {
  detach("package:bfa", unload = TRUE)
  library.dynam.unload("bfa", system.file(package = "bfa"))
}



