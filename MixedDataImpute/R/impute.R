#' Generate multiply imputed datasets 
#' 
#' @description { Imputations are generated using nonparametric Bayesian 
#' joint models (specifically the hierarchcially coupled mixture model with local dependence 
#' described in Murray and Reiter (2015); see citation(MixedDataImpute)
#'  or http://arxiv.org/abs/1410.0438). 
#'}
#' @param X A data frame of categorical variables (as factors)
#' @param Y A matrix or data frame of continuous variables
#' @param kz Number of top-level clusters
#' @param kx Number of X-model clusters
#' @param ky Number of Y-model clusters
#' @param hyperpar A list of hyperparameter values (see \code{hcmm_hyperpar})
#' @param num.impute Number of imputations
#' @param num.burnin Number of MCMC burn-in iterations
#' @param num.skip Number of MCMC iterations between saved imputations
#' @param thin.trace If negative, only save the num.impute datasets. If positive,
#' save summaries of the model state at every \code{thin.trace} iterations for
#' diagnostic purposes.
#' @param status Interval at which to print status messages
#'
#' @return A list with three elements: 
#' 
#' \code{imputations} A list of length \code{num.impute}. Each element is an imputed dataset.
#' 
#' \code{trace} MCMC output (currently the component sizes for the three mixture indices)
#' 
#' \code{model} An interface to the C++ object containing the current state
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' library(MixedDataImpute)
#' library(mice) # For the functions implementing combining rules
#' 
#' data(sipp08) 
#' 
#' set.seed(1)
#' n = 1000
#' s = sample(1:nrow(sipp08), n)
#' 
#' Y = sipp08[s,1:2]
#' Y[,1] = log(Y[,1]+1)
#' X = sipp08[s,-c(1:2,9)] # Also removes occ code, which has ~23 levels
#' 
#' # MCAR with probability 0.2, for illustration purposes (not matching the paper)
#' 
#' Y[runif(n)<0.2,1] = NA
#' Y[runif(n)<0.2,2] = NA
#' for(j in 1:ncol(X)) X[runif(n)<0.2,j] = NA
#' 
#' kz = 15
#' ky = 60
#' kx = 90
#' 
#' num.impute = 5
#' num.burnin = 10000
#' num.skip = 1000
#' thin.trace = 10
#' 
#' imp = hcmm_impute(X, Y, kz=kz, kx=kx, ky=ky, 
#'                   num.impute=num.impute, num.burnin=num.burnin, 
#'                   num.skip=num.skip, thin.trace=thin.trace)
#' 
#' # Example of getting MI estimates for a regression, using the
#' # pooling functions in mice
#' form = total_earnings~age+I(age^2) + sex*I(own_kid!=0)
#' 
#' fits = lapply(imp$imputations, function(dat) lm(form, data=dat))
#' pooled_ests = pool(as.mira(fits))
#' summary(pooled_ests)
#' 
#' # original, complete data estimates for comparison
#' comdat = sipp08[s,]
#' comdat[,1] = log(comdat[,1]+10)
#' summary(lm(form, data=comdat))
#' 
#' # true population values for comparison
#' pop = sipp08
#' pop[,1] = log(pop[,1]+10)
#' summary(lm(form, data=pop))
#' 
#'}
#' 
hcmm_impute = function(X, Y, kz, kx, ky, hyperpar=NULL, num.impute, num.burnin, num.skip, thin.trace=-1, status=50) {
  
  Y = as.matrix(Y)
  hcmmdat = prepare_data(X, Y)
  if(is.null(hyperpar)) hyperpar = hcmm_hyperpar(hcmmdat)
  
  Xmask = matrix(as.integer(is.na(X)), nrow=nrow(X))
  Ymask = matrix(as.integer(is.na(Y)), nrow=nrow(Y))
  
  t0 = Sys.time()
  model = hcmmld$new(kz, kx, ky)
  model$set_data(as.matrix(hcmmdat$Yscaled), Ymask, 
                 as.matrix(hcmmdat$Xint), Xmask, hcmmdat$cx, 
                 as.matrix(.make_lookup(hcmmdat$cx)-1))
  model$set_hyperpar(hyperpar)
  
  cat("\n%%% Burn-in", paste0('(',date(),')')," %%%\n\n")
  xx = model$mcmc(num.burnin, status, 0, -1)
  
  imputations = vector(mode="list", length = num.impute)
  traces = vector(mode="list", length = num.impute)
  cat("\n%%% MCMC %%%\n")
  for(i in 1:num.impute) {
    cat("\n%%% Imputation",i, paste0('(',date(),')'), "%%%\n\n")
    traces[[i]] = model$mcmc(num.skip, status, num.burnin+(i-1)*num.skip, thin.trace)
    tmp = remap_imputations(model$X, model$Y, hcmmdat)
    imputations[[i]] = data.frame(tmp$Y, tmp$X)
  }
  t1 = Sys.time()
  cat("\nComplete, total time",round(as.numeric(difftime(t1, t0), units="mins"), 2), "minutes\n")
  
  Z_alloc = do.call(rbind, lapply(traces, function(x) x$Z_alloc))
  X_alloc = do.call(rbind, lapply(traces, function(x) x$X_alloc))
  Y_alloc = do.call(rbind, lapply(traces, function(x) x$Y_alloc))
  tr = list(Z_alloc=Z_alloc, X_alloc=X_alloc, Y_alloc=Y_alloc)
  
  return(list(imputations=imputations, trace=tr, model=model))
  
}