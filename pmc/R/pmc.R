#' pmc
#'
#' Performs a phylogenetic monte carlo between modelA and modelB
#'
#' Simulates data under each model and returns the distribution of 
#' likelihood ratio, L(B)/L(A), under for both simulated datasets.
#' @param tree A phylogenetic tree.  Can be phylo (ape) or ouch tree
#' @param data The data matrix
#' @param modelA a model from the list, or a custom model, see details
#' @param modelB any other model from the list, or custom model, see details
#' @param nboot number of bootstrap replicates to use
#' @param optionsA additional arguments to modelA
#' @param optionsB additional arguments to modelB
#' @param ... additional arguments to both fitting methods
#' @param mc.cores number of parallel cores to use
#' @return 
#' list with the nboot likelihood ratios obtained from fitting both models
#' to data simulated by model A, and the nboot likelihood ratios obtained
#' by fitting both models to simulations from model B, and the likelihood 
#' ratio between the original MLE estimated models from the data.  
#' @import parallel
#' @importFrom dplyr bind_rows bind_cols
#' @importFrom tidyr gather_
#' @export
#' @examples
#' library("geiger")
#' geo=get(data(geospiza))
#' tmp=treedata(geo$phy, geo$dat)
#' phy=tmp$phy
#' dat=tmp$data[,1]
#' \donttest{ 
#' pmc(phy, dat, "BM", "lambda", nboot = 20, mc.cores=1)
#' }
pmc <- function(tree, data, modelA, modelB, nboot = 500, optionsA = list(), optionsB = list(), ...,  mc.cores = parallel::detectCores()){
  fit_A <- do.call(pmc_fit, c(list(tree = tree, data = data, model = modelA), c(optionsA, list(...))))
  fit_B <- do.call(pmc_fit, c(list(tree = tree, data = data, model = modelB),  c(optionsB, list(...))))

  lr_orig <- -2 * (logLik(fit_A) - logLik(fit_B))

## 1000 Simulations under each model
  A_sims <- format_sims(simulate(fit_A, nboot))
  B_sims <- format_sims(simulate(fit_B, nboot))

## here are the four fits
  AA <- mclapply(1:nboot, function(i) update(fit_A, A_sims[,i]), mc.cores = mc.cores)
  AB <- mclapply(1:nboot, function(i) update(fit_B, A_sims[,i]), mc.cores = mc.cores)
  BA <- mclapply(1:nboot, function(i) update(fit_A, B_sims[,i]), mc.cores = mc.cores)
  BB <- mclapply(1:nboot, function(i) update(fit_B, B_sims[,i]), mc.cores = mc.cores)

## which create 2 distributions
  null_dist = -2 * (sapply(AA, logLik) - sapply(AB, logLik))
  test_dist = -2 * (sapply(BA, logLik) - sapply(BB, logLik))

  suppressWarnings({
  par_dists <- dplyr::bind_rows(tidy_pars(AA), tidy_pars(AB), tidy_pars(BA), tidy_pars(BB)) 
  })     

  out <- list(lr = lr_orig, null = null_dist, test = test_dist, par_dists = par_dists, A = fit_A, B = fit_B) 
  class(out) <- "pmc"
  out
}



## Helper functions because no one returns tidy models.
format_sims <- function(s){
  if(is.list(s))
    s <- dplyr::bind_cols(s)
  s
}
tidy_pars <- function(model, label = deparse(substitute(model))){
  mtrx <- sapply(model, function(x) {
                 out <- coef(x)
                 if(is.list(out))
                   out <- unlist(out)
                 out
          })
  tmp <- data.frame(t(rbind(mtrx, rep = 1:dim(mtrx)[[2]])))
  who <- names(tmp)[which(names(tmp)!="rep")]
  data.frame(comparison = label, tidyr::gather_(tmp, "parameter", "value", who))
}




#' Fit any model used in PMC 
#'
#' The fitting function used by pmc to generalize fitting to both geiger and ouch models.
#' @param tree a phylogenetic tree. can be ouch or ape format
#' @param data trait data in ape or ouch format
#' @param model the name of the model to fit, 
#' @param ... whatever additional options would be provided 
#' to the model fit
#' @return the object returned by the model fitting routine (gfit for geiger, hansen/brown for ouch) 
#' @import geiger ouch 
#' @export
pmc_fit <- function(tree, data, model, ...){
  # Figure out if we need ape/geiger based formats or ouch formats
  fitContinuous_types <- c("BM", "OU", "lambda", "kappa", 
                           "delta", "EB", "white", "trend")

  if(model %in% fitContinuous_types){
    type <- "fitContinuous"  
  } else if(model %in% c("brown", "hansen")){
    type <- "hansen"
  } else {
    stop(paste("Model", model, "not recognized"))
  }

  ## Run a fitContinuous fit ##
  if(type == "fitContinuous"){
    object <- fitContinuous(phy = tree, dat = data, model = model, ..., ncores=1)

  
  } else if(type == "hansen"){
    if(model == "hansen")
       object <- hansen(data = data, tree = tree, ...)
    else if(model == "brown")
       object <- brown(data = data, tree = tree, ...)

  } else {
    stop("Error: format not recognized.")
  }

  object
}


#' @import ggplot2
plot.pmc <- function(x, ...){
    df <- data.frame(null = x$null, test = x$test)
    dat <- tidyr::gather_(df, "variable", "value", gather_cols = c("null", "test"))
    ggplot(dat) + geom_density(aes_string("value", fill = "variable"), alpha = .7) +
           geom_vline(xintercept = x$lr, lwd = 1, lty = 2)
}



