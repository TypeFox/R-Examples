proLik <- function(full.model, component,
	G = TRUE, negative = FALSE,
	nsample.units = 3, nse = 3,
	alpha = 0.05, tolerance = 0.001,
	parallel = FALSE, ncores = getOption("mc.cores", 2L)){
  s2 <- full.model$sigma2
  gamma.ind <- which(names(full.model$gammas) == component)
  bound.check <- names(full.model$gammas.con)[gamma.ind]
    warned <- FALSE
    if(bound.check == "Boundary"){
      warned <- TRUE
      warning("Boundary parameter: confidence interval estimation may produce strange behavior - proceed with caution)")
      }
  gamma.est <- full.model$gammas[gamma.ind][[1]]
  if(full.model$gammas[gamma.ind] == summary(full.model)$varcomp[gamma.ind, 2]) s2 <- 1
  std.err <- sqrt(diag(aiFun(full.model))[gamma.ind])
  if(is.na(std.err) | std.err == 0) std.err <- 0.1 
  full.mod2 <- asreml::update.asreml(object = full.model, start.values = TRUE)$gammas.table
  chi.val <- 0.5 * qchisq(alpha, df = 1, lower.tail = FALSE)

  proLik_keep_uniQUe_UCL <- list(gam = NULL, lambdas = NULL)
  tmpLRTU <- function(st, chi){
      proLik_keep_uniQUe_UCL[[1]] <<- c(proLik_keep_uniQUe_UCL[[1]], st)
      lrt <- constrainFun(st, full.model, full.mod2, component, G)
      proLik_keep_uniQUe_UCL[[2]] <<- c(proLik_keep_uniQUe_UCL[[2]], lrt)
      abs(chi - lrt)
  }
  proLik_keep_uniQUe_LCL <- list(gam = NULL, lambdas = NULL)
  tmpLRTL <- function(st, chi){
      proLik_keep_uniQUe_LCL[[1]] <<- c(proLik_keep_uniQUe_LCL[[1]], st)
      lrt <- constrainFun(st, full.model, full.mod2, component, G)
      proLik_keep_uniQUe_LCL[[2]] <<- c(proLik_keep_uniQUe_LCL[[2]], lrt)
      abs(chi - lrt)
  }


  Uint <- c(gamma.est, gamma.est + (nse * std.err))
  Lint <- c(gamma.est - (nse*std.err), gamma.est)
  if(!negative & Lint[1] < 0) Lint[1] <- Lint[2] / 15
  if(negative == TRUE & Uint[2] > 1) Uint[2] <- 0.999999
  if(negative == TRUE & Lint[1] < -1) Lint[1] <- -0.999999
  if(parallel){
     tmpUCL <- parallel::mcparallel(expr = expression(c(optimize(f=tmpLRTU, interval = Uint, chi = chi.val, tol = tolerance), proLik_keep_uniQUe_UCL)))

  } else{
       UCL <- optimize(f = tmpLRTU, interval = Uint, chi = chi.val, tol = tolerance)
    }
  LCL <- optimize(f = tmpLRTL, interval = Lint, chi = chi.val, tol = tolerance)
  if(parallel){
     tmpUCL.out <- parallel::mccollect(tmpUCL, wait = TRUE)[[1]]
     UCL <- list(minimum = tmpUCL.out$minimum[[1]], objective = tmpUCL.out$objective[[1]])
     proLik_keep_uniQUe_UCL <- list(gam = tmpUCL.out$gam, lambdas = tmpUCL.out$lambdas)
  }



  gamma.vec <- c(gamma.est, seq(gamma.est - (gamma.est - LCL$minimum)/2, gamma.est + (UCL$minimum - gamma.est)/2, length.out = nsample.units))

  if(parallel){
    if(length(gamma.vec) < ncores) ncores <- length(gamma.vec)
    profile <- list(lambdas = parallel::pvec(v = seq(1,length(gamma.vec),1), FUN = parConstrainFun, parameters = gamma.vec, full = full.model, fm2 = full.mod2, comp = component, G = G, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE), var.estimates = gamma.vec)
    } else{
    profile <- list(lambdas = vapply(gamma.vec, FUN = constrainFun, FUN.VALUE = vector("numeric", length = 1), full = full.model, fm2 = full.mod2, comp = component, G = G), var.estimates = gamma.vec)
      }

  profile$var.estimates <- c(profile$var.estimates, proLik_keep_uniQUe_LCL$gam, proLik_keep_uniQUe_UCL$gam)
  profile$lambdas <- c(profile$lambdas, proLik_keep_uniQUe_LCL$lambdas, proLik_keep_uniQUe_UCL$lambdas)
  ord.index <- order(profile$var.estimates)

    if(!negative & LCL$minimum < 0.01 & !warned){
      warning("Boundary parameter: confidence interval estimation may produce strange behavior - proceed with caution)")
      }


return(structure(list(lambdas = profile$lambdas[ord.index], 
	var.estimates = profile$var.estimates[ord.index] * s2, 
	UCL = UCL$minimum * s2, 
	LCL = LCL$minimum * s2, 
	component = component), class = "proLik"))
}


is.proLik <- function(x) inherits(x, "proLik")



