#####################################################
## Source code for predict method of class 'hbglm' ##
#####################################################

#############################################################################
# Args:
#   object  - a fitted class 'hbglm' object
#   newdata - a list or data.frame (same format as in hbglm())
#   grpID.col - name of column with group ID info
#               (if NULL then value from object is used)
#   type      - What kind of prediction to do? `link' is default 
#               and returns the linear predictors, `response' returns
#               predicted probabilities. `ppp' refers to the posterior
#               predictive probability. `ppp.mean' returns a distribution 
#               of the expected response (over MCMC samples). `ppp.dist'
#               returns 'times' samples per mcmc sample from the response ppp
#   times     - Used only when type = 'ppp.sample' (see its docs)
#############################################################################
predict.hbglm <- function(object, newdata = NULL, grpID.col = NULL,
                   type = c("link", "response", "ppp.mean", "ppp.dist"), 
                   nburnin = 0, times = 0, print.level = 0, ...)
{
    nburnin <- if (nburnin <= 0) 0
    if (nburnin >= object$nsamples) stop("nburnin in predict >= num samples.")
  
    samples <- object$samples # grab drawn samples from 'hbglm' object
    family <- object$family

    # Fix the type of prediction requested
    if (!any(type %in% c("link", "response", "ppp.mean", "ppp.dist")))
        stop("Invalid arg 'type' to predict.hbglm()")
    if (length(type) > 1) type <- "link"   # default option

    # Grab the old formulas 
    parsed.fm <- parseFormula(object$formula)
    if (!is.null(grpID.col)) parsed.fm$grpID.col <- grpID.col
    parsed.fixed.fm <- if (is.null(object$formula.fixed)) NULL else {
                         formula.fixed <- object$formula.fixed
                         tm <- terms(formula.fixed)
                         list(formula = formula.fixed,
                              fixed.cov = attr(tm, "term.labels"),
                              intercept = attr(tm, "intercept"),
                              response = ifelse(attr(tm, "response"),
                                all.vars(formula.fixed)[1] , NULL))
                        }
    # Make model matrices with the new data
    model <- if (is.null(newdata)) object$model else 
        model.hbglm(newdata, parsed.fm, parsed.fixed.fm, predict=TRUE)

    # Compute linear predictor, eta
    grp.indx <- lapply(1:model$J, 
                       function(j) which(model$grp.indx == model$grp.labels[j]))
    compute.linpred <- function(beta, alpha = NULL)
    {
        eta <- if(model$has.fixed) model$mat.fixed %*% alpha else 
                                   rep(0, model$N)
        for (j in 1:model$J) { # loop over groups 
          eta[grp.indx[[j]]] <- eta[grp.indx[[j]]] + 
                                model$mat.rand.split[[j]]$X %*% beta[j, ]
        }
        
        return(eta)
    }

    ##############################
    # Non-Bayesian prediction
    ##############################
    if (type == "link" || type == "response") {
      stats <- sample.stats(samples, nburn = nburnin) # Compute sample stats
      eta <- compute.linpred(matrix(stats$beta[ , 1], nrow = model$J),
                 alpha = if (model$has.fixed) stats$alpha[ , 1] else NULL)
      names(eta) <- c(1:model$N)
      if (type == "link") return(eta)
      if (type == "response") {
          if (is.null(object$family$linkinv)) {
              warning(cat("Missing inverse link func for response prediction.",
                          "Returning linear predictor instead."))
              return(eta)
          }
          return(object$family$linkinv(eta))
      }
    }  

    ######################################################################
    # Bayesian posterior predictive response calculation
    ######################################################################
    #  For GLM: expected value of response = linkinv(linear predictor)
    #  Computes distribution of the expected response.
    if (type == "ppp.mean") {
        if (is.null(object$family$linkinv))
            stop(cat("Missing inverse link func for response prediction."))
        # Computes mean & variance of response for one coeff sample
        mean.response <- function(ii) {
            eta <- compute.linpred(matrix(samples$beta[ , ii], nrow = model$J),
                alpha = if (model$has.fixed) samples$alpha[ , ii] else NULL)
            return(object$family$linkinv(eta))
        }
        response <- list(samp = sapply(c((1+nburnin):object$nsamples), 
                                          mean.response))
        rownames(response$samp) <- rownames(newdata)
        colnames(response$samp) <- c((1+nburnin):object$nsamples)
        response$stats <- get.stats(response$samp)
        rownames(response$stats) <- rownames(newdata)
        return(response)
    }
    
    #####################################################################   
    # Bayesian posterior predictive probability of response
    #####################################################################   
    if (type == "ppp.dist") {
        # Compute linear predictor for all coeff samples
        eta.mat <- sapply(c((1+nburnin):object$nsamples), 
          function(ii) {
            compute.linpred(matrix(samples$beta[ , ii], nrow = model$J),
            alpha = if (model$has.fixed) samples$alpha[ , ii] else NULL)
          }
        )
        # Evaluate log posterior predictive probability
        log.ppp <- function(resp.vec) {
          Ns <- ncol(eta.mat)
          logp <- if (family$has.tau) {
            sapply(1:Ns, function(ii) {
                eta.ii <- eta.mat[ , ii]
                sum(sapply(1:model$J, function(jj)
                    family$loglik(eta.ii[grp.indx[[jj]]], resp.vec[[jj]], 
                                  var = samples$tau[jj, ii])))
            }) 
          } else { 
            sapply(1:Ns, function(ii) family$loglik(eta.mat[ , ii], resp.vec))
          }
          return(log(sum(exp(logp)) / Ns))
        }
        # Sample from the posterior predictive probability
        times <- if(times < 1) ncol(eta.mat) else times # num samples to draw
        report.freq <- 5
        if (print.level) { cat(paste0("\nDrawing ", times, " samples from ",
                                      "posterior predictive distribution\n"))}
        ppp.samp <- matrix(rep(NA, times * nrow(eta.mat)), ncol = times)
        # Initial sample
        ppp.samp[ , 1] <- family$linkinv(rowMeans(eta.mat))
        samp <- ppp.samp[ , 1]
        # Draw samples
        for (ii in 2:times) {
          samp <- MfU.Sample(samp, f = log.ppp, uni.sampler = "slice",
              control = MfU.Control(n = length(samp), slice.m = 10)) 
              #slice.lower =
              #bounds$alpha.lo, slice.upper=bounds$alpha.hi))
          if (print.level && (ii %% report.freq == 0))
              cat(paste0("\tDrawn ", ii, " / ", times, " ppp samples.\n"))
          ppp.samp[ , ii] <- samp
        }
        response <- list(samp = ppp.samp)
        rownames(response$samp) <- rownames(newdata)
        colnames(response$samp) <- c((1+nburnin):object$nsamples)
        response$stats <- get.stats(response$samp)
        rownames(response$stats) <- rownames(newdata)
        return(response)
    }
    
    stop("Unrecognized 'type' arg in predict.hbglm().") 
} 
