###############################################################################
#                        Package: HBglm                                       #
#                                                                             #
# Hierarchical Bayesian Regression via MCMC estimation                        #
#                                                                             #
#               Scientific Computing Group, Sentrana Inc.                     #
###############################################################################

###############################################################################
#                 Main user function                                          #
###############################################################################
hbglm <- function(formula, data, formula.fixed = NULL,
                  family = "gaussian",
                  sampler.control = hbglm.sampler.control(),
                  model.control = hbglm.model.control(), 
                  ncores = 1, print.level = 0, ...)
{
    start.time <- proc.time()[3]
    call <- match.call()  # Store function call

    ##############################################################
    ## Pre-processing                                            #
    ##############################################################
    # Parse the random effects and upper level formula
    parsed.fm <- parseFormula(formula)
    grpID.col <- parsed.fm$grpID.col
    # Handle fixed effects
    parsed.fixed.fm <- if (is.null(formula.fixed)) NULL else {
                         tm <- terms(formula.fixed)
                         list(formula = formula.fixed,
                              fixed.cov = attr(tm, "term.labels"),
                              intercept = attr(tm, "intercept"),
                              response = ifelse(attr(tm, "response"),
                                all.vars(formula.fixed)[1] , NULL)
                         )    
                       }
    if (!is.null(formula.fixed)) {
        # Check fixed effects formula for compatibility
        if (any(parsed.fixed.fm$fixed.cov %in% parsed.fm$lower.cov))
            stop("Fixed effect covariates in random effects.")
        if (any(parsed.fixed.fm$fixed.cov %in% parsed.fm$upper.cov))
            stop("Fixed effect covariates in upper level.")
        fr <- parsed.fixed.fm$response
        if (!(is.null(fr) || fr == parsed.fm$response))
            stop("Response in fixed effects formula doesn't match.")
        # Include response if not provided in the fixed effects formula 
        if (is.null(fr))
            parsed.fixed.fm$response <- parsed.fm$response

       # If intercept speficied in both fixed and random effects
       # Keep intercept ONLY in fixed effects
       if (parsed.fixed.fm$intercept && parsed.fm$intercept[1]) {
           parsed.fm$intercept[1] <- FALSE
           warning("Intercept in fixed effects; excluded from random effects")
       }
    }

    # Make model matrices and infer various model parameters
    model <- model.hbglm(data, parsed.fm, parsed.fixed.fm,
                         scale.data=model.control$scale.data)
    #model <- model.hbglm(data, parsed.fm, parsed.fixed.fm, scale.data = FALSE)
    
    # Handle the 'family' argument
    fam.name <- family
    type <- ifelse(family == "gaussian", "linear", "glm")
    if      (family == "gaussian") family = hbglm.gaussian
    else if (family == "binomial") family = hbglm.binomial
    else if (family == "poisson")  family = hbglm.poisson
    else stop("Unsupported regression type requested for lower level")
    family$fam.name <- fam.name

    # Handle constraints and initialization
    bounds <- get.box.bounds(family, model, 
                             user.constraints = model.control$constraints)
    init.vals <- get.init.vals(family, model, bounds, model.control)

    # Get the number of iterations of samples to draw
    nsamples <- sampler.control$num.samples
    if (nsamples <= 0) nsamples <- default.samples(model, sampler.control)
    nitrs <- nsamples #total num of MCMC/Gibbs iterations to run

    t1 <- proc.time()[3]
    if (print.level > 0) {
      cat(paste0("\n\nModel has:\n\t", model$J, " groups\n\t", model$K,
                 " random effects covariates\n\t", model$N, " observations"))
      if (model$has.fixed)
          cat(paste0("\n\t", model$M, " fixed effects covariates"))
      if (model$has.upper.level) {
          cat(paste0("\nModel is pooled with:\n\t", model$L, " upper ",
                     "level covariates"))
      } else cat(paste0("\nModel is unpooled"))
      cat(paste0("\nTotal number of parameters in model = ", 
                 length(init.vals$vec), "\n"))
      cat(paste0("MCMC estimation:\n\t", nsamples, " sampling iterations."))
      if (print.level > 1)
          cat(paste0("\nPreprocessing for estimation took ",
              round(t1 - start.time, 3), " sec.\n"))
    }
 
    ##############################################################
    ## MCMC Estimation                                           #
    ##############################################################
    # samples is a list with keys: beta, alpha, tau, theta, Sigma
    # Each key has a matrix or NULL associated with it
    # See mcmc.R for format of the returned matrices
    if (type == "linear") {
        samples <- hblinear(model, family, bounds, init.vals, nitrs,
                            print.level = print.level)
    } else {
        samples <- mcmc(model, family, bounds, init.vals, nitrs,
                        print.level = print.level)
    }
  
    ##############################################################
    ## Post-processing                                           #
    ##############################################################
    # Make an 'hbglm' class object to return 
    fit <- structure(list(
             formula         = formula,
             formula.fixed   = formula.fixed,
             nsamples        = nsamples,
             model.control   = model.control,
             sampler.control = sampler.control,
             grpID.col       = grpID.col,
             call            = call,
             model           = model,
             family          = family,
             samples         = samples
                         ), class = "hbglm"
    )
    if (print.level) {
      cat(paste0("\nTotal time spent in hbglm = ",
                 round(proc.time()[3] - start.time, 3), " seconds.\n"))
    }
    return(fit)
} 
################################################################################

# Helper function: sets the default number of samples to take
# Equal to: samp.factor * number of parameters in model
default.samples <- function(model, sampler.control)
{
    nparam <- (model$J + model$L + 0.5 * model$K + 1.5) * model$K + model$M
    return(nparam * sampler.control$samp.factor)
}
################################################################################
