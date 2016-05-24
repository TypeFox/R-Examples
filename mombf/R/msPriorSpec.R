###
### msPriorSpec.R
###


##-----------------------------------------------------------------------------
## Invoked by validObject() method.
valid_msPriorSpec <- function(object) {

    msg <- NULL

    ## Ensure first two slots are scalar
    if (length(object@priorType) != 1) {
        msg <- c(msg,
                 sprintf("slot %s must be of length 1",
                         sQuote("priorType")))
        object@priorType <- object@priorType[1]         # process first only
    }
    if (length(object@priorDistr) != 1) {
        msg <- c(msg,
                 sprintf("slot %s must be of length 1",
                         sQuote("priorDistr")))
        object@priorDistr <- object@priorDistr[1]       # process first only
    }

    ## Ensure valid type
    valid_prior_types <- c("coefficients",
                           "modelIndicator",
                           "nuisancePars")
    found <- object@priorType %in% valid_prior_types
    if (!found) {
        msg <- c(msg,
                 sprintf("slot %s must be one of: %s",
                         sQuote("priorType"),
                         paste(dQuote(valid_prior_types), collapse=", ")))
    }

    ##-------------------------------------------------------------------------
    has_either_or <- function(this, that, among) {
        this_in <- this %in% among
        that_in <- that %in% among

        (all(this_in) & !any(that_in)) || (!any(this_in) & all(that_in))
    }


    ## Validate rest based on type
    switch(object@priorType,
           coefficients = {
               valid_coef_prior_distrs <- c("pMOM",
                                            "piMOM",
                                            "peMOM",
                                            "zellner")
               found <- object@priorDistr %in% valid_coef_prior_distrs
               if (!found) {
                   msg <- c(msg,
                            sprintf("slot %s must be one of: %s",
                                    sQuote("priorDistr"),
                                    paste(dQuote(valid_coef_prior_distrs),
                                          collapse=", ")))
               }
               if (object@priorDistr == "pMOM") {
                   reqd_pmom_prior_par <- "r"
                   found <- reqd_pmom_prior_par %in% names(object@priorPars)
                   if (!found) {
                       msg <- c(msg,
                                sprintf("slot %s must be vector with named element %s (when using %s distr)",
                                        sQuote("priorPars"),
                                        dQuote(reqd_pmom_prior_par),
                                        sQuote("pMOM")))
                   }
               }
               reqd_prior_pars <- c("tau",
                                    "a.tau",
                                    "b.tau")
               found <- has_either_or(reqd_prior_pars[1],
                                      reqd_prior_pars[-1],
                                      names(object@priorPars))
               if (!found) {
                   msg <- c(msg,
                            sprintf("slot %s must be vector with named element(s) %s, or %s",
                                    sQuote("priorPars"),
                                    dQuote(reqd_prior_pars[1]),
                                    paste(dQuote(reqd_prior_pars[-1]),
                                          collapse=" and ")))
               }
           },
           modelIndicator = {
               valid_ind_prior_distrs <- c("uniform",
                                           "binomial")
               found <- object@priorDistr %in% valid_ind_prior_distrs
               if (!found) {
                   msg <- c(msg,
                            sprintf("slot %s must be one of: %s",
                                    sQuote("priorDistr"),
                                    paste(dQuote(valid_ind_prior_distrs), collapse=", ")))
               }
               if (object@priorDistr == "binomial") {
                   reqd_prior_pars <- c("p",
                                        "alpha.p",
                                        "beta.p")
                   found <- has_either_or(reqd_prior_pars[1],
                                          reqd_prior_pars[-1], 
                                          names(object@priorPars))
                   if (!found) {
                       msg <- c(msg,
                                sprintf("slot %s must be vector with named element(s) %s, or %s (when using %s distr)",
                                        sQuote("priorPars"),
                                        dQuote(reqd_prior_pars[1]),
                                        paste(dQuote(reqd_prior_pars[-1]),
                                              collapse=" and "),
                                        sQuote("binomial")))
                   }
               }
           },
           nuisancePars = {
               valid_nuisance_prior_distrs <- "invgamma"
               found <- object@priorDistr %in% valid_nuisance_prior_distrs
               if (!found) {
                 msg <- c(msg,
                          sprintf("slot %s must be %s",
                                  sQuote("priorDistr"), 
                                  dQuote(valid_nuisance_prior_distrs)))
               }
               reqd_prior_pars <- c("alpha",
                                    "lambda")
               found <- reqd_prior_pars %in% names(object@priorPars)
               if (!all(found)) {
                   msg <- c(msg,
                            sprintf("slot %s must contain named elements %s",
                                    sQuote("priorPars"),
                                    paste(dQuote(reqd_prior_pars), 
                                          collapse=", ")))
               }
           })

    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
}

setValidity("msPriorSpec", valid_msPriorSpec)


##-----------------------------------------------------------------------------
is.msPriorSpec <- function(x) {
    extends(class(x), "msPriorSpec")
}


##-----------------------------------------------------------------------------
## Generators for this object.
msPriorSpec <- function(priorType=c("coefficients",
                                    "modelIndicator",
                                    "nuisancePars"),
                        priorDistr,
                        priorPars) {
    ## Check arguments
    stopifnot(is.character(priorType)  && length(priorType) == 1)
    stopifnot(is.character(priorDistr) && length(priorDistr) == 1)
    stopifnot(is.vector(priorPars))

    ## Create new class
    new("msPriorSpec",
        priorType=match.arg(priorType),
        priorDistr=priorDistr,
        priorPars=priorPars)
}


igprior <- function(alpha=.01, lambda=.01) {
    new("msPriorSpec",priorType='nuisancePars',priorDistr='invgamma',priorPars=c(alpha=alpha,lambda=lambda))
}

momprior <- function(tau, tau.adj=10^6, r=1) {
    new("msPriorSpec", priorType="coefficients", priorDistr="pMOM", priorPars=c(tau=tau,tau.adj=tau.adj,r=r))
}

imomprior <- function(tau, tau.adj=10^6) {
    new("msPriorSpec", priorType="coefficients", priorDistr="piMOM", priorPars=c(tau=tau, tau.adj=tau.adj))
}

emomprior <- function(tau, tau.adj=10^6) {
    new("msPriorSpec", priorType="coefficients", priorDistr="peMOM", priorPars=c(tau=tau, tau.adj=tau.adj))
}

zellnerprior <- function(tau, tau.adj=10^6) {
    new("msPriorSpec", priorType="coefficients", priorDistr="zellner", priorPars=c(tau=tau, tau.adj=tau.adj))
}

modelunifprior <- function() {
    new("msPriorSpec",priorType='modelIndicator',priorDistr='uniform',priorPars=double(0))
}

modelbinomprior <- function(p=0.5) {
    new("msPriorSpec",priorType='modelIndicator',priorDistr='binomial',priorPars=c(p=p))
}

modelbbprior <- function(alpha.p=1, beta.p=1) {
    new("msPriorSpec",priorType='modelIndicator',priorDistr='binomial',priorPars=c(alpha.p=alpha.p,beta.p=beta.p))
}


igprior <- function(alpha=.01, lambda=.01) {
    new("msPriorSpec",priorType='nuisancePars',priorDistr='invgamma',priorPars=c(alpha=alpha,lambda=lambda))
}
