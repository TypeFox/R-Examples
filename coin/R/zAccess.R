### generic method for extracting p-values from objects
setGeneric("pvalue",
    function(object, ...) {
        standardGeneric("pvalue")
    }
)

setMethod("pvalue",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        object@pvalue(q)
    }
)

setMethod("pvalue",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        pvalue(object@distribution, object@statistic@teststatistic)
    }
)

setMethod("pvalue",
    signature = "MaxTypeIndependenceTest",
    definition = function(object,
        method = c("global", "single-step", "step-down", "unadjusted"),
###     combinations = c("free", "restricted"), # placeholder
        distribution = c("joint", "marginal"),
        type = c("Bonferroni", "Sidak"), ...) {
            method <- match.arg(method,
                          choices = c("global", "single-step", "step-down",
                                      "unadjusted", "discrete"),
                          several.ok = TRUE)[1]
            if (method == "discrete")
                warning(sQuote(paste("method =", dQuote(method))),
                        " is deprecated; see ", sQuote("?pvalue"))
            distribution <- match.arg(distribution)
            type <- match.arg(type)

            C <- attr(object@statistic@xtrans, "contrast")
            if (!is.null(C) && method != "global")
                warning(paste("p-values may be incorrect due to violation of",
                              "the subset pivotality condition"))

            if (method == "global")
                pvalue(object@distribution, object@statistic@teststatistic)
            else if (method == "single-step") {
                if (distribution == "joint")
                    singlestep(object, ...)
                else {
                    if (type == "Bonferroni")
                        marginal(object, bonferroni = TRUE,
                                 stepdown = FALSE, ...)
                    else
                        marginal(object, bonferroni = FALSE,
                                 stepdown = FALSE, ...)
                }
            } else if (method == "step-down") {
                if (distribution == "joint")
                    stepdown(object, ...)
                else {
                    if (type == "Bonferroni")
                        marginal(object, bonferroni = TRUE,
                                 stepdown = TRUE, ...)
                    else
                        marginal(object, bonferroni = FALSE,
                                 stepdown = TRUE, ...)
                }
            }
            ## <DEPRECATED>
            else if (method == "discrete")
                dbonf(object, ...)
            ## </DEPRECATED>
            else
                unadjusted(object, ...)
        }
)


### generic method for extracting mid-p-values from objects
setGeneric("midpvalue",
    function(object, ...) {
        standardGeneric("midpvalue")
    }
)

setMethod("midpvalue",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        object@midpvalue(q)
    }
)

setMethod("midpvalue",
   signature = "IndependenceTest",
   definition = function(object, ...) {
       midpvalue(object@distribution, object@statistic@teststatistic)
   }
)


### generic method for extracting p-value intervals from objects
setGeneric("pvalue_interval",
    function(object, ...) {
        standardGeneric("pvalue_interval")
    }
)

setMethod("pvalue_interval",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        object@pvalueinterval(q)
    }
)

setMethod("pvalue_interval",
   signature = "IndependenceTest",
   definition = function(object, ...) {
       pvalue_interval(object@distribution, object@statistic@teststatistic)
   }
)


### generic method for the permutation distribution from objects
setGeneric("dperm",
    function(object, x, ...) {
        standardGeneric("dperm")
    }
)

setMethod("dperm",
    signature = "NullDistribution",
    definition = function(object, x, ...) {
        vapply(x, object@d, NA_real_)
    }
)

setMethod("dperm",
    signature = "AsymptNullDistribution",
    definition = function(object, x, ...) {
        object@d(x)
    }
)

setMethod("dperm",
    signature = "IndependenceTest",
    definition = function(object, x, ...) {
        dperm(object@distribution, x)
    }
)


### generic method for the permutation distribution from objects
setGeneric("pperm",
    function(object, q, ...) {
        standardGeneric("pperm")
    }
)

setMethod("pperm",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        vapply(q, object@p, NA_real_)
    }
)

setMethod("pperm",
    signature = "AsymptNullDistribution",
    definition = function(object, q, ...) {
        object@p(q)
    }
)

setMethod("pperm",
    signature = "IndependenceTest",
    definition = function(object, q, ...) {
        pperm(object@distribution, q)
    }
)


### generic method for the permutation distribution from objects
setGeneric("qperm",
    function(object, p, ...) {
        standardGeneric("qperm")
    }
)

setMethod("qperm",
    signature = "NullDistribution",
    definition = function(object, p, ...) {
        vapply(p, object@q, NA_real_)
    }
)

setMethod("qperm",
    signature = "AsymptNullDistribution",
    definition = function(object, p, ...) {
        object@q(p)
    }
)

setMethod("qperm",
    signature = "IndependenceTest",
    definition = function(object, p, ...) {
        qperm(object@distribution, p)
    }
)


### generic method for the permutation distribution from objects
setGeneric("rperm", function(object, n, ...)
    standardGeneric("rperm"))

setMethod("rperm",
          signature = "NullDistribution",
          definition = function(object, n, ...) {
              qperm(object, runif(n))
          }
)

setMethod("rperm",
          signature = "IndependenceTest",
          definition = function(object, n, ...) {
              qperm(object, runif(n))
          }
)


### generic method for the permutation distribution from objects
setGeneric("support",
    function(object, ...) {
        standardGeneric("support")
    }
)

setMethod("support",
    signature = "NullDistribution",
    definition = function(object, ...) {
        object@support(...)
    }
)

setMethod("support",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        support(object@distribution, ...)
    }
)


### generic method for extracting statistics from objects
setGeneric("statistic",
    function(object, ...) {
            standardGeneric("statistic")
    }
)

setMethod("statistic",
    signature = "IndependenceLinearStatistic",
    definition = function(object,
        type = c("test", "linear", "standardized"), ...) {
            nc <- ncol(object@ytrans)
            nr <- ncol(object@xtrans)
            type <- match.arg(type)
            dn <- statnames(object)$dimnames
            switch(type,
                "test"         = stop(sQuote(paste("type =", dQuote("test"))),
                                      " not defined for objects of class ",
                                      dQuote("IndependenceLinearStatistic")),
                "linear"       = matrix(object@linearstatistic,
                                        nrow = nr, ncol = nc, dimnames = dn),
                "standardized" = matrix(object@standardizedlinearstatistic,
                                        nrow = nr, ncol = nc, dimnames = dn)
            )
        }
)

setMethod("statistic",
    signature = "IndependenceTestStatistic",
    definition = function(object,
        type = c("test", "linear", "standardized"), ...) {
            nc <- ncol(object@ytrans)
            nr <- ncol(object@xtrans)
            type <- match.arg(type)
            dn <- statnames(object)$dimnames
            switch(type,
                "test"         = object@teststatistic,
                "linear"       = matrix(object@linearstatistic,
                                        nrow = nr, ncol = nc, dimnames = dn),
                "standardized" = matrix(object@standardizedlinearstatistic,
                                        nrow = nr, ncol = nc, dimnames = dn)
            )
        }
)

setMethod("statistic",
    signature = "IndependenceTest",
    definition = function(object,
        type = c("test", "linear", "standardized"), ...) {
            statistic(object@statistic, type)
    }
)


### generic method for extracting expectations from objects
setGeneric("expectation",
    function(object, ...) {
        standardGeneric("expectation")
    }
)

setMethod("expectation",
    signature = "IndependenceLinearStatistic",
    definition = function(object, ...) {
        object@expectation
    }
)

setMethod("expectation",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        expectation(object@statistic, ...)
    }
)


### generic method for extracting the covariance matrix from objects
setGeneric("covariance",
    function(object, ...) {
        standardGeneric("covariance")
    }
)

setMethod("covariance",
    signature = "CovarianceMatrix",
    definition = function(object, ...) {
        object@covariance
    }
)

setMethod("covariance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, ...) {
        if (!extends(class(object@covariance), "CovarianceMatrix"))
            covariance(new("IndependenceLinearStatistic",
                           object, varonly = FALSE))
        else
            covariance(object@covariance)
    }
)

setMethod("covariance",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        covariance(object@statistic, ...)
    }
)


### generic method for extracting the variances
setGeneric("variance",
    function(object, ...) {
        standardGeneric("variance")
    }
)

setMethod("variance",
    signature = "Variance",
    definition = function(object, ...) {
        object@variance
    }
)

setMethod("variance",
    signature = "CovarianceMatrix",
    definition = function(object, ...) {
        diag(object@covariance)
    }
)

setMethod("variance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, ...) {
        variance(object@covariance)
    }
)

setMethod("variance",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        variance(object@statistic, ...)
    }
)
