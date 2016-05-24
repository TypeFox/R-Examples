### Conditional Expectation and Covariance
setClass("ExpectCovar",
    representation = representation(
        expectation = "numeric",
        covariance  = "matrix",
        dimension   = "integer"
   )
)

### Expectation and Covariance of the influence function
### (+ sum of weights)
setClass("ExpectCovarInfluence",
    contains = "ExpectCovar",
    representation = representation(
        sumweights = "numeric"
    )
)

### Covariance matrix
setClass("CovarianceMatrix",
    representation = representation(
        covariance = "matrix"
    )
)

### Variance only
setClass("Variance",
    representation = representation(
        variance = "numeric"
    )
)

### Virtual class for covariance and variance
setClassUnion("VarCovar",
    members = c("CovarianceMatrix", "Variance")
)

### Class for raw data: a set of `x' variables and a set of `y' variables,
### possibly blocked and with weights
setClass("IndependenceProblem",
    representation = representation(
        x       = "data.frame",
        y       = "data.frame",
        block   = "factor",
        weights = "numeric"
    ),
    validity = function(object) {
        dims <- ((nrow(object@x) == nrow(object@y)) &&
                 (nrow(object@x) == length(object@block)))
        dims <- dims && (length(object@block) == length(object@weights))
        Wint <- max(abs(object@weights - floor(object@weights))) < eps()
        block <- all(table(object@block) > 1L)
        NAs <- all(complete.cases(object@x) & complete.cases(object@y))
        (dims && block) && (NAs && Wint)
    }
)

### Class for transformed data, the `x' variables are transformed
### to a (n x p) matrix `xtrans' and the `y' variables to `ytrans' (n x q).
### `scores' is a matrix of scores
setClass("IndependenceTestProblem",
    contains = "IndependenceProblem",
    representation = representation(
        xtrans     = "matrix",
        ytrans     = "matrix",
        xtrafo     = "function",
        ytrafo     = "function"
    ),
    validity = function(object)
        (storage.mode(object@xtrans) == "double" &&
         storage.mode(object@ytrans) == "double")
)

### Linear statistic, expectation and covariance according to
### Strasser & Weber (1999)
setClass("IndependenceLinearStatistic",
    contains = "IndependenceTestProblem",
    representation = representation(
        linearstatistic = "numeric",
        expectation     = "numeric",
        covariance      = "VarCovar"
    )
)

### Tests based on linear statistics
setClass("IndependenceTestStatistic",
    contains = c("IndependenceLinearStatistic", "VIRTUAL"),
    representation = representation(
        teststatistic               = "numeric",
        standardizedlinearstatistic = "numeric"
    )
)

### teststatistic = standardizedlinearstatistic
setClass("ScalarIndependenceTestStatistic",
    contains = "IndependenceTestStatistic",
    representation = representation(
        alternative = "character",
        paired      = "logical"
    ),
    validity = function(object)
        object@alternative %in% c("two.sided", "less", "greater")
)

### teststatistic = max(abs(standardizedlinearstatistic))
setClass("MaxTypeIndependenceTestStatistic",
    contains = "IndependenceTestStatistic",
    representation = representation(
        alternative = "character"
    ),
    validity = function(object)
        object@alternative %in% c("two.sided", "less", "greater")
)

### teststatistic = quadform(linearstatistic)
setClass("QuadTypeIndependenceTestStatistic",
    contains = "IndependenceTestStatistic",
    representation = representation(
        covarianceplus = "matrix",
        df             = "numeric",
        paired         = "logical"
    )
)

### p-values
setClass("PValue",
    representation = representation(
        pvalue         = "function",
        midpvalue      = "function",
        pvalueinterval = "function",
        p              = "function",
        name           = "character"
    ),
    prototype = list(
        pvalue         = function(q) NA,
        midpvalue      = function(q) NA,
        pvalueinterval = function(q) NA,
        p              = function(q) NA,
        name           = NA_character_
    )
)

### Null distribution
setClass("NullDistribution",
    contains = "PValue",
    representation = representation(
        q          = "function",
        d          = "function",
        support    = "function",
        parameters = "list"
    ),
    prototype = list(
        q          = function(p) NA,
        d          = function(x) NA,
        support    = function() NA,
        parameters = list()
    )
)

### There are essentially three types of null distributions:
setClass("AsymptNullDistribution",
    contains = "NullDistribution",
    representation = representation(
        seed = "integer"
    )
)

setClass("ApproxNullDistribution",
    contains = "NullDistribution",
    representation = representation(
        seed = "integer"
    )
)

setClass("ExactNullDistribution",
    contains = "NullDistribution"
)

### the "fitted" test including data and everything
setClass("IndependenceTest",
    representation = representation(
        distribution = "PValue", # was: "NullDistribution",
        statistic    = "IndependenceTestStatistic",
        estimates    = "list",
        method       = "character",
        call         = "call"
    ),
    prototype = list(method = "General Independence Test")
)

### the "fitted" test for scalar linear statistics
setClass("ScalarIndependenceTest",
    contains = "IndependenceTest",
    representation = representation(
        parameter = "character",
        nullvalue = "numeric"
    ),
    prototype = list(parameter = "mu"),
    validity = function(object)
        extends(class(object@statistic),
                "ScalarIndependenceTestStatistic")
)

### possibly with confidence intervals
setClass("ScalarIndependenceTestConfint",
    contains = "ScalarIndependenceTest",
    representation = representation(
        confint    = "function",
        conf.level = "numeric"
    )
)

### max type test statistics
setClass("MaxTypeIndependenceTest",
    contains = "IndependenceTest",
    validity = function(object)
        extends(class(object@statistic),
                "MaxTypeIndependenceTestStatistic")
)

### quad form test statistics
setClass("QuadTypeIndependenceTest",
    contains = "IndependenceTest",
    validity = function(object)
        extends(class(object@statistic),
                "QuadTypeIndependenceTestStatistic")
)

### SymmetryProblems
setClass("SymmetryProblem",
    contains = "IndependenceProblem",
    validity = function(object) {
        if (ncol(object@x) != 1L || !is.factor(object@x[[1L]]))
            stop(sQuote("x"), " slot does not contain a single factor")
        if (!is_completeblock(object))
            stop(sQuote("object"),
                 " is not a an unreplicated complete block design")
        TRUE
    }
)
