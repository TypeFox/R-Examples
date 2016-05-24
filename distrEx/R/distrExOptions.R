.distrExOptions <- list(
    GLIntegrateTruncQuantile = 10*.Machine$double.eps,
    GLIntegrateOrder = 500,
    MCIterations = 1e5,
    ElowerTruncQuantile = 1e-7,
    EupperTruncQuantile = 1e-7,
    ErelativeTolerance = .Machine$double.eps^0.25,
    m1dfLowerTruncQuantile = 0,
    m1dfRelativeTolerance = .Machine$double.eps^0.25,
    m2dfLowerTruncQuantile = 0,
    m2dfRelativeTolerance = .Machine$double.eps^0.25,
    nDiscretize = 100,
    hSmooth = 0.05,
    IQR.fac = 15
)

distrExOptions <- function(...) {
    if (nargs() == 0) return(.distrExOptions)
    current <- .distrExOptions
    temp <- list(...)
    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
            list = temp <- arg,
            character = return(.distrExOptions[arg]),
            stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")
    changed <- current[n]
    current[n] <- temp
    if (sys.parent() == 0) 
        env <- asNamespace("distrEx") 
    else 
        env <- parent.frame()
    assign(".distrExOptions", current, envir = env)

    invisible(current)
}

getdistrExOption <- function(x) distrExOptions(x)[[1]]
distrExoptions <- distrExOptions 
