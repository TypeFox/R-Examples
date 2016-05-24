makeInteger <- function(x, label, min = 1L) {
    x <- as.integer(x)
    if (is.na(x) || x < min)
        stop(label, " must be an integer greater/equal to ", min)
    x
}

##

anyNA <- function(x)
    if (!is.null(x) && !is.function(x)) any(is.na(x)) else FALSE

##

checkList <- function(passedList, defaultList, label = "'algo'") {
    ## NAs in list
    if (any(sapply(passedList,anyNA)))
        stop("NAs are not allowed in list ", label)

    ## unnamed elements in list
    if ("" %in% names(passedList))
        warning(label, " contained unnamed elements",
                call. = FALSE)

    ## ununsed elements in list
    unusedOptions <- setdiff(names(passedList), names(defaultList))
    unusedOptions <- setdiff(unusedOptions, "")
    if (length(unusedOptions))
        warning("unknown names in ", label, ": ",
                paste(unusedOptions, collapse = ", "),
                call. = FALSE)
}

##

mRU <- function(m, n)
    array(runif(m*n), dim = c(m,n))
mRN <- function(m, n)
    array(rnorm(m*n), dim = c(m,n))

##

mcList <- function(mc.control) {
    mc.settings <- list(mc.preschedule = TRUE,
                        mc.set.seed = TRUE,
                        mc.silent = FALSE,
                        mc.cores = getOption("mc.cores", 2L),
                        mc.cleanup = TRUE,
                        mc.allow.recursive = TRUE)
    checkList(mc.control, mc.settings, "'mc.control'")
    mc.settings[names(mc.control)] <- mc.control
    mc.settings
}

##

repair1c <- function(x, up, lo) {
    xadjU <- x - up
    xadjU <- xadjU + abs(xadjU)
    xadjL <- lo - x
    xadjL <- xadjL + abs(xadjL)
    x - (xadjU - xadjL)/2
}

##

Chapters <- c("Introduction",
              "Numerical analysis in a nutshell",
              "Linear equations and Least Squares problems",
              "Finite difference methods",
              "Binomial trees",
              "Generating random numbers",
              "Modeling dependencies",
              "A gentle introduction to financial simulation",
              "Financial simulation at work: some case studies",
              "Optimization problems in finance",
              "Basic methods",
              "Heuristic methods in a nutshell",
              "Portfolio optimization",
              "Econometric models",
              "Calibrating option pricing models")

## dividends until expiry (used in vanillaOption.*)

due <- function(D, tauD, tau, q) {
    if (any(q != 0) && any(D != 0))
        stop("dividend rate and dividend amount supplied")
    if (any(D != 0) && length(D) != length(tauD))
        stop("number of dividends and number times-to-dividend differ")
    if (any(D != 0)) {
        D2keep <-  tauD <= tau & tauD > 0
        D <- D[D2keep]
        tauD <- tauD[D2keep]
    }
    list(tauD = tauD, D = D)
}
