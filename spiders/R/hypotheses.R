##' @title calculate hypotheses
##'
##' @description calculates hypotheses, given a user specifed null and alternative
##'
##' @param hyp a 2-tuple specifying the null and alternative hypotheses, respectively
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param balanced boolean specifying balanced data or not
##' @param EM boolean specifying if EM algorithm should be used
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
calcHypotheses <- function(hyp, Xdst, Ydst, J, I, balanced, EM, em_maxiter) {

    ## create map of estimators
    fns <- list()
    fns[["1"]] <- est1; fns[["c"]] <- estC
    fns[["Cs"]] <- estCs; fns[["Ct"]] <- estCt
    fns[["Cst"]] <- estCst; fns[["gen"]] <- estGen

    ## calculate null and alternative hypotheses
    null <- fns[[hyp[1]]](Xdst, Ydst, J, I, EM, em_maxiter, balanced)
    alt <- fns[[hyp[2]]](Xdst, Ydst, J, I, EM, em_maxiter, balanced)        

    ## calculate degrees of freedom
    S <- ncol(Xdst)
    T <- nrow(Xdst)
    ST <- S*T
    nullDF <- length(null[["c"]])
    altDF <- ifelse(is.null(alt[["c"]]), ST, length(alt[["c"]]))
    df <- altDF - nullDF
    
    list(llH0=null[["ll"]], llH1=alt[["ll"]],
         null=null, alt=alt, df=df)    
}

##' function to check user specified hypotheses
##'
##' @param hyp a 2-tuple specifying the null and alternative hypotheses, respectively
checkHypotheses <- function(hyp) {

    hyp <- tolower((as.character(hyp)))
    
    ## initialize output
    H <- rep(0, 2)
    
    ## H0
    if ( grepl("t", hyp[1]) ) {
        H[1] <- "Ct"
    } else if ( grepl("s", hyp[1]) ) {
        H[1] <- "Cs"
    } else if ( grepl("c", hyp[1]) ) {
        H[1] <- "c"
    } else if ( grepl("1", hyp[1]) ) {
        H[1] <- "1"
    } else {
        stop("Null hypothesis specified incorrectly; please check documentation.")
    }

    ## H1
    if ( grepl("g", hyp[2]) ) {
        H[2] <- "gen"
    } else if ( grepl("st", hyp[2]) ) {
        H[2] <- "Cst"
    } else if ( grepl("t", hyp[2]) ) {
        H[2] <- "Ct"
    } else if ( grepl("s", hyp[2]) ) {
        H[2] <- "Cs"
    } else if ( grepl("c", hyp[2]) ) {
        H[2] <- "c"
    } else {
        stop("Alt hypothesis specified incorrectly; please check documentation.")
    }

    ## both
    if ( H[1] == H[2] ) {
        stop("Null and Alternative hypotheses must be different.")
    }

    H
}
