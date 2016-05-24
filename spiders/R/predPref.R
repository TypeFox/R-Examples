##' @title parameter estimation
##'
##' @description Estimates parameters of predator preferences model and calculates LRT.
##' Eaten and caught dataframes are indexed with rows across time points
##' and columns of prey species.
##'
##' @return A list of class 'predPref' with the following elements:
##'
##' null: parameters as estimated under the specified null hypothesis.
##' 
##' alt: parameters as estimated under the specified alternative hypothesis.
##' 
##' loglikH0: the null hypothesis log-likelihood, with constants not accounted for.
##' 
##' loglikH1: the alternative hypothesis log-likelihood, with constants not accounted for.
##' 
##' J: a column vector of dimension T containing the number of predators in
##' each time period.
##'
##' I: a column vector of dimension T containing the number of traps in each time period.
##'
##' LRT: the likelihood ratio test statistics.
##'
##' hypotheses: a 2-tuple of the user specified hypotheses.
##'
##' data.name: a character string giving the names of the data.
##'
##' @param eaten a dataframes of eatings preferences; TxS
##' @param caught a dataframes of caught prey species; TxS
##' @param hypotheses a 2-tuple specifying the null and alternative hypotheses, respectively
##' @param alpha LRT level of significance
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
##'
##' @seealso \code{\link{simPref}} \code{\link{summary.predPref}}
##'
##' @examples
##' # set parameters
##' Predators <- Traps <- 100
##' PreySpecies <- 2
##' Times <- 5
##' g <- matrix(sqrt(2), nrow=Times, ncol=PreySpecies)     # gamma
##' l <- matrix(seq(0.4,1.8,length.out=5)*sqrt(2), nrow=Times, ncol=PreySpecies) # ct
##'
##' # fit model
##' \dontrun{
##' fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=FALSE)
##' predPref(fdata$eaten, fdata$caught, hypotheses=c('ct', 'cst'))
##' }
##' 
##' @export
predPref <- function(eaten, caught, hypotheses = c("c", "Ct"), alpha=0.05, em_maxiter=1000) {

    ## check hypotheses specification
    hypotheses <- checkHypotheses(hypotheses)

    ## data
    dname <- paste(deparse(substitute(eaten)), "and", deparse(substitute(caught)))
    xNames <- colnames(eaten)
    yNames <- colnames(caught)

    ## ensure column named time
    if ( !("time" %in% xNames) || !("time" %in% yNames)) {
        stop("Need column named time.")
    }

    ## traps left out for different numbers of days?
    if ( !("adj" %in% xNames) ) {
        eaten$adj <- 1
    }
    if ( !("adj" %in% yNames) ) {
        caught$adj <- 1
    }

    ## prey names
    extraVars <- c("time", "adj")
    preyNames <- setdiff(xNames, extraVars)

    ## data errors
    X <- eaten[,preyNames]
    Y <- caught[,preyNames]

    if ( any(X < 0) || any(Y < 0) ) {
        stop("Count data can not be less than 0.")
    }

    ## predators (J), traps (I)
    J <- getTimeCounts(eaten, "adj")[,2]
    I <- getTimeCounts(caught, "adj")[,2] # total days traps were out each t
    
    ## data for calculations
    Xdst <- getTimeCounts(eaten, preyNames)[,preyNames, drop=F]
    Ydst <- getTimeCounts(caught, preyNames)[,preyNames, drop=F]

    ## do we run EM?
    EM <- ifelse(!any(X > 1), TRUE, FALSE)
    
    ## are data balanced
    BAL <- length(unique(J)) == 1 && length(unique(I)) == 1

    ## errors with time points
    if ( nrow(Xdst) != nrow(Ydst) ) {
        stop("Differing number of time points in eaten/caught data.")
    }

    ## check hypotheses
    if ( is.null(hypotheses) ) {
        stop("Need to speficy hypotheses.")
    }

    calcs <- calcHypotheses(hyp = hypotheses,
                      Xdst = Xdst, Ydst = Ydst, J=J, I=I,
                      balanced = BAL, EM=EM, em_maxiter = em_maxiter)
    llH0 <- calcs$llH0; llH1 <- calcs$llH1
    null <- calcs$null; alt <- calcs$alt
    df <- calcs$df

    ## LRT stat
    Lambda <- -2*(llH0 - llH1)
    lrt <- list(Lambda = Lambda, df = df,
                p.value = pchisq(Lambda, df=df, lower.tail=F))
    
    out <- list(alt=alt, null=null,
                loglikH1=llH1, loglikH0=llH0,
                numPredators=J, numTraps=I,
                LRT = lrt, hypotheses = hypotheses, data.name = dname)
    class(out) <- "predPref"
    out
}
