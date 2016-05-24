

#' Derive mutation rate estimation from outbreak's outputs
#'
#' The function \code{get.mu} is used to obtain a distribution of the mutation
#' rate from outbreaker's ouptput (functions \code{outbreaker} and
#' \code{outbreaker.parallel}). The mutation rates used in outbreaker's model
#' are expressed per generation of infection, which can be problematic to
#' interprete biologically. \code{get.mu} derives classical estimates of the
#' mutation rate per unit of time, with one value being estimated for each
#' chain of the MCMC. By default, the mutation rate is expressed in number of
#' nucleotide changes per unit time and per genome. If \code{genome.size} is
#' provided, the mutation rate is expressed in number of nucleotide changes per
#' unit time and per site.
#'
#'
#' @aliases get.mu
#'
#' @export
#'
#' @param x the output of \code{outbreaker} or \code{outbreaker.parallel}.
#' @param burnin an integer indicating the number of steps of the MCMC to be
#' discarded as burnin period. Defaults to 20,000.
#' @param genome.size the size of the genome; if not provided, mutation rate
#' will be expressed in number of mutations per unit of time and per genome.
#' @return A vector of mutation rates derived from the MCMC.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @examples
#'
#' ## load data
#' data(fakeOutbreak)
#' attach(fakeOutbreak)
#'
#' mu <- get.mu(res, genome.size=ncol(dat$dna))
#' hist(mu, col="grey",
#'      main="Inferred distribution of mu",
#'      xlab="mutations/site/day")
#' abline(v=1e-4,lty=2, lwd=4, col="royalblue")
#' mtext(side=3, "Dashed line = actual value")
#'
#' detach(fakeOutbreak)
#'
get.mu <- function(x, burnin=2e4, genome.size=NULL){

    ## GET DATA ##
    ## remove burnin
    if(!any(x$chains$step>burnin)) stop("burnin too high - no chain left")
    dat <- x$chains[x$chains$step>burnin,,drop=FALSE]
    ances <- dat[,grep("alpha", names(x$chains)),drop=FALSE] # table of ancestries
    Tinf <-  dat[,grep("Tinf", names(x$chains)),drop=FALSE]
    D <- as.matrix(x$D)
    n <- ncol(ances)


    ## AUXILIARY FUNCTIONS, GET RESULTS FOR ONE SAMPLE ##
    f1 <- function(vecAnces, vecTinf){
        ## get the number of mutations between cases ##
        vecAnces <- as.integer(vecAnces)
        vecAnces[vecAnces<1] <- NA
        nmut <- sapply(1:n, function(i) D[i,vecAnces[i]])

        ## get the time between cases ##
        vecTinf <- as.integer(vecTinf)
        deltaT <- vecTinf-vecTinf[vecAnces]

        ## get mutation rate ##
        out <- mean(nmut/deltaT, na.rm=TRUE)
        return(out)
    }

    ## GET MUTATION RATES FROM POSTERIOR SAMPLES ##
    out <- unlist(lapply(1:nrow(ances), function(i) f1(ances[i,],Tinf[i,])))

    ## rescale mutation rate if necessary ##
    if(!is.null(genome.size)) out <- out/genome.size
    return(out)
} # end get.mu
