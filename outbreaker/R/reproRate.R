
#' Derive reproduction numbers from outbreak's outputs
#'
#' These functions are used to compute reproduction numbers and derive
#' incidence curves from outbreaker's ouptput (functions \code{outbreaker} and
#' \code{outbreaker.parallel}). They all rely on the entire outbreak having
#' been sampled.  \itemize{ \item \code{get.R} derive distributions of
#' individual effective reproduction numbers.  \item \code{get.Rt} derives
#' effective reproduction numbers averaged for each time step.  \item
#' \code{get.incid} derives incidence curves for each time step.  }
#'
#' @export
#'
#' @aliases get.R get.Rt get.incid
#'
#' @rdname repro
#'
#' @param x the output of \code{outbreaker} or \code{outbreaker.parallel}.
#' @param burnin an integer indicating the number of steps of the MCMC to be
#' discarded as burnin period. Defaults to 20,000.
#' @param plot a logical indicating whether a plot should be displayed.
#' @param type a character indicating the type of plot to be used.
#' @param lines a logical indicating whether individual lines should be added
#' to the plot.
#' @param fill.col the color to be used for the boxplot.
#' @param lines.col the color to be used to the lines.
#' @param \dots further arguments to be passed to other functions.
#' @return These functions return a \code{data.frame} containing the plotted
#' information.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @examples
#'
#' ## load data
#' data(fakeOutbreak)
#' attach(fakeOutbreak)
#'
#' ## individual R
#' barplot(table(get.R(res)), main="Individual effective reproduction numbers")
#'
#' ## R(t)
#' get.Rt(res)
#'
#' ## incidence
#' get.incid(res)
#'
#' detach(fakeOutbreak)
#'
#'
get.Rt <- function(x, burnin=2e4, plot=TRUE, type=c("boxplot", "lines"), lines=FALSE,
                   fill.col="gold", lines.col=transp("grey"), ...){
    ## if(!require(adegenet)) stop("the adegenet package is required.")
    type <- match.arg(type)

    ## GET DATA ##
    ## remove burnin
    if(!any(x$chains$step>burnin)) stop("burnin too high - no chain left")
    dat <- x$chains[x$chains$step>burnin,,drop=FALSE]

    ## get ancestries
    ances <- dat[,grep("alpha", names(x$chains)),drop=FALSE] # table of ancestries
    tabAnces <- apply(ances, 1, table) # count nb of descendents per case for each chain

    ## get infection times
    Tinf <-  dat[,grep("Tinf", names(x$chains)),drop=FALSE]
    timeSpan <- range(Tinf)
    timeStep <- seq(timeSpan[1],timeSpan[2],by=1)
    emptyOut <- rep(NA, length(timeStep))
    names(emptyOut) <- timeStep

    ## function to get Rt for one chain 'i'
    f1 <- function(i){
        ## get nb of descendents per ancestor
        if(is.list(tabAnces)){
                e <- tabAnces[[i]][-1] # -1: remove '0's
        } else {
            e <- tabAnces[-1,i] # -1: remove '0's
        }

        ## create empty output
        out <- emptyOut

        ## find time steps with at least one new case, set default R to 0
        out[as.character(unique(as.integer(Tinf[i,])))] <- 0

        ## get infection times of infectors
        Tinf.temp <- Tinf[i,as.integer(names(e))]

        ## count mean nb of secondary infections created by cases infected at each time step
        meanNbCasePerTimeStep <- tapply(as.numeric(e),as.numeric(Tinf.temp),mean)

        out[names(meanNbCasePerTimeStep)] <- meanNbCasePerTimeStep
        return(out)
    }

    ## GET RT FOR ALL RELEVANT CHAINS ##
    res <- lapply(1:nrow(dat), function(i) f1(i))
    res <- t(Reduce("rbind",res))


    ## MAKE PLOT IF NEEDED ##
    if(plot){
        ## if(type=="CI"){
        ##     ## CI-based
        ##     stat <- apply(res, 1, quantile, c(CI.level,0.5, 1-CI.level), na.rm=TRUE)
        ##     xcoord <- as.numeric(colnames(stat))
        ##     matplot(xcoord, t(stat), type="n", ...)
        ##     polygon(c(xcoord,rev(xcoord)), c(stat[3,], rev(stat[1,])), col=fill.col, border=NA)
        ##     if(lines){
        ##         matplot(rownames(res),res, type="l", lty=1, col=lines.col, lwd=1, add=TRUE)
        ##     }
        ##     matplot(xcoord, t(stat), type="l", lty=1, col="black", lwd=2, add=TRUE)
        ## }



        ## boxplot-based
        if(type=="boxplot"){
            boxplot(t(res), col=fill.col, at=as.integer(rownames(res)),
                    xlab="Time", ylab="Effective reproduction number (R(t))", ...)
            if(lines){
                matplot(rownames(res),res, type="l", lty=1, col=lines.col, lwd=1, add=TRUE)
            }
        }

        ## just lines
        if(type=="lines"){
            matplot(rownames(res),res, type="l", lty=1, col=lines.col,
                    xlab="Time", ylab="Effective reproduction number (R(t))", ...)
        }
    }

    return(res)
} # end get.Rt




#' @rdname repro
#' @export
#'
get.R <- function(x, burnin=2e4, ...){

    ## GET DATA ##
    ## remove burnin
    if(!any(x$chains$step>burnin)) stop("burnin too high - no chain left")
    chains <- x$chains[x$chains$step>burnin,,drop=FALSE]

    ## get ancestries
    ances <- chains[,grep("alpha", names(x$chains)),drop=FALSE] # table of ancestries

    ## function to get R for a given step
    id <- 1:length(x$collec.dates)
    f1 <- function(vecAnces){
        return(sapply(id, function(i) sum(vecAnces==i,na.rm=TRUE)))
    }

    res <- t(apply(ances,1,f1))
    colnames(res) <- id
    return(res)
} # end get.R





#' @rdname repro
#' @export
#'
get.incid <- function(x, burnin=2e4, plot=TRUE, type=c("boxplot", "lines"), lines=FALSE,
                      fill.col="gold", lines.col=transp("grey"), ...){
    ## if(!require(adegenet)) stop("the adegenet package is required.")
    type <- match.arg(type)

    ## GET DATA ##
    ## remove burnin
    if(!any(x$chains$step>burnin)) stop("burnin too high - no chain left")
    dat <- x$chains[x$chains$step>burnin,,drop=FALSE]

    ## get infection times
    Tinf <-  as.matrix(dat[,grep("Tinf", names(x$chains)),drop=FALSE])
    timeSpan <- range(Tinf)
    timeStep <- seq(timeSpan[1],timeSpan[2],by=1)
    emptyOut <- rep(0, length(timeStep))
    names(emptyOut) <- timeStep

    ## function to get Rt for one chain
    f1 <- function(i){
        ## get nb of descendents per ancestor
        e <- Tinf[i,,drop=TRUE] # -1: remove '0's

        ## create empty output
        out <- emptyOut

        ## fill in output
        out[names(table(e))] <- table(e)

        return(out)
    }

    ## GET RT FOR ALL RELEVANT CHAINS ##
    res <- lapply(1:nrow(dat), function(i) f1(i))
    res <- t(Reduce("rbind",res))


    ## MAKE PLOT IF NEEDED ##
    if(plot){
        ## if(type=="CI"){
        ##     ## CI-based
        ##     stat <- apply(res, 1, quantile, c(CI.level,0.5, 1-CI.level), na.rm=TRUE)
        ##     xcoord <- as.numeric(colnames(stat))
        ##     matplot(xcoord, t(stat), type="n", ylab="Incidence", ...)
        ##     polygon(c(xcoord,rev(xcoord)), c(stat[3,], rev(stat[1,])), col=fill.col, border=NA)
        ##     if(lines){
        ##         matplot(rownames(res),res, type="l", lty=1, col=lines.col, lwd=1, add=TRUE)
        ##     }
        ##     matplot(xcoord, t(stat), type="l", lty=1, col="black", lwd=2, add=TRUE)
        ## }



        ## boxplot-based
        if(type=="boxplot"){
            boxplot(t(res), col=fill.col, at=as.integer(rownames(res)),
                    xlab="Time", ylab="Incidence", ...)
            if(lines){
                matplot(rownames(res),res, type="l", lty=1, col=lines.col, lwd=1, add=TRUE)
            }
        }

        ## just lines
        if(type=="lines"){
            matplot(rownames(res),res, type="l", lty=1, col=lines.col,
                    xlab="Time", ylab="Incidence", ...)
        }
    }

    return(res)
} # end get.incid
