
#' Select 'good' runs from independent MCMC chains
#'
#' The function \code{selectChains} is used to discard 'bad' MCMC chains from
#' outbreaker's ouptput (functions \code{outbreaker} and
#' \code{outbreaker.parallel}). This is useful whenever several chains were run
#' and converged towards different posterior modes or distributions.  This can
#' happen for instance when imported cases are hard to disentangle, resulting
#' in different runs identifying different imports and therefore having
#' different likelihood.
#'
#' Three modes are available, depending on the argument \code{select} (see also
#' arguments below): \itemize{ \item \code{visual}: (default) interactive mode
#' plotting the log-posterior values for the different chains and asking the
#' user to identify runs to be discarded.  \item \code{auto}: an automatic
#' procedure is used to discard 'bad' runs; see details.  \item
#' \code{[numbers]}: numbers indicating the runs to be discarded.  }
#'
#' The automatic procedure relies on the following recursive process: \itemize{
#' \item 1. Make the ANOVA of the log-posterior values as a function of the run
#' identifier.  \item 2a. If the P-value is greater than alpha
#' (non-significant), exit.  \item 2b. Otherwise, discard the run with the
#' lowest mean log-posterior value, and go back to 1.  }
#'
#' @aliases selectChains
#'
#' @export
#'
#' @param x the output of \code{outbreaker} or \code{outbreaker.parallel}.
#' @param select a character string matching \code{visual} or \code{auto}, or a
#' vector of integers indicating the runs to be discarded.
#' @param alpha the alpha threshold to be used to the automatic procedure (see
#' details)
#' @param \dots further arguments to be passed to \code{\link{plotChains}}.
#' @return These functions similar objects to the inputs, from which 'bad' runs
#' have been discarded.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
selectChains <- function(x, select="visual", alpha=0.001, ...){
    ## CHECKS ##
    if(!is.list(x)) stop("x should be a list as output by outbreaker / outbreaker.parallel")
    if(x$n.runs==1) { # return x unchanged if only one run
        return(x)
    }
    if(is.character(select)) select <- match.arg(select, c("auto","visual"))


    ## SELECTION BASED ON VISUAL INSPECTION OF THE CHAINS ##
    if(is.character(select) && select=="visual"){
        repeat{
            plotChains(x, ...)
            if(x$n.runs==1) {
                cat("Only one run left - exiting.\n")
                break
            }
            cat("Indicate the runs to remove from the results (0 for exit): ")
            answer <- suppressWarnings(as.integer(readLines(n = 1)))
            if(is.na(answer) || answer==0 || answer>x$n.runs) break
            x <- selectChains(x, select=setdiff(1:x$n.runs, answer))
            cat(paste("(removed run ", answer,")\n",sep=""))
        }
    }


    ## SELECTION BASED ON AUTOMATIC PROCEDURE ##
    ## idea:
    ## 1) make ANOVA
    ## if significant:
    ##   2a) remove the run with smallest log-post
    ##   2b) go to 1)
    ## if not: exit
    if(is.character(select) && select=="auto"){
        repeat{
            if(x$n.runs==1) {
                cat("Only one run left - exiting.\n")
                break
            }
            pval <- anova(lm(x$chains$post ~ factor(x$chains$run)))$"Pr(>F)"[1]
            if(pval >= alpha) break
            toRemove <- which.min(tapply(x$chains$post, factor(x$chains$run), mean))
            x <- selectChains(x, select=setdiff(1:x$n.runs, toRemove))
            cat(paste("(removed run ", toRemove,")\n",sep=""))
        }
    }


    ## SELECTION PROVIDED AS NUMBERS ##
    if(is.numeric(select) || is.integer(select)){
        x$chains <- x$chains[x$chains$run %in% select, , drop=FALSE]
        x$chains$run <- as.integer(factor(x$chains$run))
        x$n.runs <- length(select)
    }

    return(x)

} # end selectChains
