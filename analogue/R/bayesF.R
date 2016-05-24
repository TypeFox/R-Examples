## compute Bayes factors (syn. likelihood ratios) of
## positive and negative events
bayesF <- function(x, prior = rep(0.5, 2)) {
    ## GLS: Need to get this to work with the priors stored in each
    ## roc object, OR take a new prior, either a global one for all
    ## groups, or individual priors. both supplied by the user here
    FUN <- function(x, prior) {
        pos <- x$TPF / x$FPE
        neg <- (1 - x$TPF) / (1 - x$FPE)
        ## should we correct? if so what to?
        ##pos[is.infinite(pos)] <- 1
        ##neg[is.nan(neg)] <- 0
        ## prior probs
        prior.prob.pos <- prior[1]
        prior.prob.neg <- prior[2]
        ## prior odds
        prior.odds.pos <- prior.prob.pos / (1 - prior.prob.pos)
        prior.odds.neg <- prior.prob.neg / (1 - prior.prob.neg)
        ## posterior odds
        post.odds.pos <- pos * prior.odds.pos
        post.odds.neg <- neg * prior.odds.neg
        ## posterior probabilities
        post.probs.pos <- post.odds.pos / (1 + post.odds.pos)
        post.probs.neg <- post.odds.neg / (1 + post.odds.neg)
        ## return object
        retval <- list(bayesF = list(pos = pos, neg = neg),
                       posterior.odds = list(pos = post.odds.pos,
                       neg = post.odds.neg),
                       posterior.probs = list(pos = post.probs.pos,
                       neg = post.probs.neg),
                       prior.prob = list(pos = prior.prob.pos,
                       neg = prior.prob.neg),
                       roc.points = unique(x$roc.points),
                       optimal = x$optimal,
                       max.roc = x$max.roc)
    }
    if(is.null(prior))
        prior <- rep(0.5, 2)
    retval <- lapply(x$roc, FUN, prior = prior)
    names(retval) <- names(x$roc)
    retval$object <- deparse(substitute(x))
    retval$prior <- prior
    class(retval) <- "bayesF"
    attr(retval, "method") <- attr(x, "method")
    return(retval)
}

print.bayesF <- function(x, digits = min(3, getOption("digits") - 4),
                         ...) {
    cat("\n")
    writeLines(strwrap("Bayes factors (likelihood ratios)",
                       prefix = "\t"))
    cat("\n")
    cat(paste("Object:", x$object, "\n"))
    ## groups names
    gnames <- names(x)
    gnames <- gnames[!gnames %in% c("object","prior", "Combined")]
    ##gnames <- matrix(gnames, nrow = 1)
    ##names(gnames) <- paste("Group", seq_along(gnames))
    cat("\n")
    writeLines(strwrap(paste("Groups (N = ", length(gnames), "):",
                             sep = "")))
    writeLines(strwrap(paste(gnames, collapse = ", ", sep = ""),
                       prefix = "  "))
    cat("\n")
    cat("\nPrior probabilities:\n")
    cat(paste("Positive:", round(x$prior[1], digits),
              "  Negative:", round(x$prior[2], digits), "\n"))
    cat("\n")
    invisible(x)
}
