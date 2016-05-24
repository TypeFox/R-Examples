# Canonical discriminant analysis for all terms in an mlm

# Written by: John Fox
# Last updated: 10/29/2007 9:33:42 AM MF
# --- Changed default to ask=interactive() in plot.candiscList

# Last modified: 31 August 2007 by J. Fox


candiscList <- function(mod, ...){
    UseMethod("candiscList")
    }
    
candiscList.mlm <- function(mod, type="2", manova, ndim, ...){
    if (missing(manova)) manova <- Anova(mod, type=type)
    result <- as.list(1:length(manova$terms))
    names(result) <- manova$terms
    for (term in manova$terms){
        result[[term]] <- if (missing(ndim))
            candisc(mod, type=type, manova=manova, term=term)
            else {
                nd <- if(is.list(ndim)) ndim$term
                    else ndim 
                candisc(mod, type=type, manova=manova, ndim=nd, term=term)
                }
        }
    class(result) <- "candiscList"
    result
    }

print.candiscList <- function(x, ...){
    terms <- names(x)
    for (term in terms){
        cat("\nTerm:", term, "\n")
        print(x[[term]], ...)
        cat("\n")
        }
    }

summary.candiscList <- function(object, ...){
    terms <- names(object)
    for (term in terms){
        cat("\nTerm:", term, "\n")
        summary(object[[term]], ...)
        cat("\n")
        }
    }

plot.candiscList <- function(x, term, ask=interactive(), graphics = TRUE, ...) {
    if (!missing(term)){
        if (is.character(term)) term <- gsub(" ", "", term)
        plot(x[[term]], ...)
        return(invisible())
        }
    terms <- names(x)
    if (ask){
        repeat {
            selection <- menu(terms, graphics = graphics, title = "Select term to plot")
            if (selection == 0) break
            else plot(x[[selection]], ...)
            }
        }
    else {
        nterms <- length(x)
        for (i in 1:nterms) {
        	plot(x[[i]], ...)
        	}
        }
}
