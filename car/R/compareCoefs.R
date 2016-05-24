# 21 May 2010: small changes to output when there is just one model. J. Fox 
# 15 Aug 2010: changed name of function to compareCoefs to avoid name clash. J. Fox
# 18 May 2011: check for 'mer' objects, and handle them correctly. S. Weisberg
#  8 Sep 2011: check for 'lme' objects, and handle them correctly. S. Weisberg
# 11 Jan 2012: fix to work with any 'S4' object with a coef() method. 
#   suggested by David Hugh-Jones  University of Warwick http://davidhughjones.googlepages.com 
# 3 May 2012:  fixed bug if models are less than full rank.
# 17 Sept 2012: suppressing printing calls when there are none. J. Fox
# 22 June 2013: tweaks for lme4. J. Fox
# 26 Sept 2014: cleaned up printing of calls. J. Fox

compareCoefs <- function (..., se = TRUE, print = TRUE, digits = 3) {
    splitExpr <- function(expr, width=getOption("width") - 4, at="[ ,=]"){
        if (length(grep("\n", expr)) >0 ){
            cmds <- strsplit(expr, "\n")[[1]]
            allcmds <- character(length(cmds))
            for (i in 1:length(cmds))
                allcmds[i] <- splitExpr(cmds[i], width=width, at=at)
            return(paste(allcmds, collapse="\n"))
        }
        if (nchar(expr) <= width) return(expr)
        where <- gregexpr(at, expr)[[1]]
        if (where[1] < 0) return(expr)
        singleQuotes <- gregexpr("'", expr)[[1]]
        doubleQuotes <- gregexpr('"', expr)[[1]]
        comment <- regexpr("#", expr)
        if (singleQuotes[1] > 0 && (singleQuotes[1] < doubleQuotes[1] || doubleQuotes[1] < 0 ) && (singleQuotes[1] < comment[1] || comment[1] < 0 )){
            nquotes <- length(singleQuotes)
            if (nquotes < 2) stop("unbalanced quotes")
            for(i in seq(nquotes/2))
                where[(where > singleQuotes[2 * i - 1]) & (where < singleQuotes[2 * i])] <- NA
            where <- na.omit(where)
        }  
        else if (doubleQuotes[1] > 0 && (doubleQuotes[1] < singleQuotes[1] || singleQuotes[1] < 0) && (doubleQuotes[1] < comment[1] || comment[1] < 0 )){
            nquotes <- length(doubleQuotes)
            if (nquotes < 2) stop("unbalanced quotes")
            for(i in seq(nquotes/2))
                where[(where > doubleQuotes[2 * i - 1]) & (where < doubleQuotes[2 * i])] <- NA
            where <- na.omit(where)
        }
        else if (comment > 0){
            where[where > comment] <- NA
            where <- na.omit(where)
        }
        if (length(where) == 0) return(expr)
        where2 <- where[where <= width]
        where2 <- if (length(where2) == 0) where[1]
        else where2[length(where2)]
        paste(substr(expr, 1, where2), "\n  ", 
            Recall(substr(expr, where2 + 1, nchar(expr)), width, at), sep="")
    } 
    removeExtraQuotes <- function(string) sub("\\\"$", "", sub("^\\\"", 
        "", string))
    squeezeMultipleSpaces <- function(string) gsub(" {2,}", " ", string)
    models <- list(...)
    n.models <- length(models)
    if (n.models < 1) 
        return(NULL)
    getnames <- function(model) {
        if (inherits(model, "merMod") || inherits(model, "mer") | 
                inherits(model, "lme")) 
            names(fixef(model))
        else names(coef(model))
    }
    getcoef <- function(model) {
        if (inherits(model, "merMod") || inherits(model, "mer") | 
                inherits(model, "lme")) 
            fixef(model)
        else coef(model)
    }
    getcall <- function(model) {
        paste(deparse(if (isS4(model))  model@call else model$call), collapse="")
    }
    getvar <- function(model) {
        if (inherits(model, "merMod") || inherits(model, "mer")) 
            as.matrix(vcov(model))
        else vcov(model)
    }
    coef.names <- unique(unlist(lapply(models, getnames)))
    table <- matrix(NA, length(coef.names), n.models * (1 + se))
    rownames(table) <- coef.names
    colnames(table) <- if (se) 
        if (n.models > 1) 
            paste(rep(c("Est.", "SE"), n.models), rep(1:n.models, 
                each = 2))
    else c("Estimate", "Std. Error")
    else if (n.models > 1) 
        paste(rep("Est.", n.models), 1:n.models)
    else "Estimate"
    calls <- !any(sapply(models, getcall) == "NULL")
    if (print == TRUE && calls) 
        cat("\nCall:")
    for (i in 1:n.models) {
        model <- models[[i]]
        fout <- getcall(model)
        mod <- if (n.models > 1) 
            paste(i, ": ", sep = "")
        else ""
        if (print && calls) 
            cat(splitExpr(squeezeMultipleSpaces(paste("\n", mod, removeExtraQuotes(fout[1]), 
                sep = ""))))
        if (print && calls && length(fout) > 1) 
            for (f in fout[-1]) cat("\n", splitExpr(squeezeMultipleSpaces(removeExtraQuotes(f))))
        if (se) {
            ests <- getcoef(model)
            new <- cbind(ests, rep(NA, length(ests)))
            new[!is.na(ests), 2] <- sqrt(diag(getvar(model)))
            table[getnames(model), 2 * (i - 1) + c(1, 2)] <- new
        }
        else table[getnames(model), i] <- getcoef(model)
    }
    if (print == TRUE) {
        cat("\n")
        printCoefmat(table, na.print = "", digits = digits, tst.ind = NULL)
    }
    else table
}


