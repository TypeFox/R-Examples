# singles.R:  plotmo.singles and plotmo.pairs

#------------------------------------------------------------------------------
# Return a vector of indices of predictors for degree1 plots, e.g, c(1,3,4).
# The indices are col numbers in the x matrix.  The caller will sort the
# returned vector and remove duplicates.  The default method simply
# returns the indices of all predictors.  The object specific methods
# typically return only the predictors actually used in the model.

plotmo.singles <- function(object, x, nresponse, trace, all1)
{
    UseMethod("plotmo.singles")
}
plotmo.singles.default <- function(object, x, nresponse, trace, all1)
{
    seq_len(NCOL(x))
}
#------------------------------------------------------------------------------
# Get the pairs of predictors to be displayed in degree2 plots.
# Each row of the returned pairs matrix is the indices of two predictors
# for a degree2 plot. Example (this was returned from plotmo.pairs.rpart):
#
#    1    2
#    1    2
#    2    1
#
# The indices are col numbers in the x matrix.  The caller will remove
# duplicated pairs and re-order the pairs on the order of the predictors
# in the original call to the model function.  The above example will
# become simply
#
#    1    2
#
# It is ok to return NULL or a matrix with zero rows.

plotmo.pairs <- function(object, x, nresponse=1, trace=0, all2=FALSE)
{
    if(all2) # TODO include all1 in this?
        return(get.all.pairs(object, x, trace))
    UseMethod("plotmo.pairs")
}
# Predictors x1 and x2 are considered paired if they appear in
# the formula in forms such as x1:x2 or I(x1*x2) or s(x1,x2)
#
# We use both formula(object) and attr(terms(object), "term.labels").
# formula(object) is necessary for gam formula like "s(x,v1)" because it
# appears in attr(terms,"term.labels") as "x" "v1" (i.e. as unpaired).
# But our rudimentary parsing of the formula is not reliable, so we also
# use the term.labels.  An lm formula like Volume~(Girth*Height2)-Height
# has term.labels "Girth" "Height2" "Girth:Height2"

plotmo.pairs.default <- function(object, x, nresponse, trace, all2)
{
    formula.vars <- NULL
    formula <- try(formula(object), silent=trace < 2)
    if(is.try.err(formula) || is.null(formula))
        trace2(trace,
               "formula(object) failed for \"%s\" object in plotmo.pairs.default\n",
               class(object)[1])
    else {
        trace2(trace, "formula(object) returned %s\n", strip.space(format(formula)))
        # Note that formula() returns a formula with "." expanded.
        # After as.character: [1] is "~", [2] is lhs, and [3] is rhs
        form <- as.character(formula(object))[3] # rhs of formula
        if(grepl("\\-", form)) { # "-" in formula?
            # formula() gives "(Girth + Height)-Height" for Volume~.-Height
            form <- sub("\\-.*", "", form)    # delete "-" and all after
            form <- gsub("\\(|\\)", "", form) # delete ( and )
        }
        formula.vars <- unlist(strsplit(form, "+", fixed=TRUE))
        formula.vars <- strip.space(formula.vars)
        trace2(trace, "formula.vars %s\n", quotify.trunc(formula.vars))
    }
    term.labels <- NULL
    terms <- try(terms(object), silent=trace < 2)
    if(is.try.err(terms) || is.null(terms))
        trace2(trace,
               "terms(object) failed for \"%s\" object in plotmo.pairs.default\n",
               class(object)[1])
    else {
        term.labels <- attr(terms, "term.labels")
        if(is.null(term.labels))
            trace2(trace,
                   "attr(terms,\"term.labels\") is NULL in plotmo.pairs.default\n")
        else
            trace2(trace, "term.labels %s\n", quotify.trunc(term.labels))
    }
    if(is.null(formula.vars) && is.null(term.labels))
        return(NULL)
    plotmo.pairs.from.term.labels(c(formula.vars, term.labels), colnames(x), trace)
}
get.all.pairs <- function(object, x, trace)
{
    singles <- plotmo.singles(object, x, nresponse=1, trace, all1=TRUE)
    if(length(singles) == 0)
        return(NULL) # no pairs (must be an intercept only model)
    singles <- unique(singles)
    max <- 7 # note that 7 * 6 / 2 is 21 plots
    if(length(singles) > max) {
        warning0("too many predictors to plot all pairs,\n         ",
                 "so plotting degree2 plots for just the first ",
                 max, " predictors")
        singles <- singles[1:max]
    }
    form.pairs(singles)
}
form.pairs <- function(varnames) # return a two column matrix, each row is a pair
{
    col1 <- rep(varnames, times=length(varnames))
    col2 <- rep(varnames, each=length(varnames))
    pairs <- cbind(col1, col2)
    pairs[col1 != col2, , drop=FALSE]
}
# Given the term.labels, return a npairs x 2 matrix specifying which predictors
# are paired. The elements in the returned matrix are column indices of x.
#
# This routine is not infallible but works for the commonly used formulas.
# It works by extracting substrings in each term.label that looks like a
# predictor pair.  The following combos of x1 and x2 for example are
# considered pairs: x1*x2, x1:x2, s(x1,x2), and similar.
#
# TODO this would probably be done best by processing the parse tree

plotmo.pairs.from.term.labels <- function(term.labels, pred.names, trace)
{
    trace2(trace, "plotmo.pairs.from.term.labels\n")
    trace2(trace, "term.labels: %s\n", quotify.trunc(term.labels), "\n")
    trace2(trace, "pred.names:  %s\n", quotify.trunc(pred.names), "\n")
    pairs <- matrix(0, nrow=0, ncol=2)          # no pairs initially
    for(i in 1:length(term.labels)) {
        s <- strip.space(term.labels[i])
        s <- gsub("[+*/,]", ":", s)             # replace + * / , with :
        s <- gsub("=[^,)]+", "", s)             # delete "=any"

        # get the indices of expressions of the form "ident1:ident2"
        igrep <- gregexpr("[._$[:alnum:]]+:[._$[:alnum:]]+", s)[[1]]

        trace2(trace, "considering %s", s)

        if(igrep[1] > 0) for(i in seq_along(igrep)) {
            # extract the i'th "ident1:ident2" into pair
            start <- igrep[i]
            stop <- start + attr(igrep, "match.length")[i] - 1
            pair <- substr(s, start=start, stop=stop)
            pair <- strsplit(pair, ":")[[1]]    # pair is now c("ident1","ident2")
            # are the variables in the candidate pair in pred.names?
            ipred1 <- which(pred.names == pair[1])
            ipred2 <- which(pred.names == pair[2])
            trace2(trace, " ->%s%s",
                if(length(ipred1))
                    sprintf(" %g=%s", ipred1, pred.names[ipred1]) else "",
                if(length(ipred2))
                    sprintf(" %g=%s", ipred2, pred.names[ipred2]) else "")
            if(length(ipred1) == 1 && length(ipred2) == 1 && pair[1] != pair[2])
                pairs <- c(pairs, ipred1, ipred2)
        }
        trace2(trace, "\n")
    }
    matrix(pairs, ncol=2, byrow=TRUE)
}
