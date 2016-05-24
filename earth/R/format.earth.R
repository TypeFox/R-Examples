# format.earth.R

# Return a vector s of strings, length of vector is nresponses.
# But if there are embedded GLM model(s) then length of s is 2 * nresponses
# and strings for the GLM model(s) start at s[nresponses].
#
# For each model string, there is one term per line.  Each term (except
# the intercept) is made up of a coefficent which multiplies one or more
# hockey stick funcs.
#
# For the default style="h" the result looks like this:
#
#         23
#         +  5.7 * h(Girth-12.9)
#         -  2.9 * h(12.9-Girth)
#         + 0.72 * h(Height-76)
#
# For style="pmax" the result looks like this:
#
#         23.208244
#         +  5.7459616 * pmax(0,  Girth -   12.9)
#         -  2.8664516 * pmax(0,   12.9 -  Girth)
#         + 0.71833643 * pmax(0, Height -     76)
#
# Style="max" is the same as "pmax" but prints "max" rather than "pmax".
#
# Style="C" looks like this:
#         23.208244
#         +  5.7459616 * max(0  x[0] - 12.9)
#         -  2.8664516 * max(0  12.9 - x[0])
#         + 0.71833643 * max(0, x[1] - 76)
#
# For style="bf" the result looks like this:
#
#         bf1: h(Girth-12.9)
#         bf2: h(12.9-Girth)
#         bf3: h(Height-76)
#
#         23
#         +  5.7 * bf1
#         -  2.9 * bf2
#         + 0.72 * bf3
#
# decomp argument: see reorder.earth()
#
# The first arg is actually an object but called x for consistency with generic
#
# TODO would be nice to add an option to print in term importance order

format.earth <- function(
    x           = stop("no 'x' argument"),   # "earth" object
    style       = "h",      # see get.term.strings
    decomp      = "anova",  # see reorder.earth for legal decomp values
    digits      = getOption("digits"),
    use.names   = TRUE,
    colon.char  = ":",      # convert colons in expression to this char
    ...)                    # unused, for consistency with generic
{
    check.classname(x, substitute(x), "earth")
    warn.if.dots(...)
    use.names <- check.boolean(use.names)
    nresp <- NCOL(x$coefficients)
    s <- vector(mode = "character", length=nresp)
    if(style[1] == "C") {
        if(digits < 5)
            digits <- 5
        use.names <- -1 # tell variable.names.earth to use zero based indexing
        colon.char = "*"
    }
    for(iresp in seq_len(nresp))
        s[iresp] <- format.one.response(iresp, x, digits, use.names,
                                        decomp, style=style,
                                        colon.char=colon.char, coefs=NULL)

    if(!is.null(x$glm.list))   # embedded GLM model(s)?
        for(iresp in seq_len(nresp))
            s[nresp + iresp] <- format.one.response(iresp, x, digits, use.names,
                                        decomp, style=style,
                                        colon.char=colon.char,
                                        coefs=x$glm.list[[iresp]]$coefficients)
    s
}
format.one.response <- function( # called by format.earth
    iresp,          # response index i.e. column in y matrix
    object,         # "earth" object
    digits,
    use.names,
    decomp,         # see reorder.earth for legal decomp values
    style,          # see get.term.strings
    colon.char=":", # convert colons in output to this char
    coefs=NULL)     # if not NULL use these instead of object$coefficients
{
    new.order <- reorder.earth(object, decomp=decomp)
    if(is.null(coefs))
        coefs <- object$coefficients[, iresp]
    coefs <- coefs[new.order]
    which.terms <- object$selected.terms[new.order]
    dirs <- object$dirs
    check.which.terms(dirs, which.terms)
    term.names <- get.term.strings(object, digits, use.names, style, new.order)
    # convert colons to colon.char
    term.names <- make.unique(gsub(":", colon.char, term.names), sep="_")
    coef.width <- get.coef.width(coefs[-1], digits)
    s <- if(style[1] == "C") "" else "  "  # result goes into this string
    s <- pastef(s, "%.*g\n", digits=digits, coefs[1])
    iterm <- 2
    while(iterm <= length(which.terms)) {
        coef <- coefs[iterm]
        if(coef < 0)
            s <- pastef(s, "  - %s ",
            format(-coef, justify="left",w=coef.width,digits=digits,format="%g"))
        else
            s <- pastef(s, "  + %s ",
                        format(coef, justify="left",
                               w=coef.width,digits=digits,format="%g"))
        s <- pastef(s, "* %s", term.names[iterm])
        s <- pastef(s, "\n")
        iterm <- iterm + 1
    }
    if(pmatch(style, "bf", nomatch=0)) # append table of basis functions?
        s <- paste0(s, "\n", get.table.of.basis.functions(object, new.order))
    s
}
get.coef.width <- function(coefs, digits)   # get print width for earth coefs
{
    if(length(coefs) > 0)
        max(nchar(format(abs(coefs), digits=digits)))
    else
        10 # arbitrary width if no coefs
}
# style argument:
#   "h"    gives "h(survived-0) * h(16-age)"
#   "pmax" gives "pmax(0, survived - 0) * pmax(0, 16 - age)"
#   "max"  gives "max(0, survived - 0) * max(0, 16 - age)"
#   "C"    gives "max(0, x[0]) * max(0, 16 - x[1])"
#   "bf"   gives basis functions e.g. "bf1" or "bf1 * bf3"

get.term.strings <- function(object, digits, use.names,
                             style = c("h", "pmax", "max", "C", "bf"),
                             neworder)
{
    switch(match.arg1(style, "style"),
        "h"    = get.term.strings.h(object, digits, use.names, neworder),
        "pmax" = get.term.strings.pmax(object, digits, use.names, neworder, "pmax"),
        "max"  = get.term.strings.pmax(object, digits, use.names, neworder, "max"),
        "C"    = get.term.strings.pmax(object, digits, use.names, neworder, "max"),
        "bf"   = get.term.strings.bf(object,   digits, use.names, neworder))
}
get.term.strings.h <- function(object, digits, use.names, new.order)
{
    # digits is unused

    if(!use.names)
        warning("use.names=FALSE ignored because style=\"h\"")

    s <- colnames(object$bx)[new.order]
}
# TODO need to add factor simplification to this routine
# TODO need to add get.ndigits functionality (in get.earth.term.name) to this routine

get.term.strings.pmax <- function(object, digits, use.names, new.order, fname)
{
    # get.width returns the width for printing elements of the earth expression.
    # This is used to keep things lined up without too much white space.
    # This returns the widest of all possible printed elements.

    get.width <- function(which.terms, dirs, var.names, cuts, digits)
    {
        if(length(which.terms) == 1)
            return(10)  # return arbitrary width for intercept only model
        used.dirs <- dirs[which.terms, , drop=FALSE]
        # used.preds is a logical index vector which selects used x predictors
        used.preds <- apply(used.dirs, 2, any1)
        # as.list is needed so format treats each cut independently
        max(nchar(var.names[used.preds]),
            nchar(format(as.list(cuts[which.terms, used.preds]), digits=digits)))
    }
    which.terms <- object$selected.terms[new.order]
    cuts <- object$cuts
    var.names <- variable.names.earth(object, use.names=use.names)
    which.terms <- object$selected.terms[new.order]
    dirs <- object$dirs
    width <- get.width(which.terms, dirs, var.names, cuts, digits)
    nterms <- length(which.terms)
    s <- character(nterms)
    s[1] = "(Intercept)"
    iterm <- 2
    fname <- if(fname=="h") "h(" else paste0(fname, "(0, ")
    while(iterm <= nterms) {
        isel.term <- which.terms[iterm]
        dir <- dirs[isel.term, , drop=FALSE]
        cut <- cuts[isel.term, , drop=FALSE]
        npreds <- ncol(cuts)
        prefix <- ""
        for(ipred in seq_len(npreds)) {
            if(dir[ipred]) {
                if(dir[ipred] == 2)     # linear predictor?
                    s[iterm] <- pastef(s[iterm], "%s%-*s %*s            ",
                                    prefix, width=width,
                                    var.names[ipred], width=width, "")
                else if(dir[ipred] == -1)
                    s[iterm] <- pastef(s[iterm], "%s%s%s - %*s) ",
                                    prefix, fname,
                                    format(cut[ipred], width=width, digits=digits),
                                    width, var.names[ipred])
                else if(dir[ipred] == 1)
                    s[iterm] <- pastef(s[iterm], "%s%s%*s - %s) ",
                                 prefix, fname,
                                 width=width, var.names[ipred],
                                 format(cut[ipred], width=width, digits=digits))
                else
                    stop0("illegal direction ", dir[ipred], " in 'dirs'")

                prefix <- "* "
            }
        }
        iterm <- iterm + 1
    }
    s
}
# return a data.frame, each row has 2 elements: the original and new basis function names

get.bfs <- function(names)
{                                           # Example: start of with names:
                                            # "(Intercept)", "h(temp-58)", "h(humidity-55)*h(temp-58)", ...
    # make a single long string
    s0 <- paste(names, collapse="")         # "(Intercept)h(temp-58)h(humidity-55)*h(temp-58)..."
    # replace * with nothing
    s1 <- gsub("*", "", s0, fixed=TRUE)     # "(Intercept)h(temp-58)h(humidity-55)h(temp-58)..."
    # replace ) with )@ so @ are split points
    s2 <- gsub(")", ")@", s1, fixed=TRUE)   # "(Intercept)@h(temp-58)@h(humidity-55)@h(temp-58)..."
    # separate strings at split points
    s3 <- strsplit(s2, split="@")[[1]]      # "(Intercept)", "h(temp-58)", "h(humidity-55)" "h(temp-58)", ...
    # remove duplicate strings, result is a vector of all basis function names
    original <- unique(s3)                  # "(Intercept)", "h(temp-58)", "h(humidity-55)" "h(temp-58)", ...

    # -1 below so first term is bf1 (i.e. intercept is bf0, which is unused)
    new <- paste0("bf", seq_along(original)-1) # "bf1", "bf2", "bf3", ...

    data.frame(original, new)
}
get.term.strings.bf <- function(object, digits, use.names, new.order)
{
    # digits is unused

    if(!use.names)
        warning("use.names=FALSE ignored because style=\"h\"")

    names <- colnames(object$bx)[new.order] # "(Intercept)", "h(temp-58)", "h(humidity-55)*h(temp-58)", ...
    bfs <- get.bfs(names)

    # replace original names with names in new.bfs

    if(nrow(bfs) > 1) for(i in 2:nrow(bfs)) # start at 2 to skip intercept
        names <- gsub(bfs[i,1], bfs[i,2], names, fixed=TRUE)

    gsub("*", " * ", names, fixed=TRUE)     # put space around *
}
# Return a string like this:
#   bf1  h(temp-58)
#   bf2  h(234-ibt)
#   bf3  h(200-vis)
#   bf4  h(doy-92)

get.table.of.basis.functions <- function(object, new.order)
{
    names <- colnames(object$bx)[new.order]
    bfs <- get.bfs(names)
    s <- ""
    if(nrow(bfs) > 1) for(i in 2:nrow(bfs)) # start at 2 to skip intercept
        s <- paste0(s, sprintf("%6s  %s\n", bfs[i,2], bfs[i,1]))
    s
}
# Return a string representing the linear model.
# Example: a <- lm(Volume ~ ., data = trees); cat(format(a))
# which yields:
#
#   -58
#   +  4.71 * Girth
#   + 0.339 * Height
#
# The first arg is actually an object but called x for consistency with generic
#
# TODO this function doesn't really belong in the earth package

format.lm <- function(
    x           = stop("no 'x' argument"),   # "lm" object, also works for "glm" objects
    digits      = getOption("digits"),
    use.names   = TRUE,
    colon.char  = ":",                  # convert colons in expression to this char
    ...)                                # unused, for consistency with generic
{
    format1 <- function(coef)
    {
        format(coef, justify="left", w=coef.width, digits=digits, format="%g")
    }
    check.classname(x, substitute(x), "lm")
    use.names <- check.boolean(use.names)
    dataClasses <- attr(x$terms, "dataClasses")
    # TODO extend this function to handle factors
    if(any((dataClasses == "factor") | (dataClasses == "ordered")))
        stop0("a predictor has class \"factor\" and format.lm cannot handle that")
    coefs <- coef(x)
    stopifnot(length(coefs) > 0)
    if(!is.vector(coefs) || NCOL(coefs) > 1)
        stop0("format.lm can only handle single response models")
    pred.names <- names(coefs)
    if(is.null(pred.names) || !use.names) {
        pred.names <- paste("x[,", 0:length(coefs), "]", sep="")
        pred.names[1] <- "(Intercept)"
    }
    pred.names <- make.unique(gsub(":", colon.char, pred.names), sep="_")
    # if any coef is NA then issue warning and change the coef to 0
    if(anyNA(coefs)) {
        which <- which(is.na(coefs))
        warnf("coefficient for %s%s is NA, printing it as 0",
              pred.names[which[1]],
              if(length(which) > 1) " and others" else "")
        coefs[which] <- 0
    }
    intercept <- 0
    intercept.index <- match("(Intercept)", names(coefs), nomatch=0)
    if(intercept.index) {
        stopifnot(intercept.index == 1)
        intercept <- coefs[1]
        pred.names <- pred.names[-1] # drop intercept
        coefs <- coefs[-1]
    }
    s <- sprintf("  %.*g\n", digits=digits, intercept)
    coef.width <- get.coef.width(coefs, digits)
    for(ipred in seq_along(coefs)) {
        coef <- coefs[ipred]
        if(coef < 0)
            s <- pastef(s, "  - %s ", format1(-coef))
        else
            s <- pastef(s, "  + %s ", format1(coef))
        s <- pastef(s, "* %s", pred.names[ipred])
        s <- pastef(s, "\n")
    }
    s
}
