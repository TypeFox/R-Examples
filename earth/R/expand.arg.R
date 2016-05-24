# expand.arg.R:

# contr.earth.response returns an nlevels by nlevels diag matrix e.g.
#
#       A B C
#     A 1 0 0
#     B 0 1 0
#     C 0 0 1
#
# The base and contrasts arguments are ignored

contr.earth.response <- function(x, base, contrasts)
{
     contr <- array(0, c(length(x), length(x)), list(x, x))
     diag(contr) <- 1
     contr
}
# Return x but with factors expanded into dummy variables.
# and with all values converted to double.
# Always returns a matrix (never a vector) and always with column names.
# Factors in earth's y argument are expanded to indicator columns exactly
# like factors in the x argument, except contr.earth.response is always used.

expand.arg <- function(x,              # "x" is x or y arg to earth
                       env,            # evironment for evaluation
                       trace,          # passed to gen.colnames
                       is.y.arg=FALSE, # is.y.arg is TRUE if y arg to earth
                       xname=NULL)     # used for colnames when x has no name
{
    if(is.null(xname))
        xname <- sub(".*\\$", "", trunc.deparse(substitute(x))) # temp$y becomes y
    if(is.null(ncol(x))) # make sure x is a matrix, not a vector
        dim(x) <- c(nrow=length(x), ncol=1)
    if(is.y.arg) {
        # We must do this here else the call to model.matrix later generates
        # two columns for each logical or two-level factor column.
        x <- convert.two.level.to.numeric(x)
    }
    if(is.double(x)) {
        # Already double so no need to convert.  Note that is.double returns
        # TRUE for a matrix of doubles but always FALSE for data.frames.
        colnames(x) <-
            gen.colnames(x, xname, if(is.y.arg) "y" else "x", trace, xname)
        return(x)
    }
    if(is.y.arg) {
        # We always use contr.earth.response for the y argument (left side of formula).
        old.contrasts <- getOption("contrasts")
        on.exit(options(contrasts=old.contrasts))
        options(contrasts=c("contr.earth.response", "contr.earth.response"))
    }
    is.data.frame <- is.data.frame(x)
    if(is.data.frame)
        mf <- call("model.frame", formula = ~., data=x, na.action=na.pass)
    else
        mf <- call("model.frame", formula = ~x, na.action=na.pass)
    mf <- eval(mf, env) # this is slow
    mf.has.colnames <- !is.null(colnames(mf))
    x <- model.matrix(object=attr(mf, "terms"), data=mf)
    intercept <- match("(Intercept)", colnames(x), nomatch=0)
    if(intercept)
        x <- x[, -intercept, drop=FALSE]     # discard intercept

    # If !is.data.frame, model.matrix prepends "x" to the column names,
    # so remove the "x", but only if x had column names to begin with.

    if(!is.data.frame && mf.has.colnames)
        colnames(x) <- substr(colnames(x), 2, 61) # strip 1st char of each colname
    colnames(x) <-
        gen.colnames(x, xname, if(is.y.arg) "y" else "x", trace, xname)

    x   # all columns are now double with column names
}
# Here "two.level" means logical or two-level factor.
# These get converted to a numeric column of 0s and 1s.
# This code doesn't touch y if no changes are needed.

convert.two.level.to.numeric <- function(y)
{
    stopifnot(!is.null(dim(y)))
    if(is.data.frame(y)) {
        # Dataframe, so handle each column independently.
        # Get here if y is a dataframe in call to earth.default.
        for(icol in seq_len(ncol(y))) {
            ycol <- y[,icol]
            if(is.logical(ycol))
                y[,icol] <- as.numeric(ycol)
            else if(is.factor(ycol) && nlevels(ycol) <= 2) # two-level factor?
                y[,icol] <- as.numeric(ycol) - 1
        }
    } else {
        # Not dataframe, must be a matrix.  All columns are of the same class.
        # Can't use above code because for example y[,icol] <- as.numeric(ycol)
        # generates NAs if y is a matrix of factors.
        convert.logical <- convert.factor <- FALSE
        colnames <- colnames(y)
        for(icol in seq_len(ncol(y))) {
            ycol <- y[,icol]
            if(is.character(ycol))
                # model.matrix in expand.arg would later convert the column to two columns
                # if nbr unique strings is two, which is incorrect, so block that here
                # TODO it would be better to process this rather than issue an error message
                stop0("y is a character variable: \"",
                      ycol[1], "\", \"",
                      if(length(y) >= 2) ycol[2] else "", "\", \"",
                      if(length(y) >= 3) ycol[3] else "", "\", ...")
            else if(is.logical(ycol))
                convert.logical <- TRUE
            else if(is.factor(ycol) && nlevels(ycol) <= 2) { # two-level factor?
                convert.factor <- TRUE
                if(!is.null(colnames) || ncol(y) == 1)
                    colnames[icol] <- levels(ycol)[2]
            }
        }
        nrow <- nrow(y)
        ncol <- ncol(y)
        if(convert.logical) {
            y <- as.numeric(y)      # convert to 0s and 1s
            dim(y) <- c(nrow, ncol)
            colnames(y) <- colnames
        } else if(convert.factor) {
            y <- as.numeric(y) - 1  # minus 1 to convert to 0s and 1s
            dim(y) <- c(nrow, ncol)
            colnames(y) <- colnames
        }
    }
    y
}
