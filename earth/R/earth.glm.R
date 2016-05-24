# earth.glm.R: Generalized Linear Model support for earth

try.something.like <-
    "Try something like earth(y~x, glm=list(family=binomial))"

check.glm.model <- function(g, resp.name) # g is a model created by calling glm()
{
    # TODO following check is pointless? df is only defined for summary.glm?
    df <- if("df" %in% names(g)) g[["df"]] else NULL # avoid partial matching
    if(!is.null(df) && (nsingular <- df[3] - df[1]))
        stop0("earth glm response \"", resp.name, "\": ", nsingular,
              " coefficients not defined because of singularities")

    glm.coef <- coef(g)
    if(length(glm.coef) == 0)
        stop0("earth glm response \"", resp.name, "\": no glm coefficients")
    check.vec(glm.coef, "glm coef")
}
# Check for a common user error: specifying a family argument
# to earth that is not wrapped in glm=list(family=...))

check.no.family.arg.to.earth <- function(...)
{
    dots <- match.call(expand.dots=FALSE)$...
    if(!is.null(dots$fa)) # partial match
        stop0("illegal 'family' argument to earth\n", try.something.like)
}
# This duplicates some tests in binomial in family.R for
# better error reporting in the earth context

check.yarg.for.binomial.glm <- function(yarg, mustart, more.y.columns)
{
    if(ncol(yarg) == 1 && any(yarg < 0 | yarg > 1)) {
        cat("Error: binomial glm model with a vector y:",
            "y values must be between 0 and 1\n")
        if(more.y.columns)
            cat("Possible remedy: pair this column with",
                "the next column using the 'bpairs' arg\n")
        cat0("The first few rows of the y argument to glm are\n")
        print(head(yarg))
        stop0("glm with earth failed, see above messages")
    }
}
# Note that on entry get.glm.arg has already checked the glm argument
# Most args are direct copies of args to earth.fit.

earth.glm <- function(bx, y, weights, na.action, glm,
                      trace, glm.bpairs, resp.name, env)
{
    hack.intercept.only.glm.model <- function(g)
    {
        # Skullduggery for intercept-only glm models.
        # Get the model into a form usable by later functions like predict().
        # We need to remove all references to EarthIntercept else predict()
        # will try to find it and complain because it cannot.

        g$coefficients <- g$coefficients[1]
        g$R            <- g$R[1,1]
        g$qr$qr        <- g$qr$qr[,1,drop=FALSE]
        g$model        <- g$model[,1,drop=FALSE]
        g$x            <- g$x[,1,drop=FALSE]
        g$data         <- NULL

        # yarg ~ EarthIntercept becomes yarg ~ yarg
        # TODO this approach used because won't allow just yarg~
        g$terms[[3]] <- g$terms[[2]]

        # list(yarg, EarthIntercept) ecomes list(yarg)
        attr(g$terms, "variables") <- call("list", quote(yarg))

        #        EarthIntercept
        # yarg                0   becomes an empty matrix

        attr(g$terms, "factors") <- matrix(nrow=0, ncol=0)

        # "EarthIntercept" becomes an empty character vector
        attr(g$terms, "term.labels") <- character(0)

        # list(yarg, EarthIntercept) becomes list(yarg)
        attr(g$terms, "predvars") <- call("list", quote(yarg))

        g
    }
    #--- earth.glm starts here ---
    if(trace >= 2)
        cat("\n")
    ncases <- nrow(bx)
    intercept.only <- ncol(bx) == 1
    if(intercept.only) {
        # glm() requires something on the rhs of the formula.
        # But this is an intercept-only model, so actually nothing on the rhs.
        # To work around that, give glm() the earth intercept, which will have
        # no effect on the glm model but will cause an extra coefficient etc. in
        # the value returned by glm.  We remove that extra data later (in
        # hack.intercept.only.glm.model).
        # Actually the fake intercept does have a small effect on the
        # model: dof is off by one (which also affects vals derived from dof).
        trace1(trace, "earth.glm: intercept-only earth model, faking the glm model\n")
        bx.data.frame <- as.data.frame(bx) # bx has a single column, the earth intercept
        colnames(bx.data.frame) <- "EarthIntercept" # for sanity checking
    } else {
        # default operation: drop intercept with -1
        bx.data.frame <- as.data.frame(eval.parent(bx[, -1, drop=FALSE]))
    }
    # Convert args to form expected by glm().
    # We need to convert glm() args whose default is not NULL.

    control <- glm$control
    if(is.null(control))
       control <- glm.control()
    # FIXED (earth 2.3-5): get control params
    if(!is.null(glm$epsilon))
       control$epsilon <- glm$epsilon
    if(!is.null(glm$maxit))
       control$maxit <- glm$maxit
    if(!is.null(glm$trace))
       control$trace <- glm$trace
    family <- get.glm.family(glm$family, env=env)
    is.binomial <- is.binomial(family)
    stopifnot(is.null(glm$weights))

    # Fit a glm model for each y column.  Except that if there are
    # paired y columns then fit a single glm for each pair.
    # Note that we don't need to look at earth's wp argument here
    # because each glm model is built independently.

    iycol <- 1          # y column index
    imodel <- 1         # model index
    glm.list <- list()  # returned list of glm models

    while(iycol <= ncol(y)) {
        # get yarg, the response for call to glm()
        # it will be a single column or a a paired binomial column

        if(is.null(glm.bpairs))
            yarg <- y[, iycol, drop=FALSE]      # single y column
        else {
            if(!glm.bpairs[iycol])
                stop0(
"unmatched FALSE value in 'bpairs' for y column ", iycol, "\n",
"       Each FALSE in 'bpairs' should be preceded by a TRUE\n",
"       Your bpairs is ", paste.collapse(glm.bpairs))

            if(iycol + 1 <= ncol(y) && !glm.bpairs[iycol+1]) {
                yarg <- y[, c(iycol,iycol+1)]   # paired y columns
                iycol <- iycol + 1
            } else                              # single y column
                yarg <- y[, iycol, drop=FALSE]

        }
        if(is.binomial)
            check.yarg.for.binomial.glm(yarg, glm$mustart, iycol < ncol(y))
        iycol <- iycol + 1
        stopifnot(!is.null(colnames(yarg)))
        if(trace >= 4) {
            print_summary(yarg, "y arg to glm()", trace=2)
            printf("\n")
        }
        # FIXED (earth 2.3-4): removed offset etc. arguments because of
        # difficulties evaluating them later in the correct environment
        # (get.glm.arg has already checked if such args were supplied by the user).
        # TODO if weights are used, glm gives warning "non-integer #successes in a binomial glm"

        g <- glm(yarg ~ ., family=family, data=bx.data.frame,
                weights=weights, na.action=na.action,
                control=control, model=TRUE, trace=(trace>=2),
                method="glm.fit", x=TRUE, y=TRUE, contrasts=NULL)

        if(intercept.only)
            g <- hack.intercept.only.glm.model(g)

        if(trace == 0 && !g$converged) # give a message specific to this response
            cat0("earth glm ", resp.name,
                 ": did not converge after ", g$iter," iterations\n")
        if(trace >= 1) {
            cat0("GLM ", colnames(yarg)[1], ": ")
            print.one.earth.glm(g, digits=getOption("digits"))
        }
        check.glm.model(g, colnames(yarg)[1])
        glm.list[[imodel]] <- g
        imodel <- imodel + 1
    }
    glm.list
}
# process family here instead of in glm() so can give relevant error message

get.glm.family <- function(family, env)
{
    if(is.null(family))
        family <- gaussian
    if(is.character(family))
        family <- get(family, mode="function", envir=env)
    if(is.function(family))
        family <- family()
    if(is.null(family$family))
        stop0("earth: illegal 'family' in 'glm' argument\n",
              try.something.like)
    family
}
# This returns the glm argument but with abbreviated names
# expanded to their full name.  It also checks that the glm argument is
# valid. Called before calling glm().  We want to make sure that the
# user hasn't specified, say, subset as a glm argument. The subset
# should only be specified as an earth argument so the subset is the
# same for earth and glm.
# FIXED (earth 2.3-4): disallow offset etc. arguments because of
# difficulties evaluating them later in the correct environment.

get.glm.arg <- function(glm)    # glm arg is earth's glm arg
{
   # return glm but with abbreviated names expanded to their full name.

    match.glm.arg <- function(glm)
    {
        glm.args <- c("formula", "family", "data", "weights", "subset",
            "na.action", "control", "model", "method", "x", "y",
            "contrasts", "epsilon", "maxit", "trace", "bpairs")

        for(i in seq_along(glm)) {
            j <- pmatch(names(glm)[[i]], glm.args, nomatch=0)
            if(j == 0)
                stop0("earth: '", names(glm)[[i]],
                      "' is not supported in glm argument to earth")
            names(glm)[[i]] <- glm.args[j]
        }
        # expand family argument if it is a string

        if(is.character(glm$family)) {
            family.strings <-
              c("binomial", "gaussian",  "Gamma",  "inverse.gaussian",
                "poisson",  "quasi",  "quasibinomial",  "quasipoisson")

            i <- pmatch(glm$family, family.strings, nomatch=0)
            if(i == 0)
                stop0("earth: illegal family '", glm$family, "' in glm argument\n",
                      try.something.like)
            glm$family <- family.strings[i]
        }
        glm
    }
    #--- get.glm.arg starts here ---

    if(!is.list(glm))
        stop0("earth: 'glm' argument must be a list\n", try.something.like)
    if(length(glm) == 0)
        stop0("earth: 'glm' argument list is empty\n", try.something.like)
    argnames <- names(glm)
    if(length(argnames) == 0)
        stop0("earth: no argument names in 'glm' argument list\n", try.something.like)
    glm <- match.glm.arg(glm)  # expand argument names to their full name
    if(is.null(glm$family))
        stop0("earth: 'glm' argument must have a 'family' parameter\n",
              try.something.like)

    always.true.args <- c("x", "y", "model")

    imatch <- pmatch(always.true.args, argnames)
    imatch <- imatch[!is.na(imatch)]
    if(any(imatch))
        stop0("earth: illegal '", argnames[imatch[1]], "' in 'glm' argument\n",
              "These are always effectively TRUE")

    earths.args <- c("formula", "subset", "weights")
    imatch <- pmatch(earths.args, argnames)
    imatch <- imatch[!is.na(imatch)]
    if(any(imatch)) {
        stop0("earth: illegal '", argnames[imatch[1]], "' in 'glm' argument\n",
              "Use earth's '", argnames[imatch[1]], "' argument instead")
    }
    glm
}
# get.glm.coefs returns a ncoeffs * nresponses matrix

get.glm.coefs <- function(glm.list, nresp, selected.terms, term.names, resp.names)
{
    coefs <- matrix(nrow=length(selected.terms), ncol=nresp)
    col.names <- character(length=nresp)
    for(iresp in seq_len(nresp)) {
        coefs[,iresp] <- glm.list[[iresp]]$coefficients
        col.names[iresp] <- resp.names[iresp]
    }
    colnames(coefs) <- col.names
    rownames(coefs) <- term.names[selected.terms]
    coefs
}
# Return a boolean vector saying which cols in y must be passed
# on to the C earth routine.
# If the user explicitly specified bpairs then we use that bpairs.
# Else we try to figure out bpairs automatically.
# Returns NULL if all y columns should be used
# (and returns NULL if family is not binomial).

get.glm.bpairs <- function(y, glm)
{
    check.no.na.in.mat(y)
    bpairs <- glm$bpairs
    if(!is.null(bpairs)) {              # bpairs provided by user?
        if(ncol(y) == 1)
            stop0("'bpairs' argument is not allowed because y has only one column")
        if(!is.binomial(glm$family))
            stop0("'bpairs' argument is not allowed because the family ",
                  "is not binomial or quasibinomial")
        bpairs <- check.index(bpairs, "bpairs", y, is.col.index=TRUE)
        bpairs <- to.logical(bpairs, NCOL(y))
    } else {
        bpairs <- repl(TRUE, ncol(y))
        if(is.binomial(glm$family)) {
            # If two adjacent columns both have values <0 or >1 then
            # assume that the columns are paired

            if(nrow(y) < 1)
                return(NULL) # later check will issue err msg for too short y
            i <- 1           # column number
            repeat {
               if(i + 1 > ncol(y))
                   break
               if(any(y[,i] < 0 | y[,i] > 1) && any(y[,i+1] < 0 | y[,i+1] > 1)) {
                   bpairs[i+1] <- FALSE # columns are paired, so discard 2nd column
                   i <- i + 1
               }
               i <- i + 1
            }
        }
    }
    if(all(bpairs))
        NULL                                # all y columns used
    else
        bpairs
}
is.binomial <- function(family) # return true if family is binom or quasibinom
{
    (
        (is.character(family) &&            # e.g. "binomial"
            (substr(family, 1, 1) == "b" ||
             substr(family, 1, 6) == "quasib"))
        ||
        (class(family)[1] == "function" &&  # e.g. binomial
            (identical(body(family), body(binomial)) ||
             identical(body(family), body(quasibinomial))))
        ||
        (class(family)[1] == "family" &&    # e.g. binomial()
            (family$family == "binomial" ||
             family$family == "quasibinomial"))
    )
}
is.poisson <- function(family) # return true if family is poisson or quasipoisson
{
    (
        (is.character(family) &&            # e.g. "poisson"
            (substr(family, 1, 1) == "p" ||
             substr(family, 1, 6) == "quasip"))
        ||
        (class(family)[1] == "function" &&  # e.g. poisson
            (identical(body(family), body(poisson)) ||
             identical(body(family), body(quasipoisson))))
        ||
        (class(family)[1] == "family" &&    # e.g. poisson()
            (family$family == "poisson" ||
             family$family == "quasipoisson"))
    )
}
# called from print.summary.earth

print.earth.glm <- function(object, digits, fixed.point)    # object is an earth object
{
    glm.list <- object$glm.list
    nresp <- length(glm.list)

    cat("\nGLM ")
    if(nresp == 1)
        print.one.earth.glm(glm.list[[1]], digits)
    else {                                  # create a matrix and print that
        cat0("(family ", glm.list[[1]]$family$family, ", link ",
             glm.list[[1]]$family$link, ")\n")

        a <- matrix(nrow=nresp, ncol=6)
        colnames(a) <- c("null.deviance", "df", "deviance", "df", "iters", "converged")
        rownames(a) <- colnames(object$fitted.values)
        for(iresp in seq_len(nresp)) {
            g <- glm.list[[iresp]]
            a[iresp,] <- c(g$null.deviance, g$df.null,
                           g$deviance, g$df.residual,
                           g$iter, g$converged)
        }
        if(fixed.point)
            a <- my.fixed.point(a, digits)
        print(a, digits=digits)
    }
}
# Called from print.summary.earth
# g is a glm object
# Most of the following was lifted from print.summary.glm
# but tweaked to include response names (necessary for multiple
# response glm earth models).

print.glm.details <- function(g, nresp, digits, fixed.point, resp.name)
{
    if(nresp > 1)
        prefix <- paste("GLM", resp.name)
    else
        prefix <- paste("GLM")
    sumg <- summary(g)
    cat0(prefix, " deviance residuals:\n")
    if(sumg$df.residual > 5) {
        sumg$deviance.resid <- quantile(sumg$deviance.resid,na.rm=TRUE)
        names(sumg$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
    }
    print.default(sumg$deviance.resid, digits=digits, na.print="", print.gap=2)
    df <- if("df" %in% names(sumg)) sumg[["df"]] else NULL
    cat0("\n", prefix)
    cat0(" coefficients (family ", g$family$family, ", link ", g$family$link, ")\n")
    if(!is.null(df) && (nsingular <- df[3] - df[1]))
        cat0(nsingular,  # should never happen for earth glm
             " coefficients not defined because of singularities\n")
    aliased <- is.na(coef(g))
    stopifnot(length(aliased) > 0, all(!aliased)) # already checked in check.glm.model
    coefs <- sumg$coefficients
    rownames(coefs) <- spaceout(rownames(coefs))
    # TODO can't use fixed.point here, would like to
    printCoefmat(coefs, digits=digits, signif.stars=FALSE, na.print="NA")
    if(sumg$dispersion != 1) # only show dispersion if it is not 1
        cat0("\n", prefix, " dispersion parameter for ", sumg$family$family,
             " family taken to be ", format(sumg$dispersion), "\n")
    cat("\n")
    NULL
}
print.one.earth.glm <- function(g, digits)
{
    cat0("null.deviance ",  format(g$null.deviance, digits=digits),
         " (",             g$df.null, " dof)",
         "   deviance ",   format(g$deviance, digits=digits),
         " (",             g$df.residual, " dof)",
         "   iters ",      g$iter,
         if(!g$converged) " did not converge" else "",
         "\n")
}
