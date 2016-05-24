# plot.glmnetx.R:
#
# Same as plot.glmnet except
# (i)   new type "rlambda"  ("r" for reverse) with new args s and nresponse
# (ii)  use dots package so can change xlab, ylab, col, etc. using dot args
# (iii) label is now the nbr of labels to display, special value TRUE means all
# (iv)  new argument grid.col to optionally add a grid
# (v)   simplicity of original plot.glmnet code is gone :(
#
# This code is based on glmnet version 2.0-1.

plot.glmnetx <- function(x, xvar=c("norm","lambda","rlambda","dev"), label=20,
                         grid.col=NA, s=NA, nresponse=1, norm=NULL, ...)
{
    object <- x
    beta <- object$beta
    if(is.list(beta)) { # multiple response model?
        check.index(nresponse, "nresponse", beta)
        beta <- beta[[nresponse]]
    }
    ibeta <- nonzeroCoef(beta) # ibeta is a vector of coefficient indices
    if(length(ibeta) == 0) {
        plot(0:1, 0:1, col=0)
        legend("topleft", legend="all coefficients are zero", bty="n")
        return()
    }
    # following was in original plot.glmnet code but seems unnecessary?
    # if(length(ibeta) == 1) {
    #     plot(0:1, 0:1, col=0)
    #     legend("topleft", legend="only one coefficient is nonzero", bty="n")
    #     return()
    # }
    beta <- as.matrix(beta[ibeta, , drop=FALSE])
    xvar <- match.arg(xvar)
    switch(xvar,
    "norm"= {
        if(inherits(object, "multnet") || inherits(object, "mrelnet")) {
            # we don't (yet) precalc norm or support type.coef, so have to stop here
            stop0("xvar=\"norm\" is not supported by plotres for multiple responses")
        }
        stopifnot(is.null(norm))
        x <- if(is.null(norm)) apply(abs(beta), 2, sum) else norm
        xlim <- c(min(x), max(x))
        xlab <- "L1 Norm"
        approx.f <- 1
    },
    "lambda"= {
        x <- log(object$lambda)
        xlim <- c(min(x), max(x))
        xlab <- "Log Lambda"
        approx.f <- 0
    },
    "rlambda"= {
        x <- log(object$lambda)
        xlim <- c(max(x), min(x)) # backwards
        xlab <- "Log Lambda"
        approx.f <- 0
    },
    "dev"= {
        x <- object$dev.ratio
        xlim <- c(min(x), max(x))
        xlab <- "Fraction Deviance Explained"
        approx.f <- 1
    })

    # named index of varnames to be printed on right of plot, NULL if none
    iname <- get.iname(beta, ibeta, label)

    opar <- par("mar", "cex.axis", "cex.lab")
    on.exit(par(mar=opar$mar, cex.axis=opar$cex.axis, cex.lab=opar$cex.lab))
    mar4 <- opar$mar[4] # right hand margin
    if(length(iname)) {
        cex.names <- min(1, max(.5, 2 / sqrt(length(iname)))) # seems reasonable
        # ensure right margin is big enough for the varnames
        # can't use strwidth because no plot yet, so just estimate
        mar4 <- max(opar$mar[4] + 1,
                    .75 * cex.names * par("cex") * max(nchar(names(iname))))
    }
    # set mar[3] for top axis label and mar[4] for right hand labels
    par(mar=c(opar$mar[1:2], max(opar$mar[3], 3), mar4))
    par(cex.axis=.8)
    # line colors, max is 6 in default to avoid yellow, rep_len to recycle
    col <- rep_len(dot("col", DEF=1:6, ...), length(ibeta))
    ylab <- "Coefficients"
    if(is.list(object$beta)) # multiple response model?
        ylab <- paste0(ylab, ": Response ", rownames(object$dfmat)[nresponse])

    # any arg prefixed with def. can be overridden by a user-specifed arg in dots
    # main="" because we will later manually add a top axis instead of main
    call.plot(graphics::matplot, force.x=x, force.y=t(beta), force.main="",
        force.col=col, def.xlim=xlim, def.xlab=xlab, def.ylab=ylab,
        def.lty=1, def.lwd=1, def.type="l", ...)

    abline(h=0, col="gray", lty=3) # zero axis line
    maybe.grid(x=x, beta=beta, grid.col=grid.col, col1=col, ...)
    if(xvar == "rlambda") {
        # args are named below to prevent potential clash with argnames in dots
        annotate.rlambda(lambda=object$lambda, x=x, beta=beta, s=s,
                         grid.col=grid.col, col1=col, ...)
        main <- "Lambda"
    } else {
        top.axis(object, x, nresponse, approx.f)
        main <- "Degrees of Freedom"
    }
    mtext(dot("main", DEF=main, ...), side=3,
          line=1.8, cex=par("cex") * par("cex.lab"))
    if(length(iname))
        right.labs(beta, iname, cex.names, col)
}
get.iname <- function(beta, ibeta, label)
{
    iname <- NULL
    check.integer.scalar(label, min=0, logical.ok=TRUE, na.ok=TRUE)
    if(!is.na(label) && label) { # allow label=NA, treat as FALSE
        names <- if(is.null(rownames(beta))) paste(ibeta)
                 else                        rownames(beta)
        names[!nzchar(names)] <- paste(ibeta)[!nzchar(names)]
        iname <- order(abs(beta[, ncol(beta)]), decreasing=TRUE)
        if(label != 1 && length(iname)>label) # label=1 or TRUE is special meaning all
            iname <- iname[1:label]
        names(iname) <- abbreviate(names[iname], minlength=8)
    }
    iname # named index of varnames to be printed, NULL if none
}
maybe.grid <- function(x, beta, grid.col, col1, ...)
{
    if(is.specified(grid.col[1])) {
        grid(col=grid.col, lty=1)
        # replot over the grid (using add=TRUE)
        # col1 is not called col, to prevent a clash with col which may be in dots
        call.plot(graphics::matplot, force.x=x, force.y=t(beta), force.main="",
            force.add=TRUE, force.col=col1,
            def.lty=1, def.lwd=1, def.type="l", ...)
    }
}
right.labs <- function(beta, iname, cex.names, col) # varnames on right of plot
{
    usr <- par("usr")
    text(x=usr[2] + .01 * (usr[2] - usr[1]),
         y=TeachingDemos::spread.labs(beta[iname, ncol(beta)],
                                      mindiff=1.2 * cex.names * strheight("X")),
         labels=names(iname), cex=cex.names, col=col[iname], adj=0, xpd=NA)
}
top.axis <- function(object, x, nresponse, approx.f)
{
    at <- pretty(x)
    # use is.list(object$beta) to determine if multiple response model
    df <- if(is.list(object$beta)) object$dfmat[nresponse,] else object$df
    # compute df by interpolating to df at next smaller lambda
    # thanks to Yunyang Qian
    prettydf <- approx(x=x, y=df, xout=at,
                       rule=2, method="constant", f=approx.f)$y
    axis(3, at=at, labels=prettydf)
}
# Draw the top axis of an rlambda plot.  Also draw a labeled vertical
# line at lambda=s, if s isn't NA.  Dot arguments prefixed with "s". can
# be used set the annotation attributes e.g. s.col=NA or s.col=0 for no vert line.
annotate.rlambda <- function(lambda, x, beta, s, grid.col, col1, ...)
{
    check.numeric.scalar(s, null.ok=TRUE, na.ok=TRUE)
    s.col <- dot("s.col", DEF=1, ...)
    add.s.line <- !is.null(s) && !is.na(s) && is.specified(s.col)

    # top axis
    at <- pretty(x)
    labs <- signif(exp(at), digits=2)
    # hack: delete confusing rightmost lab (if any) with a value greater
    # than s but drawn to the right of the vertical line at s
    if(add.s.line && s <= labs[1])
        labs[1] <- ""
    axis(3, at=at, labels=labs)

    if(add.s.line)
        add.s.line(lambda=lambda, x=x, beta=beta, s=s,
                   grid.col=grid.col, col1=col1, s.col=s.col, ...)
}
add.s.line <- function(lambda, x, beta, s, grid.col, col1, s.col, ...)
{
    line.col <- "gray"
    line.lty <- 1
    if(is.specified(grid.col)) {
            line.col <- 1
            line.lty <- 3
    }
    log.s <- log(max(lambda[length(lambda)], s))

    abline(v=log.s, col=line.col, lty=line.lty) # vertical line at s

    # replot over the vertical line (using add=TRUE)
    call.plot(graphics::matplot, force.x=x, force.y=t(beta), force.main="",
              force.add=TRUE, force.col=col1,
              def.lty=1, def.lwd=1, def.type="l", ...)

    # add s label on vertical line
    # to minimize overplotting, y coord of label is biggest gap between matplot lines
    usr <- par("usr") # xmin, xmax, ymin, ymax
    col.index <- which.min(abs(lambda-s)) # lambda column corresponding to s
    y <- sort(c(usr[3], beta[, col.index], usr[4])) # include plot edges, and sort
    which <- which.max(diff(y))
    call.plot(text.on.white, PREFIX="s.",
            force.x=log.s, force.y=(y[which]+y[which+1]) / 2,
            force.label= # gsub below drops leading and trailing zeros for compactness
                if(s == 0) "s=0"
                else        paste0("s=", gsub("^0|0$|\\.0*$", "", signif(s,2))),
            force.col=s.col, def.cex=.8, def.srt=90, def.xpd=NA, ...)
}
# return NULL or an integer vector
# reproduced here so don't have to import glmnet
nonzeroCoef = function (beta, bystep = FALSE)
{
### bystep = FALSE means which variables were ever nonzero
### bystep = TRUE means which variables are nonzero for each step
  nr=nrow(beta)
  if (nr == 1) {#degenerate case
    if (bystep)
      apply(beta, 2, function(x) if (abs(x) > 0)
            1
      else NULL)
    else {
      if (any(abs(beta) > 0))
        1
      else NULL
    }
  }
  else {
    beta=abs(beta)>0 # this is sparse
    which=seq(nr)
    ones=rep(1,ncol(beta))
    nz=as.vector((beta%*%ones)>0)
    which=which[nz]
    if (bystep) {
      if(length(which)>0){
        beta=as.matrix(beta[which,,drop=FALSE])
        nzel = function(x, which) if (any(x))
          which[x]
        else NULL
        which=apply(beta, 2, nzel, which)
        if(!is.list(which))which=data.frame(which)# apply can return a matrix!!
        which
      }
      else{
        dn=dimnames(beta)[[2]]
        which=vector("list",length(dn))
        names(which)=dn
        which
      }

    }
    else which
  }
}
