if.R(s={},
     r={
as.multicomp <- function (x, ...)
  UseMethod("as.multicomp")

as.glht <- function (x, ...)
  UseMethod("as.glht")

as.multicomp.glht <-
  function(x,       ## glht object
           focus=x$focus,
           ylabel=deparse(terms(x$model)[[2]]),
           means=model.tables(x$model, type="means", cterm=focus)$tables[[focus]],
           height=rev(1:nrow(x$linfct)),
           lmat=t(x$linfct),
           lmat.rows=lmatRows(x, focus),
           lmat.scale.abs2=TRUE,
           estimate.sign=1,
           order.contrasts=TRUE,
           contrasts.none=FALSE,
           level=0.95,
           calpha=NULL,
           method=x$type,
           df,
           vcov.,
           ...
           ) {
    focus.tmp <- focus ## force evaluation

    dimnames(x$linfct)[[1]] <- gsub(" ", "", dimnames(x$linfct)[[1]]) ## remove blanks
    if (dimnames(x$linfct)[[2]][1] == "") dimnames(x$linfct)[[2]][1] <- "(Intercept)"
    if (!missing(vcov.)) x$vcov <- vcov.(x$model)

    confint.x <-
      if (is.null(calpha))
        confint(x, level=level)  ##, ...)
      else
        confint(x, level=level, calpha=calpha)  ##, ...)

    result <- list(table=cbind(
                     estimate=confint.x$confint[,"Estimate"], #
                     stderr=0, ## placeholder
                     lower=confint.x$confint[,"lwr"], #
                     upper=confint.x$confint[,"upr"]  #
                     ),
                   alpha=1-attr(confint.x$confint,"conf.level"), #
                   error.type=NA,
                   method=method,                       #
                   crit.point=attr(confint.x$confint,"calpha"), #
                   Srank=NULL,
                   simsize=NULL,
                   ylabel=ylabel,
                   call=sys.call(),         #
                   lmcall=x$model,          #
                   focus=focus,             #
                   ## lmat=lmat,         #
                   glht=x
                   )
    result$table[,"stderr"] <-  ## correct stderr for the contrast
      (result$table[,"upper"] - result$table[,"lower"]) / (2*result$crit.point)
    if (is.null(dimnames(result$table)[[1]]))
      dimnames(result$table)[[1]] <- dimnames(confint.x$confint)[[1]]
    tmp <- lmat[lmat.rows, , drop=FALSE]
    if (contrasts.none) {
      first.row <- tmp[1,,drop=FALSE]
      first.row[] <- 0
      first.row[1,1] <- 1
    }
    else
      first.row <- -apply(tmp, 2, sum)
    lmat.subscript <- rbind(first=first.row, tmp)
    lmat.factor <- lmat.subscript

    lmat.factor <- sweep(lmat.factor, 2, apply(abs(lmat.factor), 2, sum)/2, "/")
    if (length(means) != nrow(lmat.factor))
      stop("Please specify lmat.rows with mmc on a design with more than one factor.")
##    result$height <- (means %*% abs(lmat.factor))[1,]
    result$height <- if (is.matrix(height)) height[1,] else height

    result$lmat <-
      if (lmat.scale.abs2 && !contrasts.none)
        sweep(lmat, 2, apply(abs(lmat.subscript), 2, sum)/2, "/")
      else
        lmat
    if (order.contrasts)
      result <- multicomp.order(result)

    result$bounds <- switch(x$alternative,
                            "two.sided"="both",
                            "greater"="lower",
                            "less"="upper")
    class(result) <- c("multicomp.hh", "multicomp")

    result <- multicomp.reverse(result, estimate.sign)
    result$glht$linfct <- t(result$lmat)
    result
  }


as.glht.multicomp <- function(x, ...) x$glht

glht.mmc <- function(...)
  .Defunct("mmc", package="HH")
mmc <- function (model, ...)
  UseMethod("mmc")

## mmc.lm <- function (model,
##            linfct=NULL,
##            focus=
##            if (is.null(linfct))
##            {
##              if (length(model$contrasts)==1) names(model$contrasts)
##              else stop("focus or linfct must be specified.")
##            }
##            else
##            {
##              if (is.null(names(linfct)))
##                stop("focus must be specified.")
##              else names(linfct)
##            },
##            focus.lmat,
##            ylabel=deparse(terms(model)[[2]]),
##            lmat=if (missing(focus.lmat)) {
##              t(linfct)
##            } else {
##              lmatContrast(t(none.glht$linfct), focus.lmat)
##              },
##            lmat.rows=lmatRows(model, focus),

##            lmat.scale.abs2=TRUE,
##            estimate.sign=1,
##            order.contrasts=TRUE,
##            level=.95,
##            calpha=NULL,
##            alternative = c("two.sided", "less", "greater"),
##            ...
##            )
##   NextMethod("glht.mmc")

## glht.mmc.default ## old name
glht.mmc.default <- function(...)
  .Defunct("mmc.default", package="HH")
mmc.default <-     ## this works for model inherits from "lm"
  function(model,  ## It needs work for lme objects
           linfct=NULL,
           focus=
           if (is.null(linfct))
           {
             if (length(model$contrasts)==1) names(model$contrasts)
             else stop("focus or linfct must be specified.")
           }
           else
           {
             if (is.null(names(linfct)))
               stop("focus must be specified.")
             else names(linfct)
           },
           focus.lmat,
           ylabel=deparse(terms(model)[[2]]),
           lmat=if (missing(focus.lmat)) {
             t(linfct)
           } else {
             lmatContrast(t(none.glht$linfct), focus.lmat)
             },
           lmat.rows=lmatRows(model, focus),

           lmat.scale.abs2=TRUE,
           estimate.sign=1,
           order.contrasts=TRUE,
           level=.95,
           calpha=NULL,
           alternative = c("two.sided", "less", "greater"),
           ...
           ) {

    mmm.data <- model$model
    if (inherits(model, "lme")) mmm.data <- model$data
    factors <- sapply(mmm.data, inherits, "factor")
    is.contr.treatment <- function(x) {
      cx <- contrasts(x)
      tx <- contr.treatment(nrow(cx))
      all(cx==tx)
    }
    contrasts.are.treatment <- sapply(mmm.data[, factors, drop=FALSE], is.contr.treatment)
    if (!all(contrasts.are.treatment))
      stop("mmc requires an aov in which ALL factors use treatment contrasts.")

    result <- list(mca=NULL)

##  none.glht <- glht(model, linfct=mcp(focus.value="Means"))

    if (TRUE)
      {
        if (length(focus) > 1) stop("mmc requires no more than one focus factor.")
        focus.linfct <-
          ## multcomp:::meanslinfct(model, focus, formula=terms(model),
          multcomp.meanslinfct(model, focus, formula.in=terms(model),
                                 contrasts.arg=model$contrasts)
        none.glht <- glht(model, linfct=focus.linfct,
                          alternative=alternative, ...)
        if (is.null(none.glht$focus)) none.glht$focus <- focus
      }
    else
      {
        mcp.args <- list("Means")
        names(mcp.args) <- focus
        none.glht <- glht(model, linfct=do.call("mcp", mcp.args),
                          alternative=alternative)
      }

    means <-
      if (is.null(calpha)) {
        confint(none.glht, calpha=1.96)$confint[,"Estimate"] ## fake 1.96 ## , ...
      }
      else
        confint(none.glht, calpha=calpha)$confint[,"Estimate"] ## , ...


    if (is.null(linfct)) {
      linfct.focus <- mcalinfct(model, focus, linfct.Means=none.glht$linfct)
      method="Tukey"
    }
    else {
      if (match(focus, names(linfct), 0) != 0) {
        method <- attr(linfct[[focus]], "type")
        linfct.focus <-
##           if (is.matrix(linfct[[focus]]) &&
##               ncol(linfct[[focus]]) == length(coef(model)))
##             linfct[[focus]]
##           else
            linfct
      }
      else {
        linfct.focus <- linfct
        method <- NULL
      }
    }
## recover()
    mca.glht <- glht(model, linfct=linfct.focus,
                     alternative=alternative, ...)
    if (!is.null(method)) mca.glht$type <- method

    height.mca <-
      if (is.null(method) || method=="Tukey")
        means %*% abs(t(contrMat(table(mmm.data[[focus]]), "Tukey")))
      else
        means %*% abs(t(linfct.focus[[focus]])) ## fixme, this works for Dunnett

    result$mca <- as.multicomp(mca.glht, focus=focus, ylabel=ylabel,
                               means=means,
                               height=height.mca,
                               lmat.rows=lmat.rows,
                               lmat.scale.abs2=lmat.scale.abs2,
                               estimate.sign=estimate.sign,
                               order.contrasts=order.contrasts,
                               calpha=calpha,
                               level=level, ...)

    result$none <- as.multicomp(none.glht, focus=focus, ylabel=ylabel,
                                means=means,
                                height=means*2,
                                lmat=t(none.glht$linfct), lmat.rows=lmat.rows,
                                contrasts.none=TRUE,
                                estimate.sign=0,  ## observed means, contrasts: no reversal!
                                order.contrasts=order.contrasts,
                                level=1-result$mca$alpha,
                                calpha=result$mca$crit.point,
                                method=result$mca$method, ...)
    ## Workaround for incorrect handling of vcov for group means
    ## in multcomp:::lme in R.
    if (inherits(result$none$glht$model, "lme"))
        result$none$table[, c("stderr","lower","upper")] <- NA
    ## End workaround
    if (!missing(lmat) || !missing(focus.lmat)) {
      if (lmat.scale.abs2) {
        tmp <- lmat[lmat.rows, , drop=FALSE]
        first.row <- -apply(tmp, 2, sum)
        lmat.subscript <- rbind(first=first.row, tmp)
        lmat <- sweep(lmat, 2, apply(abs(lmat.subscript), 2, sum)/2, "/")
      }
      lmat.glht <- glht(model, linfct=t(lmat),
                        alternative=alternative, ...)
      if (is.null(lmat.glht$focus)) lmat.glht$focus <- focus
      if (missing(focus.lmat)) stop("'focus.lmat' is missing.")
      result$lmat <-
        as.multicomp(lmat.glht, focus=focus, ylabel=ylabel,
                     means=means,
                     height=means %*%
                     abs(sweep(focus.lmat, 2,
                               apply(abs(focus.lmat), 2, sum)/2, "/"))[names(means),],
                     lmat=lmat, lmat.rows=lmat.rows,
                     estimate.sign=estimate.sign,
                     order.contrasts=order.contrasts,
                     level=1-result$mca$alpha,
                     calpha=result$mca$crit.point,
                     method=result$mca$method, ...)
    }

    class(result) <- "mmc.multicomp"
    result
  }
## assignInNamespace("mmc.default", mmc.default, "HH")


mmc.glht <- function(model, ...) {
  ##  do.call("mmc", c(list(model=model$model), list(...)))
  NextMethod("mmc", model$model)
}

print.glht.mmc.multicomp <- function(...)
  .Defunct("print.mmc.multicomp", package="HH")
## ## prints glht components of mmc.multicomp object
## print.glht.mmc.multicomp <- function (x, ...) {
##   cat(paste("Fit:", deparse(x$mca$glht$model$call, width.cutoff=500), "\n"))
##   cat("Focus =", x$mca$focus, "\n")
##   cat("Estimated Quantile =", x$mca$crit.point, "\n")
##   cat(round((1-x$mca$alpha)*100), "% family-wise confidence level\n", sep="")
##   tmp <- list(mca = x$mca$glht, none = x$none$glht)
##   if (is.null(tmp$none)) tmp$none <- x$none$table
##   if (!is.null(x$lmat))
##     tmp$lmat <- x$lmat$glht
##   print(tmp)
##   invisible(x)
## }


## prints table and height components of multicomp object
print.multicomp <- function (x, ...) {
  if (inherits(x, "TukeyHSD")) ## protect against "multicomp" class in stats:::TukeyHSD output
    NextMethod("print")
  else
    print(cbind(x$table, height=x$height/2))
  invisible(x)
}

## print.multicomp.hh is in print.multicomp.hh.s


## prints table and height components of each multicomp object in a mmc object
print.mmc.multicomp <- function (x, ..., width.cutoff=options()$width-5) {
  cat(paste(x$mca$method, "contrasts\n"))
  mmc.call <- deparse(x$mca$glht$model$call, width.cutoff=width.cutoff)
  mmc.call[1] <- paste("Fit:", mmc.call[1], "\n")
  if (length(mmc.call) > 1)
    mmc.call[-1] <- paste("    ", mmc.call[-1], "\n")
  cat(mmc.call)
  ## cat(paste("Fit:", deparse(x$mca$glht$model$call,
  ##                           width.cutoff=width.cutoff), "\n"))
  cat("Estimated Quantile =", x$mca$crit.point, "\n")
  cat(round((1-x$mca$alpha)*100), "% family-wise confidence level\n", sep="")
  cat("$mca\n")
  print(x$mca)
  cat("$none\n")
  print(x$none)
  if (!is.null(x$lmat)) {
    cat("$lmat\n")
   print(x$lmat)
  }
  invisible(x)
}

plot.multicomp <- function (x, ...) {
  n.contrasts <- dim(x$lmat)[2]
  plot(confint(as.glht(x)), ylim=c(.5, n.contrasts+.5), ...)
}

plot.multicomp.adjusted <- function (x, ...) {
  ## We used the standard plot function
  ##    multcomp:::plot.confint.glht
  ## with possibly different confidence bounds
  n.contrasts <- dim(x$lmat)[2]
  x.confint <- confint(as.glht(x))
  x.adjusted <- x.confint
  x.adjusted$confint[,c("lwr","upr")] <-
    x$table[, c("lower","upper"), drop=FALSE]
  plot(x.adjusted, ylim=c(.5, n.contrasts+.5), ...)
}
## assignInNamespace("plot.multicomp.adjusted", plot.multicomp.adjusted, "HH")


## plot.multicomp.hh is in file plot.multicomp.R
})


## source("~/HH-R.package/HH/R/glht.mmc.R")

## c.mmc <- glht.mmc(catalystm1.aov, linfct = mcp(catalyst = "Tukey"))




##            focus.columns=
##            match(focus, attr(terms(model), "term.labels")) ==
##            attr(model$qr$qr, "assign"),
