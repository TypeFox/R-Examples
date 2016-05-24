##' @export
##' @method ascii survdiff
ascii.survdiff <- function (x, include.rownames = TRUE, include.colnames = TRUE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = TRUE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
    # From print.survdiff
    if (length(x$n) == 1) {
        z <- sign(x$exp - x$obs) * sqrt(x$chisq)
        temp <- c(x$obs, x$exp, z, signif(1 - pchisq(x$chisq,
            1), digits))
        names(temp) <- c("Observed", "Expected", "Z", "p")
        temp <- t(temp)
        include.rownames = FALSE
    } else {
        if (is.matrix(x$obs)) {
            otmp <- apply(x$obs, 1, sum)
            etmp <- apply(x$exp, 1, sum)
        } else {
            otmp <- x$obs
            etmp <- x$exp
        }
        df <- c((sum(1 * (etmp > 0))) - 1, rep(NA, length(x$n) - 1))
        p <- c(1 - pchisq(x$chisq, df[!is.na(df)]), rep(NA, length(x$n) - 1))
        temp <- cbind(x$n, otmp, etmp, ((otmp - etmp)^2)/etmp,
            ((otmp - etmp)^2)/diag(x$var), df, p)
        dimnames(temp) <- list(names(x$n), c("N", "Observed",
            "Expected", "(O-E)^2/E", "(O-E)^2/V", "df", "p"))
    }

  temp <- as.data.frame(temp, checknames = FALSE)

  obj <- asciiTable$new(x = temp, include.rownames = include.rownames,
      include.colnames = include.colnames, rownames = rownames, colnames = colnames,
      format = format, digits = digits, decimal.mark = decimal.mark, na.print = na.print,
      caption = caption, caption.level = caption.level, width = width, frame = frame,
      grid = grid, valign = valign, header = header, footer = footer, align = align,
      col.width = col.width, style = style,
      tgroup = tgroup, n.tgroup = n.tgroup, talign = talign,
      tvalign = tvalign, tstyle = tstyle,
      bgroup = bgroup, n.bgroup = n.bgroup, balign = balign,
      bvalign = bvalign, bstyle = bstyle,
      lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign,
      lvalign = lvalign, lstyle = lstyle,
      rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign,
      rvalign = rvalign, rstyle = rstyle)
  return(obj)
}

##' @param scale A numeric value to rescale the survival time, e.g., if the
##'   input data to survfit were in days, \code{scale=365} would scale the
##'   printout to years (see \code{print.survfit()} in package
##'   \code{survival}).
##' @param print.rmean Option for computation and display of the restricted mean (see \code{print.survfit()} in package \code{survival}).
##' @param rmean Option for computation and display of the restricted mean (see \code{print.survfit()} in package \code{survival}).
##' @export
##' @method ascii survfit
##' @rdname ascii
ascii.survfit <- function (x, scale = 1, print.rmean = getOption("survfit.print.rmean"), rmean = getOption("survfit.rmean"), include.rownames = TRUE, include.colnames = TRUE, header = TRUE, ...) {
    omit <- x$na.action
    na <- NULL
    if (length(omit))
         na <- ascii(list(naprint(omit)), list.type = "none")
    if (!missing(print.rmean) && is.logical(print.rmean) && missing(rmean)) {
        if (print.rmean)
            rmean <- "common"
        else rmean <- "none"
    }
    if (is.null(rmean)) {
        if (is.logical(print.rmean)) {
            if (print.rmean)
                rmean <- "common"
            else rmean <- "none"
        }
        else rmean <- "none"
    }
    if (is.numeric(rmean)) {
        if (is.null(x$start.time)) {
            if (rmean < min(x$time))
                stop("Truncation point for the mean is < smallest survival")
        }
        else if (rmean < x$start.time)
            stop("Truncation point for the mean is < smallest survival")
    }
    else {
        rmean <- match.arg(rmean, c("none", "common", "individual"))
        if (length(rmean) == 0)
            stop("Invalid value for rmean option")
    }
    temp <- survival:::survmean(x, scale = scale, rmean)
    mat <- ascii(temp$matrix, include.rownames = include.rownames, include.colnames = include.colnames, header = header, ...)

    restrm <- NULL
    if (rmean != "none") {
        if (rmean == "individual")
            restrm <- ascii(list("* restricted mean with variable upper limit"))
        else restrm <- ascii(list(paste("* restricted mean with upper limit = ",
            format(temp$end.time[1]))))
    }
    obj <- asciiMixed$new(args = list(na, mat, restrm))
    obj
}

##' @export
##' @method ascii summary.survfit
ascii.summary.survfit <- function (x, include.colnames = TRUE, header = TRUE, digits = c(0, 0, 0, 3, 3, 3, 3), ...) {
  omit <- x$na.action
  na <- NULL
  if (length(omit))
    na <- ascii(list(naprint(omit)), list.type = "none")
  if (x$type == "right" || is.null(x$n.entered)) {
    mat <- cbind(x$time, x$n.risk, x$n.event, x$surv)
    cnames <- c("time", "n.risk", "n.event")
  } else if (x$type == "counting") {
    mat <- cbind(x$time, x$n.risk, x$n.event, x$n.enter,
                 x$n.censor, x$surv)
    cnames <- c("time", "n.risk", "n.event", "entered", "censored")
  }
  if (is.matrix(x$surv)) {
    ncurve <- ncol(x$surv)
  } else ncurve <- 1
  if (ncurve == 1) {
    cnames <- c(cnames, "survival")
    if (!is.null(x$std.err)) {
      if (is.null(x$lower)) {
        mat <- cbind(mat, x$std.err)
        cnames <- c(cnames, "std.err")
      } else {
        mat <- cbind(mat, x$std.err, x$lower, x$upper)
        cnames <- c(cnames, "std.err", paste("lower ",
                                             x$conf.int * 100, "% CI", sep = ""), paste("upper ",
                                                                                        x$conf.int * 100, "% CI", sep = ""))
      }
    }
  } else cnames <- c(cnames, paste("survival", seq(ncurve), sep = ""))
  if (!is.null(x$start.time)) {
    mat.keep <- mat[, 1] >= x$start.time
    mat <- mat[mat.keep, , drop = FALSE]
    if (is.null(dim(mat)))
      stop(paste("No information available using start.time =",
                 x$start.time, "."))
  }
  if (!is.matrix(mat))
    mat <- matrix(mat, nrow = 1)
  if (!is.null(mat)) {
    dimnames(mat) <- list(NULL, cnames)
    if (is.null(x$strata)) {
      res <- ascii(mat, include.colnames = include.colnames, header = header, digits = digits, ...)
    } else {
      strata <- x$strata
      if (!is.null(x$start.time))
        strata <- strata[mat.keep]
      res <- NULL
      for (i in levels(strata)) {
        who <- (strata == i)
        res <- asciiMixed$new(args = list(res, ascii(mat[who, ], caption = i, include.colnames = include.colnames, header = header, digits = digits, ...)))
      }
    }
  } else stop("There are no events to print.  Please use the option ",
            "censored=TRUE with the summary function to see the censored ",
            "observations.")
  obj <- asciiMixed$new(args = list(na, res))
  return(obj)
}

# based on xtable package

##' @export
##' @method ascii coxph
ascii.coxph <- function (x, include.rownames = TRUE, include.colnames = TRUE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = TRUE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){

    cox <- x
    beta <- cox$coef
    se <- sqrt(diag(cox$var))
    tmp <- cbind(beta, exp(beta), se, beta/se, 1 - pchisq((beta/se)^2, 1))
    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", "se(coef)", "z", "p"))

    obj <- asciiTable$new(x = as.data.frame(tmp), include.rownames = include.rownames,
      include.colnames = include.colnames, rownames = rownames, colnames = colnames,
      format = format, digits = digits, decimal.mark = decimal.mark, na.print = na.print,
      caption = caption, caption.level = caption.level, width = width, frame = frame,
      grid = grid, valign = valign, header = header, footer = footer, align = align,
      col.width = col.width, style = style,
      tgroup = tgroup, n.tgroup = n.tgroup, talign = talign,
      tvalign = tvalign, tstyle = tstyle,
      bgroup = bgroup, n.bgroup = n.bgroup, balign = balign,
      bvalign = bvalign, bstyle = bstyle,
      lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign,
      lvalign = lvalign, lstyle = lstyle,
      rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign,
      rvalign = rvalign, rstyle = rstyle)
  return(obj)
}
