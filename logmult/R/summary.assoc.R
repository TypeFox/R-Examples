assocTable <- function(object, ass, ...) {
  if(ass$covtype != "none" && length(ass$covmat) > 0)
      errs <- se(ass)
  else
      errs <- list(phi=NA, row=NA, col=NA)

  paste2 <- function(x, y) paste(y, x, sep=":")

  if(nrow(ass$phi) > 1) {
      phi <- cbind("Normalized"=c(ass$phi), "Adjusted"=NA,
                   "Std. error"=c(errs$phi), "Pr(>|z|)"=c(2 * pnorm(-abs(ass$phi/errs$phi))))

      row <- cbind("Normalized"=c(ass$row), "Adjusted"=NA,
                   "Std. error"=c(errs$row), "Pr(>|z|)"=c(2 * pnorm(-abs(ass$row/errs$row))))

      col <- cbind("Normalized"=c(ass$col), "Adjusted"=NA,
                   "Std. error"=c(errs$col), "Pr(>|z|)"=c(2 * pnorm(-abs(ass$col/errs$col))))

      rownames(phi) <- outer(paste(names(dimnames(object$data))[3], rownames(ass$phi), sep=""),
                             colnames(ass$phi), paste2)
  }
  else {
      phi <- cbind("Normalized"=c(ass$phi), "Adjusted"=NA,
                   "Std. error"=c(errs$phi), "Pr(>|z|)"=c(2 * pnorm(-abs(ass$phi/errs$phi))))

      row <- cbind("Normalized"=c(ass$row), "Adjusted"=c(sweep(ass$row, 2, sqrt(abs(ass$phi)), "*")),
                   "Std. error"=c(errs$row), "Pr(>|z|)"=c(2 * pnorm(-abs(ass$row/errs$row))))

      col <- cbind("Normalized"=c(ass$col), "Adjusted"=c(sweep(ass$col, 2, sqrt(abs(ass$phi)), "*")),
                   "Std. error"=c(errs$col), "Pr(>|z|)"=c(2 * pnorm(-abs(ass$col/errs$col))))

      rownames(phi) <- colnames(ass$phi)
  }


  if(inherits(ass, "assoc.symm")) {
      varnames <- paste(names(dimnames(object$data))[1:2], collapse="|")

      if(dim(ass$row)[3] > 1) {
          rownames(row) <- outer(paste(varnames, rownames(ass$row), sep=""),
                                 outer(colnames(ass$row),
                                       paste(names(dimnames(object$data))[3], dimnames(ass$col)[[3]], sep=""),
                                       paste2),
                                 paste2)
      }
      else {
          rownames(row) <- outer(paste(varnames, rownames(ass$row), sep=""),
                                 colnames(ass$row), paste2)
      }

      rbind(phi, row)
  }
  else {
      if(dim(ass$row)[3] > 1) {
          rownames(row) <- outer(paste(names(dimnames(object$data))[1], rownames(ass$row), sep=""),
                                 outer(colnames(ass$row),
                                       paste(names(dimnames(object$data))[3], dimnames(ass$col)[[3]], sep=""),
                                       paste2),
                                 paste2)
      }
      else {
          rownames(row) <- outer(paste(names(dimnames(object$data))[1], rownames(ass$row), sep=""),
                                 colnames(ass$row), paste2)
      }

      if(dim(ass$col)[3] > 1) {
          rownames(col) <- outer(paste(names(dimnames(object$data))[2], rownames(ass$col), sep=""),
                                 outer(colnames(ass$col),
                                       paste(names(dimnames(object$data))[3], dimnames(ass$col)[[3]], sep=""),
                                       paste2),
                                 paste2)
      }
      else {
          rownames(col) <- outer(paste(names(dimnames(object$data))[2], rownames(ass$col), sep=""),
                                 colnames(ass$col), paste2)
      }

      rbind(phi, row, col)
  }
}

summary.assocmod <- function(object, weighting, ...) {
  if(missing(weighting))
      weighting <- if(length(object[["assoc"]]) > 0) object$assoc$weighting
                   else if(length(object$assoc.hmskew) > 0) object$assoc.hmskew$weighting
                   else if(length(object$assoc.yrcskew) > 0) object$assoc.yrcskew$weighting
                   else "marginal"
  else
      weighting <- match.arg(weighting, c("marginal", "uniform", "none"))

  coefficients <- matrix(NA, 0, 4)

  if(length(object[["assoc"]]) > 0) {
      if(weighting == object$assoc$weighting)
          coefficients <- assocTable(object, object$assoc)
      else
          coefficients <- assocTable(object, assoc(object, weighting=weighting))

      diagonal <- object$assoc$diagonal
  }

  if(length(object$assoc.hmskew) > 0) {
      if(weighting == object$assoc.hmskew$weighting)
          coefficients <- rbind(coefficients, assocTable(object, object$assoc.hmskew))
      else
          coefficients <- rbind(coefficients, assocTable(object, assoc(object, weighting=weighting)))

      diagonal <- object$assoc.hmskew$diagonal
  }

  if(length(object$assoc.yrcskew) > 0) {
      if(weighting == object$assoc.yrcskew$weighting)
          coefficients <- rbind(coefficients, assocTable(object, object$assoc.yrcskew))
      else
          coefficients <- rbind(coefficients, assocTable(object, assoc(object, weighting=weighting)))

      diagonal <- object$assoc.yrcskew$diagonal
  }

  # Not filled when several layers are present
  if(all(is.na(coefficients[,"Adjusted"])))
      coefficients <- coefficients[,!colnames(coefficients) == "Adjusted"]

  summ.gnm <- summary(asGnm(object))

  res <- list(call=object$call,
              deviance.resid=residuals(object, type="deviance"),
              chisq=sum(na.omit(c(residuals(object, "pearson")^2))),
              dissim=sum(na.omit(abs(c(residuals(object, "response")))))/sum(na.omit(abs(c(fitted(object)))))/2,
              coefficients=coefficients,
              diagonal=diagonal,
              weighting=weighting,
              deviance=object$deviance, df.residual=object$df.residual,
              bic=object$deviance - log(sum(na.omit(c(object$data)))) * object$df.residual,
              aic=object$deviance - 2 * object$df.residual,
              # Needed by R functions like anova()
              family=object$family, dispersion=summ.gnm$dispersion, df=summ.gnm$df)

  class(res) <- "summary.assocmod"

  res
}

    

print.summary.assocmod <- function(x, digits=max(3, getOption("digits") - 4), ...) {
  printModelHeading(x, digits)

  cat("\nAssociation coefficients:\n")

  if(all(is.na(x$coefficients[,"Std. error"]))) {
      printCoefmat(x$coefficients[, !colnames(x$coefficients) %in% c("Std. error", "Pr(>|z|)"), drop=FALSE],
                  print.gap=2, ...)
      cat('\nNo standard errors available\n(jackknife and bootstrap disabled, or weighting changed since fitting).\n')
  }
  else {
      printCoefmat(x$coefficients, has.Pvalue=TRUE, print.gap=2, ...)
  }

  if(length(x$diagonal) > 0) {
      cat("\nDiagonal coefficients:\n")
      print(format(x$diagonal[1:nrow(x$diagonal),], digits=digits, ...), quote=FALSE)
  }

  cat("\nNormalization weights: ", x$weighting,
      "\nDeviance:              ", format(x$deviance, digits),
      "\nPearson chi-squared:   ", format(x$chisq, digits),
      "\nDissimilarity index:   ", format(x$dissim * 100, digits), "%",
      "\nResidual df:           ", x$df.residual,
      "\nBIC:                   ", x$bic,
      "\nAIC:                   ", x$aic, "\n", sep="")
}

