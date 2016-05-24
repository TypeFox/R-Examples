## UNIDIFF model (Erikson & Goldthorpe, 1992), or uniform layer effect model (Xie, 1992)

unidiff <- function(tab, diagonal=c("included", "excluded", "only"),
                    constrain="auto",
                    weighting=c("marginal", "uniform", "none"), norm=2,
                    family=poisson,
                    tolerance=1e-8, iterMax=5000,
                    trace=FALSE, verbose=TRUE,
                    checkEstimability=TRUE, ...) {
  diagonal <- match.arg(diagonal)
  weighting <- match.arg(weighting)

  tab <- as.table(tab)

  if(length(dim(tab)) < 3)
      stop("tab must have (at least) three dimensions")

  if(diagonal != "included" &&
     (nrow(tab) != ncol(tab) || !all(rownames(tab) == colnames(tab))))
     stop(sprintf("diagonal = %s is only supported for square tables with identical row and column categories",
                  diagonal))

  if(length(dim(tab)) > 3)
      tab <- margin.table(tab, 1:3)

  tab <- prepareTable(tab, TRUE)
  vars <- names(dimnames(tab))


  f <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s",
               vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3])

  if(diagonal == "included") {
      f <- sprintf("%s + Mult(Exp(%s), %s:%s)", f, vars[3], vars[1], vars[2])

      if(identical(constrain, "auto"))
          constrain <- sprintf("^Mult\\(Exp\\(.\\), \\Q%s:%s\\E\\)\\Q.%s%s\\E$",
                               vars[1], vars[2], vars[3], dimnames(tab)[[3]][1])
  }
  else if(diagonal == "excluded") {
      f <- sprintf("%s + %s:Diag(%s, %s) + Mult(Exp(%s), %s:%s)",
                   f, vars[3], vars[1], vars[2], vars[3], vars[1], vars[2])

      if(identical(constrain, "auto"))
          constrain <- sprintf("^Mult\\(Exp\\(.\\), \\Q%s:%s\\E\\)\\Q.%s%s\\E$",
                               vars[1], vars[2], vars[3], dimnames(tab)[[3]][1])
  }
  else if(diagonal == "only") {
      f <- sprintf("%s + Mult(Exp(%s), Diag(%s, %s))", f, vars[3], vars[1], vars[2])

      if(identical(constrain, "auto"))
          constrain <- sprintf("^Mult\\(Exp\\(.\\), Diag\\(\\Q%s, %s\\E\\)\\)\\Q.%s%s\\E$",
                               vars[1], vars[2], vars[3], dimnames(tab)[[3]][1])
  }

  eliminate <- eval(parse(text=sprintf("quote(%s:%s)", vars[1], vars[3])))

  # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
  args <- list(formula=as.formula(f), data=tab, constrain=constrain,
               family=family, eliminate=eliminate, tolerance=tolerance, iterMax=iterMax,
               trace=trace, verbose=verbose)

  model <- do.call("gnm", c(args, list(...)))

  if(is.null(model))
      return(NULL)

  model$unidiff <- list()
  
  model$unidiff$layer <- getContrasts(model, pickCoef(model, sprintf("Mult\\(Exp\\(\\.\\)", vars[3])),
                                      "first", check=checkEstimability)


  nr <- nrow(tab)
  nc <- ncol(tab)

  if(diagonal == "included") {
      if(weighting == "marginal") {
          p <- prop.table(tab)
          rp <- margin.table(p, 1)
          cp <- margin.table(p, 2)
      }
      else if(weighting == "none") {
          rp <- rep(1, nr)
          cp <- rep(1, nc)
      }
      else {
          rp <- rep(1/nr, nr)
          cp <- rep(1/nc, nc)
      }

      rp1 <- prop.table(rp)
      cp1 <- prop.table(cp)

      mat <- matrix(0, nr, nc)
      con <- matrix(0, length(coef(model)), length(mat))
      ind <- pickCoef(model, sprintf("Mult\\(Exp\\(\\Q%s\\E\\)", vars[3]))
      for(i in 1:nr) {
          for(j in 1:nc) {
              mat[] <- 0
              mat[i,] <- -cp1
              mat[,j] <- -rp1
              mat[i,j] <- 1 - rp1[i] - cp1[j]
              mat <- mat + rp1 %o% cp1
              con[ind, (j - 1) * nr + i] <- mat
          }
      }

      colnames(con) <- names(ind)
      model$unidiff$interaction <- gnm::se(model, con, checkEstimability=checkEstimability)

      if(weighting == "marginal") {
          model$unidiff$phi <- maor(fitted(model)[,,1], TRUE, norm=norm,
                                    row.weights=rp, col.weights=cp)
          model$unidiff$maor <- maor(fitted(model)[,,1], norm=norm,
                                    row.weights=rp, col.weights=cp)
      }
      else {
          model$unidiff$phi <- maor(fitted(model)[,,1], TRUE, norm=norm, weighting=weighting)
          model$unidiff$maor <- maor(fitted(model)[,,1], norm=norm, weighting=weighting)
      }
  }
  else if(diagonal == "only"){
      if(weighting != "uniform")
          warning("uniform weigthing is always used when diagonal=\"only\"")

      # Quasi-variances cannot be computed for these coefficients, so hide the warning
      # Also skip the reference level
      suppressMessages(model$unidiff$interaction <- getContrasts(model, pickCoef(model, sprintf("Mult\\(Exp\\(\\Q%s\\E\\)", vars[3])),
                                                                 ref="first", check=checkEstimability)$qvframe[-1,])

      # Diag() sorts levels alphabetically, which is not practical
      model$unidiff$interaction[c(1 + order(rownames(tab))),] <- model$unidiff$interaction
  }
  else {
     # Interaction coefficients cannot be identified for now
     model$unidiff$interaction <- NULL
  }

  rownames(model$unidiff$layer$qvframe) <- dimnames(tab)[[3]]

  if(diagonal != "excluded")
      rownames(model$unidiff$interaction) <- gsub(sprintf("Mult\\(Exp\\(\\Q%s\\E\\), \\.\\)\\.(Diag\\(%s, %s\\))?",
                                                          vars[3], vars[1], vars[2]), "",
                                                  rownames(model$unidiff$interaction))

  model$unidiff$diagonal <- diagonal
  model$unidiff$weighting <- weighting

  class(model) <- c("unidiff", class(model))

  model$call.gnm <- model$call
  model$call <- match.call()

  model
}

print.unidiff <- function(x, digits=max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  cat("\nLayer coefficients:\n")
  print(setNames(exp(x$unidiff$layer$qvframe$estimate), row.names(x$unidiff$layer$qvframe)),
        digits=digits, print.gap=2, ...)

  if(length(x$unidiff$phi) > 0) {
      cat("\nLayer intrinsic association coefficients:\n")
      print(setNames(exp(x$unidiff$layer$qvframe$estimate) * x$unidiff$phi, row.names(x$unidiff$layer$qvframe)),
            digits=digits, print.gap=2, ...)
  }

  if(x$unidiff$diagonal == "included") {
      cat("\nFull two-way interaction coefficients:\n")
      interaction <- x$data[,,1]
      interaction[] <- x$unidiff$interaction[,1]
  }
  else if(x$unidiff$diagonal == "only") {
      cat("\nDiagonal interaction coefficients:\n")
      interaction <- x$unidiff$interaction[,1]
      names(interaction) <- rownames(x$unidiff$interaction)
  }
  else {
      cat("\nDiagonal interaction coefficients:\n")
      interaction <- "Not supported."
  }

  print.default(format(interaction, digits=digits, ...), quote=FALSE, print.gap=2)

  cat("\nNormalization weights:", x$unidiff$weighting)

  printModelStats(x)

  invisible(x)
}

summary.unidiff <- function(object, ...) {
  layer <- object$unidiff$layer$qvframe[,-4]
  phi <- object$unidiff$phi
  interaction <- object$unidiff$interaction

  layer <- cbind(exp(layer[,1]), layer, 2 * pnorm(-abs(layer[,1]/layer[,2])))
  colnames(layer) <- c("Exp(Estimate)", "Estimate", "Std. Error", "Quasi SE", "Pr(>|z|)")
  rownames(layer) <- paste(names(dimnames(object$data))[3], rownames(layer), sep="")

  if(length(phi) > 0) {
      phi.layer <- layer
      phi.layer[,1] <- layer[,1] * phi
      phi.layer[,2] <- layer[,2] + log(phi)
      phi.layer[,3:4] <- layer[,3:4] * phi
  }

  if(object$unidiff$diagonal != "excluded") {
      interaction <- cbind(interaction, 2 * pnorm(-abs(interaction[,1]/interaction[,2])))
      colnames(interaction) <- c("Estimate", "Std. Error", "Pr(>|z|)")
  }

  summ.gnm <- summary(asGnm(object))

  res <- list(call=object$call,
              deviance.resid=residuals(object, type="deviance"),
              chisq=sum(na.omit(c(residuals(object, "pearson")^2))),
              dissim=sum(na.omit(c(abs(residuals(object, "response")))))/sum(na.omit(c(abs(fitted(object)))))/2,
              layer=layer, phi.layer=phi.layer, interaction=interaction,
              diagonal=object$unidiff$diagonal,
              weighting=object$unidiff$weighting,
              deviance=object$deviance, df.residual=object$df.residual,
              bic=object$deviance - log(sum(na.omit(c(object$data)))) * object$df.residual,
              aic=object$deviance - 2 * object$df.residual,
              # Needed by R functions like anova()
              family=object$family, dispersion=summ.gnm$dispersion, df=summ.gnm$df)

  class(res) <- "summary.unidiff"

  res
}

print.summary.unidiff <- function(x, digits=max(3, getOption("digits") - 4), ...) {
  printModelHeading(x, digits)

  cat("\nLayer coefficients:\n")
  printCoefmat(x$layer, digits, signif.legend=FALSE, print.gap=2, ...)

  if(length(x$phi) > 0) {
      cat("\nLayer phi association coefficients:\n")
      printCoefmat(x$phi.layer, digits, signif.legend=FALSE, print.gap=2, ...)
  }

  if(x$diagonal != "only")
      cat("\nFull two-way interaction coefficients:\n")
  else
      cat("\nDiagonal interaction coefficients:\n")

  if(x$diagonal != "excluded")
      printCoefmat(x$interaction, digits, has.Pvalue=TRUE, print.gap=2, ...)
  else
      cat("Not supported.\n")


  cat("\nNormalization weights: ", x$weighting,
      "\nDeviance:              ", format(x$deviance, digits),
      "\nPearson chi-squared:   ", format(x$chisq, digits),
      "\nDissimilarity index:   ", format(x$dissim * 100, digits), "%",
      "\nResidual df:           ", x$df.residual,
      "\nBIC:                   ", x$aic,
      "\nAIC:                   ", x$bic, "\n", sep="")
}

plot.unidiff <- function(x, what=c("layer.coef", "phi", "maor"),
                         se.type=c("quasi.se", "se"), conf.int=.95,
                         numeric.auto=TRUE, type="p",
                         xlab=names(dimnames(x$data))[3], ylab=NULL,
                         add=FALSE, ylim, ...) {
  if(!inherits(x, "unidiff"))
      stop("x must be a unidiff object")

  what <- match.arg(what)
  se.type <- match.arg(se.type)

  if(is.null(ylab)) {
      if(what == "layer.coef")
          ylab <- "Log odds ratio coefficient"
      else if(what == "phi")
          ylab <- "Intrinsic association coefficient"
      else
          ylab <- "Mean absolute odds ratio"
  }

  qv <- x$unidiff$layer$qvframe

  w <- qnorm((1 - conf.int)/2, lower.tail=FALSE)

  coefs <- qv$estimate

  if(se.type == "quasi.se") {
      tops <- qv$estimate + w * qv$quasiSE
      tails <- qv$estimate - w * qv$quasiSE
   }
   else {
      tops <- qv$estimate + w * qv$SE
      tails <- qv$estimate - w * qv$SE
   }

  coefs <- exp(coefs)
  tops <- exp(tops)
  tails <- exp(tails)

  if(what == "phi") {
      coefs <- coefs * x$unidiff$phi
      tops <- tops * x$unidiff$phi
      tails <- tails * x$unidiff$phi
  }
  else if(what == "maor") {
      coefs <- x$unidiff$maor^coefs
      tops <- x$unidiff$maor^tops
      tails <- x$unidiff$maor^tails
  }

  range <- max(tops) - min(tails)

  if(missing(ylim))
      ylim <- c(min(tails) - range/20, max(tops) + range/20)

  # plot() converts x coordinates to numeric if possible, but segments
  # needs a real coordinate, so convert directly
  if(numeric.auto && !any(is.na(as.numeric(rownames(qv)))))
      at <- as.numeric(rownames(qv))
  else
      at <- factor(rownames(qv), levels=rownames(qv))

  if(!isTRUE(add)) {
      plot.default(at, coefs, type=type, xaxt="n",
                   ylim=ylim, xlab=xlab, ylab=ylab, ...)
      axis(1, at, labels=at)
  }
  else {
      points.default(at, coefs, type=type, ...)
  }

  segments(as.numeric(at), tails, as.numeric(at), tops)

  invisible(x)
}
