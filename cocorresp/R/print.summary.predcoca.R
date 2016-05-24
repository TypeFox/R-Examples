"print.summary.predcoca" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nPredictive Co-Correspondence Analysis\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    if(!is.null(x$lambda))
      {
        cat("\nEigenvalues:\n")
        print(round(x$lambda, digits), ..., print.gap = 3)
      } else {
        pcentX <- (x$varianceExp$Xblock / x$totalVar$Xblock) * 100
        Xvar.mat <- rbind(pcentX, cumsum(pcentX))
        pcentY <- (x$varianceExp$Yblock / x$totalVar$Yblock) * 100
        Yvar.mat <- rbind(pcentY, cumsum(pcentY))
        rownames(Xvar.mat) <- rownames(Yvar.mat) <- c("Individual:", "Cumulative:")
        cat("\nPercentage Variance Explained:\n")
        cat("\nY-block: variance explained in", x$namY,
            "(response) \n", sep = " ")
        print(round(Yvar.mat, digits), ..., print.gap = 2)
        cat("\nX-block: variance explained in", x$namX,
            "(predictor) \n", sep = " ")
        print(round(Xvar.mat, digits), ..., print.gap = 2)
      }
    if (!is.null(x$cocaScores$species)) {
      cat("\n")
      writeLines(strwrap("Species scores are the regression weights for Y1 and Y2"))
      cat("\nSpecies scores:", x$namY, "\n")
      print(x$cocaScores$species$Y, digits = digits, ..., print.gap = 2)
      cat("\nSpecies scores:", x$namX, "\n")
      print(x$cocaScores$species$X, digits = digits, ..., print.gap = 2)
    }
    if (!is.null(x$cocaScores$site)) {
      cat("\n")
      writeLines(strwrap("Site scores are weighted averages of the species scores"))
      cat("\nSite scores:", x$namY, "\n")
      print(x$cocaScores$site$Y, digits = digits, ..., print.gap = 2)
      cat("\nSite scores:", x$namX, "\n")
      print(x$cocaScores$site$X, digits = digits, ..., print.gap = 2)
    }
    if(!is.null(x$loadings)) {
      cat("\n")
      writeLines(strwrap("Loadings are weighted regression coefficients"))
      cat("\n")
      writeLines(strwrap("\nLoadings Y: coefficients of weighted regression of
Y on site scores for X"))
      print(x$loadings$Y, digits = digits, ..., print.gap = 2)
      cat("\n")
      writeLines(strwrap("\nLoadings X: coefficients of weighted regression of
X on site scores for X"))
      print(x$loadings$X, digits = digits, ..., print.gap = 2)
    }
    invisible(x)
  }

