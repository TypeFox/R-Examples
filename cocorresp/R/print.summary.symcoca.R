"print.summary.symcoca" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nSymmetric Co-Correspondence Analysis\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat("\nEigenvalues:\n")
    print(round(x$lambda, digits), ..., print.gap = 3)
    cat("\nInertia:\n")
    inert <- rbind(unlist(x$inertia$total),
                   unlist(x$inertia$residual))
    explained <- inert[1,] - inert[2,]
    inert <- rbind(inert, explained)
    rownames(inert) <- c("Total:", "Residual:", "Explained:")
    colnames(inert) <- c(x$namY, x$namX)
    printCoefmat(inert, digits = digits, na.print = "")
    if (!is.null(x$cocaScores$species)) {
      if (x$scaling == 2) {
        cat("\n")
        writeLines(strwrap("Species scores for axis i are rescaled by the
quarter root of the eigenvalue for axis i."))
      }
      cat("\nSpecies scores for", x$namY, "are weighted averages of site scores for", x$namX, "\n")
      cat("\nSpecies Scores: ", x$namY, "\n")
      print(x$cocaScores$species$Y, digits = digits, ..., print.gap = 2)
      cat("\nSpecies scores for", x$namX, "are weighted averages of site scores for", x$namY, "\n")
      cat("\nSpecies Scores: ", x$namX, "\n")
      print(x$cocaScores$species$X, digits = digits, ..., print.gap = 2)
    }
    if (!is.null(x$cocaScores$site)) {
      cat("\nSite scores for", x$namY, "are weighted averages of species scores for", x$namX, "\n")
      cat("\nSite Scores: ", x$namY, "\n")
      print(x$cocaScores$site$Y, digits = digits, ..., print.gap = 2)
      cat("\nSite scores for", x$namX, "are weighted averages of species scores for", x$namY, "\n")
      cat("\nSite Scores: ", x$namX, "\n")
      print(x$cocaScores$site$X, digits = digits, ..., print.gap = 2)
    }
    if (!is.null(x$cocaScores$loadings)) {
      cat("\nLoadings for", x$namY, "are coefficients of weighted regression of", x$namY, "on the X-matrix\n")
      cat("\nLoadings: ", x$namY, "\n")
      print(x$cocaScores$loadings$Y, digits = digits, ..., print.gap = 2)
      cat("\nLoadingsfor", x$namX, "are coefficients of weighted regression of", x$namX, "on the X-matrix\n")
      cat("\nLoadings: ", x$namX, "\n")
      print(x$cocaScores$loadings$X, digits = digits, ..., print.gap = 2)
    }
    if (!is.null(x$cocaScores$xmatrix)) {
      cat("\n")
      writeLines(strwrap(paste("X-Matrix = (X1 + X2) / 2, where X1 and X2 are the site scores for", x$namY, "and", x$namX, "respectively.")))
      cat("\n")
      print(x$cocaScores$xmatrix, digits = digits, ..., print.gap = 2)
    }
    cat("\n")
    invisible(x)
  }

