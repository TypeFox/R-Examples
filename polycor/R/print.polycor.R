# last modified 24 June 04 by J. Fox

"print.polycor" <-
function(x, digits = max(3, getOption("digits") - 3), ...){
  if (x$type == "polychoric"){
    se <- sqrt(diag(x$var))
    se.rho <- se[1]
    est <- if (x$ML) "ML est." else "2-step est."
    cat("\nPolychoric Correlation, ", est, " = ", signif(x$rho, digits),
      " (", signif(se.rho, digits), ")", sep="")
    if (x$df > 0)
        cat("\nTest of bivariate normality: Chisquare = ", 
        signif(x$chisq, digits), ", df = ", x$df, ", p = ", 
        signif(pchisq(x$chisq, x$df, lower.tail=FALSE), digits), "\n", sep="")
    else cat("\n")
    r <- length(x$row.cuts)
    c <- length(x$col.cuts)
    if (r == 0) return(invisible(x))
    row.cuts.se <- se[2:(r+1)]
    col.cuts.se <- se[(r+2):(r+c+1)]
    rowThresh <- signif(cbind(x$row.cuts, row.cuts.se), digits)
    if (r > 1) cat("\n  Row Thresholds\n") 
    else cat("\n  Row Threshold\n") 
    colnames(rowThresh) <- c("Threshold", "Std.Err.")
    rownames(rowThresh) <- if (r > 1) 1:r else " "
    print(rowThresh)
    colThresh <- signif(cbind(x$col.cuts, col.cuts.se), digits)
    if (c > 1) cat("\n\n  Column Thresholds\n")
    else cat("\n\n  Column Threshold\n")
    colnames(colThresh) <- c("Threshold", "Std.Err.")
    rownames(colThresh) <- if (c > 1) 1:c else " "
    print(colThresh)
    }
  else if (x$type == "polyserial"){
    se <- sqrt(diag(x$var))
    se.rho <- se[1]
    est <- if (x$ML) "ML est." else "2-step est."
    cat("\nPolyserial Correlation, ", est, " = ", signif(x$rho, digits),
      " (", signif(se.rho, digits), ")", sep="")
    cat("\nTest of bivariate normality: Chisquare = ", signif(x$chisq, digits),
      ", df = ", x$df, ", p = ", signif(pchisq(x$chisq, x$df, lower.tail=FALSE), digits),
      "\n\n", sep="")
    if (length(se) == 1) return(invisible(x))
    cuts.se <- se[-1]
    thresh <- signif(rbind(x$cuts, cuts.se), digits)
    colnames(thresh) <- 1:length(x$cuts)
    rownames(thresh) <- c("Threshold", "Std.Err.")
    print(thresh)
    }
  else print(unclass(x))
  invisible(x)
  }
