print.summary.lqa <-
function (x, digits = max (3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption ("show.signif.stars"), ...)
{
    cat ("\nCall:\n")
    cat (paste (deparse (x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat ("Deviance Residuals: \n")

    if (x$df.residual > 5) 
    {
	x$deviance.resid <- quantile (x$deviance.resid, na.rm = TRUE)
	names (x$deviance.resid) <- c ("Min", "1Q", "Median", "3Q", "Max")
    }

    print.default (x$deviance.resid, digits = digits, na.print = "", print.gap = 2)

    df <- if ("df" %in% names(x)) 
            x[["df"]] 
          else 
            NULL

    if (!is.null (df) && (nsingular <- df[3L] - df[1L]))
       cat("\nCoefficients: (", nsingular," not defined because of singularities)\n", sep = "")
    else 
       cat("\nCoefficients:\n")

    coefs <- x$coefficients

    printCoefmat(coefs, P.value = FALSE, has.Pvalue = FALSE)
    
    ##
    cat("\n(Dispersion parameter for ", x$family$family," family taken to be ", format(x$dispersion), ")\n\n",
	apply(cbind(paste(format(c("Null","Residual"), justify="right"),
                          "deviance:"),
		    format(unlist(x[c("null.deviance","deviance")]),
			   digits= max(5, digits+1)), " on",
		    format(unlist(x[c("df.null","df.residual")])),
		    " degrees of freedom\n"),
	      1L, paste, collapse=" "), sep="")

    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")

    cat("AIC: ", format(x$aic, digits= max(4, digits+1)),", trace (hat matrix) = ", format(x$rank, digits= max(4, digits+1)), "\n\n")
    if (x$method == "pirls")
	cat ("Number of Fisher Scoring iterations: ", x$n.iter, "\n")
    else 
    {
       cat ("Converged: ", x$converged, " after ", x$n.iter, "boosting iterations \n")
       cat ("Minimum AIC reached after ", x$best.iter, " iterations. \n")
    }

   
    cat("\n")
    invisible(x)
}

