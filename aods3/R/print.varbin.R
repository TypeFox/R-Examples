print.varbin <- function(x, ...) {                   
	tab <- x$tab
  nam <- rownames(tab)
  tab$varmu <- sqrt(tab$varmu)	
  names(tab)[2] <- "se"
	feat <- x$features
	cat("N = ", feat["N"], " clusters, n = ", feat["n"], " subjects, m = ", feat["m"],
		" cases.\n", sep = "")
	cat("\nProportion:\n")
  List <- lapply(tab,
                 function(x) ifelse(is.na(x),
                                    "",
                                    format(round(x, digits = 4), nsmall = 3)))
	summ <- as.data.frame(t(do.call("rbind", List)))
  rownames(summ) <- nam
  print(summ)

	cat("\nalpha = ", x$alpha, " for the CIs; R = ", length(x$muboot),
      " samples for the bootstrap estimates.\n", sep = "")
  invisible(tab)
  }
