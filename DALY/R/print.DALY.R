print.DALY <-
function(x, relative = FALSE, outcomes = FALSE, prob = .95, digits = 0, ...){
  ## Summarize total DALYs
  y <- aggregate(x, by = "total")

  total <- rbind(summarize(y$DALY,   .prob = prob, .digits = digits),
                 summarize(y$YLD,    .prob = prob, .digits = digits),
                 summarize(y$YLL,    .prob = prob, .digits = digits),
                 summarize(y$cases,  .prob = prob, .digits = digits),
                 summarize(y$deaths, .prob = prob, .digits = digits))
  rownames(total) <- c("DALY", "YLD", "YLL", "cases", "deaths")
  colnames(total)[1:2] <- c("Mean", "Median")

  ## Contribution YLD & YLL
  pctYLD <- round(100 * mean(y$YLD / y$DALY, na.rm = TRUE))
  pctYLL <- round(100 * mean(y$YLL / y$DALY, na.rm = TRUE))

  ## Summarize DALYs per outcome
  y <- aggregate(x, by = "outcome")
  
  nOutcomes <- length(x) - 2
  out <- vector("list", nOutcomes)
  
  for (i in seq(nOutcomes)){
    out[[i]] <- rbind(summarize(y[[i]]$DALY,   .prob = prob, .digits = digits),
                      summarize(y[[i]]$YLD,    .prob = prob, .digits = digits),
                      summarize(y[[i]]$YLL,    .prob = prob, .digits = digits),
                      summarize(y[[i]]$cases,  .prob = prob, .digits = digits),
                      summarize(y[[i]]$deaths, .prob = prob, .digits = digits))
    names(out)[i] <- x[[i]]$name
    rownames(out[[i]]) <- c("DALY", "YLD", "YLL", "cases", "deaths")
    colnames(out[[i]])[1:2] <- c("Mean", "Median")
  }

  ## Denominator
  if (relative){
    denom <- sum(x$pop) / 1000
  } else {
    denom <- 1
  }
  
  ## Print summaries
  cat("\nDALY Calculator: ", x$name, "\n\n")
  
  if (relative)
    cat("Total population: ", sum(x$pop), "\n\n")
  
  if (outcomes && nOutcomes > 1){
    for (i in seq(nOutcomes)){
      cat(x[[i]]$name, "\n")
      print(round(out[[i]] / denom, digits))
	  cat("\n")
    }
  } else {
	print(round(total / denom, digits))
	cat("\nYLD/DALY =", sprintf("%1.0f%%", pctYLD))
    cat("\nYLL/DALY =", sprintf("%1.0f%%", pctYLL), "\n")
  }

  cat("\n")

  ## Return total & outcome-wise DALYs, and contributions
  return(invisible(list(total = total, outcomes = out,
                        pct = c(pctYLD, pctYLL))))
}
