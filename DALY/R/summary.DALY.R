## S3 'summary' method for 'DALY' class

summary.DALY <-
function(object, relative = FALSE, outcomes = FALSE, digits = 0, ...){
  ## Summarize total DALYs
  y <- aggregate(object, by = "class")

  total <- vector("list", 5)
  nameList <- c("DALY", "YLD", "YLL", "cases", "deaths")
  
  total[[1]] <- round(colMeans(y$DALY), digits)
  total[[2]] <- round(colMeans(y$YLD), digits)
  total[[3]] <- round(colMeans(y$YLL), digits)
  total[[4]] <- round(colMeans(y$cases), digits)
  total[[5]] <- round(colMeans(y$deaths), digits)
  names(total) <- nameList

  ## Summarize DALYs per outcome
  nOutcomes <- length(object) - 2
  out <- vector("list", nOutcomes)
  
  for (i in seq(nOutcomes)){
    out[[i]] <- vector("list", 5)
	names(out)[i] <- object[[i]]$name
	names(out[[i]]) <- nameList
    for (j in seq(5)){
      out[[i]][[j]] <-
        round(colMeans(object[[i]][nameList[j]][[1]]), digits)
    }
  }
  
  ## Denominator
  if (relative){
    denom <- object$pop / 1000
  } else {
    denom <- 1
  }
  
  ## Print summaries
  cat("\nDALY Calculator: ", object$name, "\n\n")
  
  if (relative)
    cat("Total population: ", sum(object$pop), "\n\n")
  
  if (outcomes && nOutcomes > 1){
    for (i in seq(nOutcomes)){
      cat(object[[i]]$name, "\n\n")
	  for (j in seq(5)){
	    cat(nameList[j], "\n")
        print(round(out[[i]][[j]] / denom, digits))
	    cat("\n")
	  }
      cat("\n")
    }
  } else {
    for (i in seq(5)){
	  cat(nameList[i], "\n")
	  print(round(total[[i]] / denom, digits))
	  cat("\n")
	}
  }

  cat("\n")
  
  return(invisible(list(total = total, outcomes = out)))
}
