## S3 'aggregate' method for object of class 'DALY'
## aggregate by outcomes, by age/sex, or by both

aggregate.DALY <-
function(x, by = c("total", "class", "outcome"), ...){
  ## Evaluate 'by'
  by <- match.arg(by)

  ## Sum of lists
  lsum <-
  function(l, sub){
    s <- l[[1]][sub][[1]]
    if (nOutcomes > 1)
      for (i in seq(2, nOutcomes))
        s <- s + l[[i]][sub][[1]]
    return(s)
  }
  
  ## Helper objects
  nOutcomes <- length(x) - 2
  nameList <- c("DALY", "YLD", "YLL", "cases", "deaths")
  
  ## Aggregate by outcome
  if (by == "outcome"){
    out <- vector("list", nOutcomes)
    out["pop"] <- x["pop"]
    out["name"] <- x["name"]
    for (i in seq(nOutcomes)){
      out[[i]] <- vector("list", 5)
      names(out[[i]]) <- nameList
      out[[i]]["name"] <- x[[i]]["name"]
      for (j in seq(5)){
        out[[i]][[j]] <- rowSums(x[[i]][[j]])
      }
    }
  }

  ## Aggregate by class
  if (by == "class"){
    out <- vector("list", 5)
    names(out) <- nameList
    out["pop"] <- x["pop"]
    out["name"] <- x["name"]
    for (i in seq(5)){
      out[[i]] <- lsum(x, nameList[i])
    }
  }

  ## Aggregate by outcome & class
  if (by == "total"){
    y <- aggregate(x, by = "class")
    out <- vector("list", 5)
    names(out) <- nameList
    out["pop"] <- x["pop"]
    out["name"] <- x["name"]
    for (i in seq(5)){
      out[[i]] <- rowSums(y[[i]])
    }
  }

  ## Return output
  return(invisible(out))
}
