normTable <- 
  function (sumscores, from, to, statistics = "PR", group = NA, as.list = FALSE) {
    if(!is.na(group[1])) {
      sumscores.list <- by(sumscores, group, function(x) norms(x, from = from, to = to, statistics = statistics))
      if(as.list == TRUE) {
        return(sumscores.list)
      }
      sumscores.data <- do.call("cbind", sumscores.list)
      if (length(dim(sumscores.list)) > 1) {
        names(sumscores.data) <- paste(rep(levels(interaction(dimnames(sumscores.list)[[1]], dimnames(sumscores.list)[[2]])), each=length(statistics)+1), names(sumscores.data))
      }
      return(sumscores.data)
    } else {
      norms(sumscores = sumscores, from = from, to = to, statistics = statistics)
    }
  }