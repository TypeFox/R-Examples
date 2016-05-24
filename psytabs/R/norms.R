norms <- 
  function (sumscores, from, to, statistics = "PR") {
    sumscores.range <- from:to
    xecdf <- stats::ecdf(sumscores)
    sumscores.z <- (sumscores.range - mean(sumscores, na.rm=TRUE))/stats::sd(sumscores, na.rm=TRUE)
    norm.table <- data.frame(Score = sumscores.range) 
    
    if (!is.na(statistics[1])) {
      if ("PR" %in% statistics) {
        sumscores.percentranks <- round(xecdf(sumscores.range)*100, 1)
        norm.table <- cbind(norm.table, PR = sumscores.percentranks)
      }
      if ("T" %in% statistics) {
        sumscores.t <- round(50 + 10*sumscores.z, 1) 
        norm.table <- cbind(norm.table, T = sumscores.t)
      }
      if ("Stanine" %in% statistics) {
        sumscores.stanine <- trunc(5 + sumscores.z*2)
        sumscores.stanine[sumscores.stanine<1] <- 1
        sumscores.stanine[sumscores.stanine>9] <- 9
        norm.table <- cbind(norm.table, STANINE = sumscores.stanine)
      }
      if ("IQ" %in% statistics) {
        sumscores.iq <- round(100 + 15*sumscores.z, 1)
        norm.table <- cbind(norm.table, IQ = sumscores.iq)
      }
      if ("Z" %in% statistics) {
        sumscores.Z <- round(100 + 10*sumscores.z, 1)
        norm.table <- cbind(norm.table, Z = sumscores.Z)
      }
      if ("z" %in% statistics) {
        norm.table <- cbind(norm.table, z = sumscores.z) 
      }
    } else {
      print("Warning")
    }
    return(norm.table)
  }