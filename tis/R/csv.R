csv <- function(z, file = "", noDates = FALSE, row.names = !is.tis(z), ...){
  if(!inherits(file, "connection")) {
    if(file == "") 
      file <- paste(deparse(substitute(z)), ".csv", sep = "")
    file <- gsub(".csv.csv$", ".csv", paste(file, ".csv", sep = ""))
  }
  zz <- stripTis(as.matrix(z))
  colList <- columns(zz)
  n.cl <- length(colList)
  if(!is.numeric(zz)){ ## Try turning columns into numbers
    for(i in 1:n.cl){
      zi <- colList[[i]]
      naSpots <- is.na(zi)
      zn <- as.numeric(zi)
      if(!any(is.na(zn[!naSpots])))
        colList[[i]] <- zn
    }
  }
  if(length(names(colList)) != n.cl){
    if(n.cl == 1) cnames <- substr(deparse(substitute(z)), 1, 12)
    else {
      div26 <- ceiling(n.cl/26)
      cnames <- paste(letters, rep(1:div26, each = div26), sep = "")[1:n.cl]
    }
    names(colList) <- cnames
  }
  if(is.tis(z) && !noDates)
    colList <- c(list(Date = unclass(ssDate(ti(z)))), colList)

  df <- do.call("data.frame", c(colList, list(stringsAsFactors = FALSE)))
  write.table(df, file, sep = ",", row.names = row.names,
              qmethod = "double", ...)
}
