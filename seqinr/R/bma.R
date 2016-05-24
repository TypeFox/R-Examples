bma <- function(nucl, warn.non.IUPAC = TRUE, type = c("DNA", "RNA")){
  if(nchar(nucl[1]) != 1) stop("vector of single chars expected")
  type <- match.arg(type)
  nucl <- tolower(nucl)
  nucl <- unlist(sapply(nucl, amb, checkBase = FALSE))
  iupac <- sapply(amb(), amb)
  if(type == "DNA"){
    iupac$u <- NULL
  } else {
    iupac$t <- NULL
  }
  idx <- unlist(lapply(iupac, setequal, nucl))
  if(all(idx == FALSE)){
    if(warn.non.IUPAC){
      warning(paste("Undefined IUPAC symbol with:", paste(nucl, collapse = " ")))
    }
    return(NA)
  }
  return(names(iupac)[idx])
}
