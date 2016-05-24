setseq <- function(seq, levs){
  name <- deparse(substitute(seq))
  ##
  if (length(levs)!=length(seq))
    stop("Length of argument '", name,
         "' different from number of treatments.", call.=FALSE)
  ##
  if (length(unique(seq)) != length(seq))
    stop("Values for argument '", name,
         "' must all be disparate.", call.=FALSE)
  ##
  if (is.numeric(seq)){
    if (any(is.na(seq)))
      stop("Missing values not allowed in argument '",
           name, "'.", call.=FALSE)
    if (any(!(seq %in% (1:length(levs)))))
      stop(paste("Argument '", name,
                 "' must be a permutation of the integers from 1 to ",
                 length(levs), ".", sep=""), call.=FALSE)
    res <- levs[seq]
  }
  else if (is.character(seq)){
    if (length(unique(levs)) == length(unique(tolower(levs))))
      idx <- charmatch(tolower(seq), tolower(levs), nomatch=NA)
    else
      idx <- charmatch(seq, levs, nomatch=NA)
    if (any(is.na(idx)) || any(idx==0))
      stop(paste("Argument '", name,
                 "' must be a permutation of the following values:\n  ",
                 paste(paste("'", levs, "'", sep=""),
                       collapse=" - "), sep=""), call.=FALSE)
    res <- levs[idx]
  }
  else
    stop("Argument '", name, "' must be either a numeric or character vector.",
         call.=FALSE)
  
  res
}
