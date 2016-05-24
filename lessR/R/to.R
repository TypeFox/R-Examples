to <-
function(prefix, until, from=1, same.size=TRUE) {


  if (missing(prefix)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Specify the characters to precede each name with: prefix\n\n")
  }

  if (missing(until)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Specify the last name generated in the sequence with:  until\n\n")
  }

  if (from > until) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
       "The value of  until  must be greater than the value of  from .\n\n")
  }

  cstr <- character(length=0)
  for (ichar in (from:until)) {

    if (same.size) inum <- until else inum <- ichar
    nc <- nchar(as.character(inum))

    cc <- as.character(.fmtc(ichar, w=nc))
    for (i in 1:nchar(cc)) if (substr(cc,i,i) == " ") substr(cc,i,i) <- "0"
    cc <- paste(prefix, cc, sep="")
    cstr <- paste(cstr, cc)
  }

  nc <- nchar(cstr)
  cstr2 <- character(length=0)
  for (i in 2:nc) cstr2 <- paste(cstr2, substr(cstr,i,i), sep="")
  cstr2 <- strsplit(cstr2, " ")[[1]]

  return(cstr2)

}

