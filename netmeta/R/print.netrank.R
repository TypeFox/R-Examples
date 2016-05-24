print.netrank <- function(x,
                          sort=TRUE,
                          digits=max(4, .Options$digits - 3),
                          ...){
  
  meta:::chklogical(sort)
  meta:::chknumeric(digits, single=TRUE)
  
  if (sort)
    res <- as.data.frame(round(x$Pscore[order(-x$Pscore)], digits))
  else
    res <- as.data.frame(round(x$Pscore, digits))
  ##
  colnames(res) <- "P-score"
  ##
  matitle(x)
  ##
  print(res, digits=digits, ...)
  
  invisible(NULL)
}
