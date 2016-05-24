`mu.AND` <- 
function(GE, frml=NULL)
{
  if (!is.matrix(GE)) GE <- as.matrix(GE)
  if (ncol(GE)<2) return(GE)
  if (is.null(frml))
  {
    if (ncol(GE)>100)
    stop("number of variables in one comparison should be ",
         "no more than 100. now it is ", ncol(GE))
    GE  <- sq.array(GE)
    AND <- apply(GE, 1:2, any)
    nNA <- AND[,1]*0
    for (i in 1:dim(GE)[3])
    {
      nNA <- nNA + diag(GEi <- GE[,,i])
      AND <- AND * (GEi + (1-GEi)*(1-t(GEi)))
    }
    return(as.numeric(AND))
  }

  if(substring(frml, 1, 1)!="(")
    frml <- paste("(", frml, ")", sep="")
  tkn <- unlist(strsplit(frml,""))
  nok <- attr(regexpr("[0-9,()]+",frml),"match.length")
  if ( nok < nchar(frml)) print(paste(collapse="",
    "illegal character \"", tkn[nok+1], "\" in 'frml'"))

  tmp <- matrix(0, dim(GE)[1]+1, sum(tkn=="(")+1 )
  FstFree <- function(tmp) match(TRUE, tmp[1,]==0, nomatch=0)

  level <- i <- 0
  while ((i <- i+1) <= nok) {
    switch( tkn[i],
      "(" = level <- level + 1,
      "," = next,
      ")" = { tmp[1, use <- (tmp[1,]==level)] <- 0       # flag for reuse
              tmp[, FstFree(tmp)] <- c(level <- level-1, mu.AND(tmp[-1, use])) },
      { num <- as.numeric(substring(
               frml,i,i<-i-1+regexpr("[,)]",substring(frml,i+1))))
        if ((FstTmp <- FstFree(tmp)) == 0)
          FstTmp <- ncol(tmp <- cbind(tmp, 0))                        # add col
          tmp[,FstTmp] <- c(level, GE[, num]) }

    )
  }
  return(tmp[-1,1])
}
