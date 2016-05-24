#' @encoding UTF-8
#' @title Frequency Table
#'
#' @description Simulating the FREQ procedure of SPSS.
#'
#' @param .data The data.frame.
#' @param x A column for which a frequency of values is desired.
#' @param verbose A logical value, if \code{TRUE}, extra statistics are also provided.
#' @param \dots Additional arguements (currently ignored)
#'
#' @seealso \code{\link{freq}}, \code{\link{crosstable}}.
#'
#' @examples
#' data(cathedrals)
#'
#' Frequency(cathedrals, Type)
#'
#' cathedrals %>% Frequency(Height)
#'
#' @importFrom stats sd
#' @rdname Frequency
#' @aliases Freq
#' @export
`Frequency` <- function(.data, x, verbose=TRUE, ...) UseMethod("Frequency")


#' @rdname Frequency
#' @export
`Frequency.default` <- function(.data, x, verbose=TRUE, ...) {
  vec <-eval(substitute(x), .data, parent.frame())
  nmiss=sum(is.na(vec))
  fsum=summary(factor(vec))
  ftab=cbind(fsum,100*fsum/sum(fsum))
  if (nmiss==0) {
    ftab=cbind(ftab,100*cumsum(fsum)/sum(fsum))
    colnames(ftab)=c("Frequency"," Valid Percent"," Cum Percent")
    ftab[,2] <- round(ftab[,2],2)
    ftab[,3] <- round(ftab[,3],2)
    if(verbose==FALSE){
      return(ftab)
    }
    print(ftab)
  }
  else
  {
    ftab=cbind(ftab,100*fsum/sum(fsum[1:(length(fsum)-1)]),100*cumsum(fsum)/sum(fsum[1:(length(fsum)-1)]))
    ftab[length(fsum),3:4]=NA
    ftab[,2] <- round(ftab[,2],2)
    ftab[,3] <- round(ftab[,3],2)
    ftab[,4] <- round(ftab[,4],2)
    cat("\n")
    cat("--------------------------------------------------------\n")
    colnames(ftab)=c("Frequency","   Percent"," Valid Percent"," Cum Percent")
    if (dim(ftab)[1]==length(levels(vec)))
    {
      rownames(ftab)[1:length(levels(factor(vec)))]=levels(factor(vec))
    }
    if(verbose==FALSE){
      return(ftab)
    }
    print(ftab)
  }
  cat("--------------------------------------------------------\n")
  cat("Total",rep(" ",8-trunc(log10(sum(fsum)))),sum(fsum),"\n",sep="")
  cat("\n")

  if (length(attributes(vec)$class) != 0)
  {
    if ("factor" %in% attributes(vec)$class)
    {
      cat("Warning: Statistics may not be meaningful for factors!\n\n")
    }
  }
  s1=cbind(mean(as.numeric(vec),na.rm=TRUE),stats::sd(as.numeric(vec),na.rm=TRUE))
  rownames(s1)=" "
  colnames(s1)=c("       Mean","        Std dev")
  print(s1)
  s2=cbind(min(as.numeric(vec),na.rm=TRUE),max(as.numeric(vec),na.rm=TRUE))
  rownames(s2)=" "
  colnames(s2)=c("    Minimum","        Maximum")
  print(s2)
  s3=cbind(sum(!is.na(vec)),nmiss)
  rownames(s3)=" "
  colnames(s3)=c("Valid cases",  "Missing cases")
  print(s3)
  cat("\n")
}#--end of Freq
NULL

#' @export
#' @rdname Frequency
Freq <- Frequency
