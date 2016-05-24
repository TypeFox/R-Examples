# $Id: print.LD.R 395 2005-10-04 23:43:31Z warnes $

print.LD <- function(x, digits=getOption("digits"), ...)
  {
    saveopt <- options("digits")
    options(digits=digits)
    cat("\n")
    cat("Pairwise LD\n")
    cat("-----------\n")

    est <- t(as.matrix( c(D=x$"D","D'"=x$"D'","Corr"=x$"r")))
    rownames(est) <- "Estimates:"
    print(est)
    cat("\n")

    test <- t(as.matrix( c("X^2"=x$"X^2", "P-value"=x$"P-value",
                           "N"=x$"n") ) )
    rownames(test) <- "LD Test:"
    print(test)
    cat("\n")

    options(saveopt)
    invisible(x)
  }


summary.LD.data.frame <- function(object, digits=getOption("digits"),
                                which=c("D", "D'", "r", "X^2",
                                        "P-value", "n", " "),
                                rowsep, show.all=FALSE,
                                ...)
  {

    if(missing(rowsep))
      if(length(which)==1)
        rowsep <- NULL
      else
        rowsep <- " "

    if(is.null(rowsep))
      blank <- NULL
    else
      blank <- matrix(rowsep, ncol=ncol(object$"D"), nrow=nrow(object$"D"))
    


    saveopt <- options("digits")
    options(digits=digits)

    
    pdat <- list()
    for(name in which)
         pdat[[name]] <- object[[name]]
    
    tab <- interleave(
                      "D" = if('D' %in% names(pdat)) pdat$D else NULL,
                      "D'" = pdat$"D'",
                      "Corr." = pdat$"r",
                      "X^2"= pdat$"X^2",
                      "P-value" = pdat$"P-value",
                      "n" = pdat$"n",
                      " "=blank,
                      sep=" "
                      )

    statlist <- which[ ! (which %in% c("P-value", "n", " ") ) ]
    statlist[statlist=="X^2"] <- "X\\^2"

    formatlist <- sapply( statlist, function(object) grep(object, rownames(tab) ) )
    formatlist <- unique(sort(unlist(formatlist)))
    
    pvallist   <- grep( "P-value", rownames(tab) )
    
    tab[formatlist,] <- formatC(as.numeric(tab[formatlist,]), digits=digits,
                                format="f")
    tab[pvallist,] <- apply(object$"P-value", c(1,2),
                            function(object)trim(format.pval(object, digits=digits)))
    
    tab[trim(tab)=="NA"] <- NA

    if(!show.all)
      {
         # drop blank row/column
        entrylen <- nrow(tab)/nrow(object$n)
        tab <- tab[1:(nrow(tab) - entrylen),-1]
      }

    
    options(saveopt)
    class(tab) <- "summary.LD.data.frame"
    tab
  }

print.summary.LD.data.frame <- function(x, digits=getOption("digits"), ...)
{
  cat("\n")
  cat("Pairwise LD\n")
  cat("-----------\n")

  print(as.matrix(unclass(x)), digits=digits, quote=FALSE, 
        na.print="    ", right=TRUE) 
        
  cat("\n")

  invisible(x)

  
}


print.LD.data.frame <- function(x, ...)
  print(summary(x))
