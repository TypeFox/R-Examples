
# generic function to print an S3 object of class "lm.br"


print.lm.br  <-  
  function ( x, digits = max(3L, getOption("digits") - 3L), ... )  {

# add 'type' after 'formula'
  chr <- paste(deparse(x$call), sep = "\n", collapse ="\n")
  chrs <- strsplit(chr,NULL)[[1]]
  nc <- length(chrs)
  commas <- which(chrs == ",")
  xcall <- character(0)
  if(length(commas)==0) { 
    xcall <- paste( c(chrs[1:(nc-1)], ", type = \"LL\")"), collapse="")
  } else {
    c1 <- commas[1]
    if(chrs[c1+2] == "t") {
      chrs[c1+10] <- toupper(chrs[c1+10])
      chrs[c1+11] <- toupper(chrs[c1+11])
      xcall <- paste(chrs[1:nc],collapse="")
    } else {
      xcall <- 
        paste( c(chrs[1:(c1-1)], ", type = \"LL\"", chrs[c1:nc]), collapse="")
    }
  }
  cat("\nCall:\n", xcall, "\n\n", sep = "")

  if (length(coef(x))>1  && !is.na(x$coef[2]))  {

# print coefficients unless 'sety' has been called
    par <- x$CppObj$param()
    if( !par[6] )  {
      cat( "Changepoint and coefficients:\n" )
      print.default( round(x$coef, 5) )
    }
    else  cat( "Use 'mle()' for parameter estimates after a call to 'sety'\n" )

    cat( "\nSignificance Level of H0:\"no changepoint\" vs",
      "H1:\"one changepoint\"\n" )
    cat("  ")
    mx1 <- max( x$x1 )
    mn1 <- min( x$x1 )
    ai <- (mx1-mn1)/( length( x$x1 ) - 1 )
    type <- x$type
    if( type=='LT' && x$xint )
      x$sl( round( min( mx1 + ai*1.5 ), 2) )
    else
      if( type=='TL' && !x$xint )
        x$sl( -Inf )
      else
        if( type=='LT' && !x$xint )
          x$sl( Inf )
        else
          x$sl( round( max( mn1 - ai*1.5 ), 2) )
    cat("\n")
    x$ci()
  }
  else  
    cat( "No coefficients\n" )

  cat("\n")
  invisible(x)
}

