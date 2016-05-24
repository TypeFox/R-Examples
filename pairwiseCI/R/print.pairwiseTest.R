"print.pairwiseTest" <-
function( x , digits=4, ...)
{
 byout <- x$byout
 bynames <- x$bynames
 method <- x$method
 pargs<-list(...)
 pargs$digits<-digits
 METHOD <- attr(byout[[1]], "method")


if(length(byout) != length(bynames)) {stop("INTERNAL: bynames and byout of different length!! ")}

cat("P-values calculated using","\n", METHOD, "\n")
cat("\n")

 for (i in 1:length(byout))
  {
   pdat <- data.frame(
    p.value=round( byout[[i]][,"p.value", drop=FALSE], digits=digits)#,
#    groupx=byout[[i]]$groupx,
#    groupy=byout[[i]]$groupy,
#    by=rep(bynames[i],length(byout[[i]]$p.value))
   )

     cat(" ","\n")
     if(length(byout)>1)
      {cat(" BY GROUP ", bynames[i], "\n")}
     pargs$x<-pdat
     do.call("print", pargs)
     cat(" ","\n")
  }

invisible(x)
}

