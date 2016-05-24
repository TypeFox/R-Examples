"print.pairwiseCI" <-
function( x , digits=4, ...)

{
byout <- x$byout
bynames <- x$bynames
method <- x$method
pargs<-list(...)
pargs$digits<-digits

METHOD<-x$byout[[1]]$method

  cat(" ","\n")
  cat(format(x$conf.level*100, digits=digits), "%-confidence intervals", "\n")
  cat(" Method: ", METHOD, "\n")
  cat(" ","\n")

if(length(byout) != length(bynames)) {stop("INTERNAL: bynames and byout of different length!! ")}

 for (i in 1:length(byout))
  {
   CItable <- round( cbind(byout[[i]]$estimate, byout[[i]]$lower, byout[[i]]$upper), digits=digits)
   colnames(CItable) <- c("estimate", "lower", "upper")
   rownames(CItable ) <- byout[[i]]$compnames

     cat(" ","\n")
     if(length(byout)>1)
      {cat(" BY GROUP ", bynames[i], "\n")}
     pargs$x<-CItable
     do.call("print", pargs)
     cat(" ","\n")
     cat(" ","\n")
  }

invisible(x)
}

