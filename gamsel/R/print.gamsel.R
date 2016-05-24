print.gamsel=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")
      out=cbind(summarynz(x)[,c(3,1,2)],signif(x$dev.ratio,digits),signif(x$lambda,digits))
      colnames(out)=c("NonZero","Lin","NonLin","%Dev","Lambda")
      print(out,...)
      invisible(out)
    }
