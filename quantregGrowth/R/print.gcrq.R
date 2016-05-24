print.gcrq <-
function(x, digits = max(3, getOption("digits") - 3), ...){
      n<-nrow(as.matrix(x$fitted.values))
#      cat("\n***Noncrossing regression quantiles via P-splines***\n")
#      cat("\nCall:\n")
      cat("Call:\n")
      print(x$call)
      n.tau<-length(x$taus)
      if(n.tau>1) {
          p<-nrow(as.matrix(x$coefficients))
          cat("\nNo. of obs:", n, "  No. of estimated parameters:", p,"(for each curve);", p*n.tau,"(total)\n")
          cat("Quantile curves (",length(x$taus),") at percentiles: ", x$taus, "\n")
          cat("Check functions:", round(x$rho,2), "\n")
          sic<- sum(log(x$rho/n)) +log(n)*sum(x$df)/(2*n)
          cat("Overall Check function =", round(sum(x$rho),digits), "  SIC =", round(sic,digits),"\n")
          } else {
          p<-length(x$coefficients)
          cat("\nNo. of obs:", n, "  No. of estimated parameters:", p,"\n")
          cat("Quantile curves at percentile: ", x$taus, "\n")
          sic<- sum(log(x$rho/n)) +log(n)*sum(x$df)/(2*n)
          cat("Check function =", round(sum(x$rho),digits), "  SIC =", round(sic,digits),"\n")
          }
      }
