summary.gcrq <-
function(object, digits = max(3, getOption("digits") - 3), ...){
      cat("\n***Noncrossing regression quantiles via P-splines***\n")
      cat("\nCall:\n")
      print(object$call)

      n<-nrow(as.matrix(object$fitted.values))
      p<-nrow(as.matrix(object$coefficients ))
      n.tau<-ncol(as.matrix(object$coefficients))
      sic<- sum(log(object$rho/n)) +log(n)*sum(object$df)/(2*n)
      

      if(!is.null(object$boot.coef)) list.vcov<-vcov.gcrq(object)
      for(j in 1:n.tau){
          est<-as.matrix(object$coefficients)[,j]
          if(!is.null(object$boot.coef)) {
              se<-sqrt(diag(list.vcov[[j]])) 
              ris<-cbind(Est=est, StErr=se, ratio=round(est/se,2))
              } else {
              ris<-cbind(Est=est)
              }
          rownames(ris)<-rownames(as.matrix(object$coefficients))
          cat("\nPercentile:", object$taus[j], "  Check function: ", object$rho[j] , "\n")
#          cat("Coefficients:\n")
          print(ris)
          }          
          cat("\nNo. of obs:", n, "  Check function =", round(sum(object$rho),digits), "  SIC =", round(sic,digits),"\n")
          cat("No. of est. params:", p,"(for each curve);", p*n.tau,"(total)\n")       
          }
       