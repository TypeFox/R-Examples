print.cpg <-
function(x,...) {
  if(nrow(x$results) <10) {
     cat("\nThe results are:\n")
      print(x$results)
      }
     else{
     cat("\nThe top ten CpG sites were:\n") 
     if(sum(order(x$results$P.value)==1:nrow(x$results))==nrow(x$results)) {
         print(x$results[1:10,])
          }
     else {
         print(x$results[order(x$results$P.value),][1:10,])        
              }
     
    cat("\nTo access results for all ", nrow(x$results)," CpG sites use object$results",
        "\nor sort(object)$results to obtain results sorted by p-value.\n")
      
         }
      cat("\nGeneral info:\n")

      temp<-data.frame(x$info[c(1:5,7)],num.Holm=nrow(x$Holm.sig),num.fdr=nrow(x$FDR.sig))
      print(temp)
      cat("\n")
      cat(temp$num.Holm, "sites were found significant by the Holm method\n")
      cat(temp$num.fdr, "sites were found significant by", x$info$FDR.method,"method\n")
      cat("\nThe beta values were taken from:",as.character(x$info$betainfo),"\n")
      if(!is.factor(x$indep)) {
        cat("Effect sizes and standard error can be accessed using $coefficients\n")
      }
      cat("Other attributes are:",paste(names(x)[1:(length(x)-1)],collapse=", "), "\n",
            "They can be accessed using the $\n")
     invisible(x)
      }
