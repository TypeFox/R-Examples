print.cpg.perm <-
function(x,...) {
  cat("\nThe permutation P-values, number of permutations and seed:\n")
  print(x$perm.p.values)
  cat("\nOther information:\n")
  temp<-data.frame(x$info[1,1:4],num.real.Holm=nrow(x$Holm.sig))
  if(!is.null(x$FDR.sig)) {
    temp<-data.frame(temp,num.real.fdr=nrow(x$FDR.sig))
    }
  print(temp)
  if(nrow(x$results)<10) {
    cat("The original analysis:\n")
    print(x$results) }
  else {
     cat("\nThe top ten CpG sites were:\n") 
     if(sum(order(x$results$P.value)==1:nrow(x$results))==nrow(x$results)) {
         print(x$results[1:10,])
          }
     else {
         print(x$results[order(x$results$P.value),][1:10,])        
              }
     cat("\nTo access results for all ", nrow(x$results)," CpG sites use object$results\n",
        "or sort(object)$results to obtain results sorted by p-value.\n")             
        }
  cat("\n")
  cat(temp$num.real.Holm, "sites were found significant by the Holm method\n")
  cat(temp$num.real.fdr, "sites were found significant by", x$info$FDR.method,"method\n")
  cat("\nThe beta values were taken from:",as.character(x$info$betainfo),"\n")
  stopholder<-which(names(x)=="coefficients")
  cat("Other attributes are:",paste(names(x)[1:floor(stopholder/2)],collapse=", "), ",\n",
         paste(paste(names(x)[(floor(stopholder/2)+1):stopholder],collapse=", "),".",sep=""),"\nThey can be accessed using the $\n")
     invisible(x)
             }
