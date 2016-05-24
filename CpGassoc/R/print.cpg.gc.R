print.cpg.gc <-
function(x,...) {

  cat("Using genomic control adjustment the top sites are:\n")
      print(x$gc.results[order(x$gc.results$Adjust.P.value),][1:10,])
      
  cat("\nGeneral info:\n")

  print(x$gc.info)
  cat("\n")
  cat(x$gc.info$num.holm, "sites were found significant by the Holm method\n")
  cat(x$gc.info$num.fdr, "sites were found significant by", x$gc.info$FDR.method,"method\n")
  }
