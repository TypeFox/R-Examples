print.cooccur <-
function(x, ...){
  ptab <- x$results
  cat("Call:\n")
  print(x$call)
  if("omitted" %in% names(x)){
    cat(paste("\nOf ",x$pot_pairs," species pair combinations, ",x$omitted," pairs (",round(x$omitted / x$pot_pairs * 100,2)," %) were removed from the analysis because expected co-occurrence was < 1 and ",sep=""))
    cat(paste(x$pairs," pairs were analyzed","\n",sep=""))
  }else{
  cat(paste("\n",x$pairs," pairs were analyzed","\n",sep=""))
  }
  cat("\nCooccurrence Table:\n")
  print(ptab[ptab$p_gt <= 0.05 | ptab$p_lt <= 0.05,])
  #NextMethod("print")
}
