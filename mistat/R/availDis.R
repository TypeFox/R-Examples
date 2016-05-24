availDis <- function(ttf, ttr, n, seed=NA, printSummary=TRUE){
  ttf <- as.numeric(ttf)
  ttr <- as.numeric(ttr)
  
  if(!is.na(seed))
    set.seed(seed)
  
  as1 <- boot::boot(data=ttf, statistic=function(x,i) mean(x[i]), R=n)$t
  as2 <- boot::boot(data=ttr, statistic=function(x,i) mean(x[i]), R=n)$t
  bootavail <- as1/(as1 + as2)
  if(printSummary){
    av1 <- mean(ttf)
    av2 <- mean(ttr)
    avail <- av1 / (av1 + av2)
    cat(paste(" The estimated MTTF from ttf is", round(av1, 2)))
    cat(paste("\n The estimated MTTR from ttr is", round(av2, 2)))
    cat(paste("\n The estimated asymptotic availability is", round(avail, 4)))
    cat("\n")
    colnames(bootavail) <- "availability EBD"
    print(summary(bootavail)) 
  } 
  invisible(as.numeric(bootavail))
}