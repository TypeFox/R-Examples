print.mitml.anova <- function(x,...){
# print method for anova method

  cl <- x$call
  test <- x$test
  fml <- x$formula
  mth <- x$method
  use <- x$use
  reml <- x$reml
  m <- x$test[[1]]$m

  # header
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")

  cat("\nModel comparison calculated from",m,"imputed data sets.")
  cat("\nCombination method:",mth,
    if(mth=="D2"){paste("(",use,")",sep="")},"\n")

  # model formulas
  cat("\n")
  for(mm in 1:length(fml)) cat("Model ",mm,": ",fml[mm],"\n", sep="")

  # check for very large values
  out <- sapply(test, function(z) z$test)
  fmt <- c("%.3f","%.0f","%.3f","%.3f","%.3f")
  ln <- apply(out, 1, function(z) any(z>=10^5))
  fmt[ln] <- "%.3e"
  out <- apply(out, 2, function(z) sprintf(fmt,z))

  # model comparisons
  cat("\n")
  nt <- length(test)
  out <- matrix(out,ncol=nt)
  comp <- paste0(1:nt, " vs ", 2:(nt+1),":")
  out <- rbind(comp, out)
  w <- max(sapply(c(out,colnames(test[[1]]$test)),nchar))
  cat("",format(c("",colnames(test[[1]]$test)),justify="right",width=w),"\n")
  for(mm in 1:ncol(out)) cat("",format(out[,mm],justify="right",width=w),"\n")

  if(reml){
    cat("\nModels originally fit with REML were automatically refit using ML.\n")
  }

  cat("\n")

  invisible()
}

