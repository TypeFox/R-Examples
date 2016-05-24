print.GPCE.sparse <- function(x, ...) {
  #cat("\nModel runs:", deparse(length(x$Output)), "\n", sep = "")
  if (!is.null(x$Sensitivity)) {
    Names = character(0)
    jmax=x$Args$jmax
    ValueS=character(0)
    ValueST=character(0)
    for (nn in 1:jmax){
      tmp = combn(jmax,nn)
      for (ii in 1:length(tmp[1,])){
        Names = append(Names,paste(as.character(tmp[,ii]),collapse = ""))
      }
    }
    for (i in 1:jmax){
      ValueS=cbind(ValueS,t(x$Sensitivity$S[[i]][1,]))
      ValueST=cbind(ValueST,t(x$Sensitivity$ST[[i]][1,]))  
    }
    cat("Sensitiviy indices:","\n\n")
    print(rbind(ValueS,Names))
    cat("\n\n")
    cat("Total Sensitivity indices:","\n\n")
    print(rbind(ValueST,Names))
  }
}
