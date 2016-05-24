print.GPCE.lar <- function(x, ...) {
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
    print("Sensitiviy indices:",...)
    print(rbind(ValueS,Names),...)
    print("Total Sensitivity indices:",...)
    print(rbind(ValueST,Names),...)
  }
}
