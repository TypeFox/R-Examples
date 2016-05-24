###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
print.mvmeta <-
function(x, digits=4, ...) {
#
################################################################################
# HEADING AND SUB-HEADING
#
  # HEADING
  cat("Call:  ",paste(deparse(x$call),sep="\n",collapse="\n"),"\n\n",sep="")
#
  # SUB-HEADING
  cat("Fixed-effects coefficients:","\n",sep="")
  table <- formatC(x$coefficients,digits=digits,format="f")
  print(table,quote=FALSE,right=TRUE,print.gap=2)
  cat("\n")
#
################################################################################
# FIT STATS
#
  cat(x$dim$m," studies, ",x$df$nall," observations, ",x$df$fixed," fixed and ",
    x$df$random," random-effects parameters","\n",sep="")
  if(na <- length(x$na.action)) cat(" (",na," stud",ifelse(na>1L,"ies","y"),
    " removed due to missingness",")\n",sep="")
  if(!x$method%in%c("mm","vc")) {
    table <- c(x$logLik,AIC(x),BIC(x))
    names(table) <- c("logLik","AIC","BIC")
    table <- formatC(table,digits=digits,format="f")
    print(table,quote=FALSE,right=TRUE,print.gap=2)
  }
  cat("\n")
#
}

