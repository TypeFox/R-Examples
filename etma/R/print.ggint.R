print.ggint <-
function (x, ...) {
  cat("Epistasis Test in Meta-Analysis (ETMA)\n")
  cat("A MCMC algorithm for detecting gene-gene interaction in meta-analysis.\n\n")
  cat(paste0("This analysis include ",x$df+3," studies. (df = ",x$df,") \n\n"))

  out.matrix=data.frame(matrix(NA,3,7))
  colnames(out.matrix)=c("b","se","OR",paste0(round((1-x$sig.level)*100,2),"%ci.l"),paste0(round((1-x$sig.level)*100,2),"%ci.u"),"t value","p value")
  rownames(out.matrix)=c("SNP1(mutation)","SNP2(mutation)","Interaction")
  out.matrix[,1]=formatC(x$b,5,format="f")
  out.matrix[,2]=formatC(x$se,5,format="f")
  out.matrix[,3]=formatC(x$OR,3,format="f")
  out.matrix[,4]=formatC(x$ci.l,3,format="f")
  out.matrix[,5]=formatC(x$ci.u,3,format="f")
  out.matrix[,6]=formatC(x$t,4,format="f")
  out.matrix[,7]=.ROUND.p(x$pval,4)
  print(out.matrix)
}
