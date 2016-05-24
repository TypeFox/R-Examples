summary.ggint <-
  function(object, ...) {
    x<-object
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

    cat("\n")

    out.matrix2=out.matrix[c(2,1,3),]
    row.names(out.matrix2)=c("SNP1(wild type) & SNP2(mutation)","SNP1(mutation) & SNP2(wild type)","SNP1(mutation) & SNP2(mutation)")
    #out.matrix2[3,1]=formatC(sum(x$b),5,format="f")
    #out.matrix2[3,2]=formatC(sqrt(t(rep(1,3))%*%x$vcov%*%t(t(rep(1,3)))),5,format="f")
    out.matrix2[3,3]=formatC(exp(sum(x$b)),3,format="f")
    out.matrix2[3,4]=formatC(exp(sum(x$b)-qt(1-x$sig.level/2,x$df)*sqrt(t(rep(1,3))%*%x$vcov%*%t(t(rep(1,3))))),3,format="f")
    out.matrix2[3,5]=formatC(exp(sum(x$b)+qt(1-x$sig.level/2,x$df)*sqrt(t(rep(1,3))%*%x$vcov%*%t(t(rep(1,3))))),3,format="f")
    out.matrix2[3,6]=formatC(sum(x$b)/sqrt(t(rep(1,3))%*%x$vcov%*%t(t(rep(1,3)))),4,format="f")
    out.matrix2[3,7]=.ROUND.p(pt(abs(sum(x$b)/sqrt(t(rep(1,3))%*%x$vcov%*%t(t(rep(1,3))))),df=x$df,lower.tail=FALSE),4)
    print(out.matrix2[,3:7])
  }
