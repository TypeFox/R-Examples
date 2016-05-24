ROC_plot <-
function(matrixZ1,matrixZ2,matrixZ_chain,plot_name,result_file_name,burn){

#  library(ROCR)


   G1<-dim(matrixZ1)[1]
   G2<-dim(matrixZ2)[1]
   K<-dim(matrixZ1)[2]

   total_iter<-dim(matrixZ_chain[[1]])[1]
   est_Z1<-apply(matrixZ_chain[[1]][burn:total_iter,],2,mean)
   est_Z2<-apply(matrixZ_chain[[2]][burn:total_iter,],2,mean)


   pdf(plot_name,20,10)
   par(mfrow=c(1,2))
   pred1 <- prediction(est_Z1, as.vector(matrixZ1))
   perf1 <- performance(pred1, "tpr", "fpr")
   auc1 <- performance(pred1, "auc")
   plot(perf1, avg='threshold', spread.estimate='stddev', colorize=TRUE, main="matrixZ1")
   pred2 <- prediction(est_Z2, as.vector(matrixZ2))
   perf2 <- performance(pred2, "tpr", "fpr")
   auc2 <- performance(pred2, "auc")
   plot(perf2, avg='threshold', spread.estimate='stddev', colorize=TRUE, main="matrixZ2")
   dev.off()


   save(est_Z1,est_Z2,pred1,perf1,pred2,perf2,auc1,auc2,file=result_file_name)

}

