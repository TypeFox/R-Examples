summary.threshpt <-
function(object, ...)
{
      summary_results<-list()
      thresh_best<-object$best.fit
      dev<-thresh_best[7]
      threshold<-thresh_best[8]
      parm_coef=object$parm.coef
      myform<-object$formula

      Best_fit<-matrix(0,2,3)
      rownames(Best_fit)<-c("<Threshold", ">=Threshold")
      colnames(Best_fit)<-c("Estimate", "Std. Error","Pr(>|z|)")

      Best_fit[1,1]<-thresh_best[1]
      Best_fit[1,2]<-thresh_best[2]
      Best_fit[1,3]<-thresh_best[3]
      Best_fit[2,1]<-thresh_best[4]
      Best_fit[2,2]<-thresh_best[5]
      Best_fit[2,3]<-thresh_best[6]

	summary_results$call<-object$call
	summary_results$formula<-myform
	summary_results$deviance<-dev
	summary_results$threshold<-threshold
	summary_results$coefficients<-rbind(parm_coef[1,],rbind(Best_fit,parm_coef[-1,]))

    	class(summary_results)<-"summary.threshpt"
	summary_results
}
