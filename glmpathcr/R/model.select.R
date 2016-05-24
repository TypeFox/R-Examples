model.select <-
function(object,which="BIC") {
	summary.glmpath<-summary(object)
	if (c("BIC","AIC")[charmatch(which,c("BIC","AIC"))]=="BIC") {
		min<-which.min(summary.glmpath[,which])
		}
	if (c("BIC","AIC")[charmatch(which,c("BIC","AIC"))]=="AIC") {
		min<-which.min(summary.glmpath[,which])
		}
	best<-as.numeric(gsub("Step ","",rownames(summary.glmpath)[min]))
   best
}

