print.imperfect.trees<-function(x,...){
      if (!is.element("imperfect.trees", class(x))) 
        stop("Argument 'x' must be an object of class \"imperfect.trees\".")
cat("There are",x$ntrees,"independent trees and", dim(x$coefficients)[1],"parameters in this model. \n")

accuracy<-matrix(c(x$Seny,x$Senx,x$Espy,x$Espx),nrow=1)
colnames(accuracy)<-c("Se_R","Se_T","Sp_R","Sp_T")
prevalences<-x$Prevl
   for(i in 1:length(prevalences)){
      eval(parse(text=paste0("Prevalence_",i,"<-",prevalences[i]))) 
      }
Prevalencia<-prevalences[1]
for(i in 2:length(prevalences)){
      eval(parse(text=paste0("Prevalencia<-data.frame(Prevalencia,Prevalence_",i,")"))) 
      }
names(Prevalencia)[1]<-"Prevalence_1"
cat("\n")
cat("Estimation of the accuracy indices of both reference and the test of interest:\n")
print(accuracy,print.gap=2)
cat("\n")
cat("Estimation of prevalence in each primary study:\n")
print(Prevalencia)
cat("\n")
cat("Model fit statistics:\n")
cat("AIC=",x$aic,"\n")
cat("\n")
print(x$gof)
cat("\n")
}


