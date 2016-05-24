print.logforest <-
function(x, ...)
{
 if (class(x)!= "logforest")
    stop("n not of class logforest")
 num<-x$numout
 pscore<-sort(x$Predictor.importance, decreasing=TRUE)
 pred.freq<-x$Predictor.frequency[names(pscore)[1:num]]
 piscore<-sort(x$PI.importance, decreasing=TRUE)
 vimp.freq<-x$PI.frequency[names(piscore)[1:num]]
 p.1st<-paste("Top ", num, " Predictors", sep="")
 pi.1st<-paste("Top ", num, " Interactions", sep="")
 if (x$norm==TRUE)
   {
   pred.score<-round(pscore[1:num]/pscore[1], digits=4)
   vimp.score<-round(piscore[1:num]/piscore[1], digits=4)
   p.cnames<-c(p.1st, "Normalized Predictor Importance","Frequency")
   pi.cnames<-c(pi.1st, "Normalized Interaction Importance","Frequency")
   }
 else {
   pred.score<-round(pscore[1:num], digits=2)
   vimp.score<-round(piscore[1:num], digits=2)
   p.cnames<-c(p.1st, "Predictor Importance","Frequency")
   pi.cnames<-c(pi.1st, "Interaction Importance","Frequency")
   }
 prds<-cbind(names(pred.score), pred.score, pred.freq)  
 pis<-cbind(names(vimp.score), vimp.score, vimp.freq)
 colnames(prds)<-p.cnames
 colnames(pis)<-pi.cnames
 rownames(pis)<-rownames(prds)<-c(1:num)
 cat("Number of logic regression trees =", length(x$AllFits), sep=" ")
 cat("\n")
 cat("\n")
 cat(num, " most important predictors \n", sep="")
 cat("\n")
 print.default(prds, quote=FALSE, print.gap=3)
 cat("\n")
 cat(num, " most important interactions \n", sep="")
 cat("\n")
 print.default(pis, quote=FALSE, print.gap=3)
}
