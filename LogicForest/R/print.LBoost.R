print.LBoost <-
function(x, num, ...)
{
 if(class(x)!="LBoost")
   stop("not of class LBoost")
 if(missing(num)) {num<-5}
 k<-length(x$CVmod)
 nBS2<-length(x$CVmod[[1]]$AllFits)
 cat(k, "-fold training datasets, ", nBS2, " trees per training dataset")
 cat("\n")
 cat("Number of logic regression trees =", k*nBS2)
 cat("\n")
 cat("\n")
 CVerr<-x$CVmisclass$CVmiss
 cat("CV model error rate = ", CVerr)
 cat("\n")
 cat("\n")
 if(x$PredImp==T) 
   {
   toppreds<-sort(x$Pred.import, decreasing=T)[1:num]
   norm.toppred<-round(toppreds/toppreds[1], digits=2)
   topp.nms<-names(toppreds)
   pfreq<-x$Pred.freq[topp.nms]
   pc.names<-c("Predictor", "Importance", "Frequency")
   prd<-cbind(topp.nms, norm.toppred, pfreq)
   colnames(prd)<-pc.names
   rownames(prd)<-c(1:num)
   print.default(prd, quote=FALSE, print.gap=3)
   cat("\n")
   }
 if(x$PIimp=="AddRemove")
   {
   topPI<-sort(x$AddRemove.PIimport, decreasing=T)[1:num]
   norm.topPI<-round(topPI/topPI[1], digits=2)
   topPInms<-names(topPI)
   PIfreq<-x$PI.freq[topPInms]
   PIc.names<-c("Prime Implicant", "Add-in/Leave-out Importance", "Frequency")
   pis<-cbind(topPInms, norm.topPI, PIfreq)
   colnames(pis)<-PIc.names
   rownames(pis)<-c(1:num)
   print.default(pis, quote=FALSE, print.gap=3)
   }
 if(x$PIimp=="Permutation")
   {
   topPI<-sort(x$Perm.PIimport, decreasing=T)[1:num]
   topPInms<-names(topPI)
   norm.topPI<-round(topPI/topPI[1], digits=2)
   PIfreq<-x$PI.freq[topPInms]
   PIc.names<-c("Prime Implicant", "Permutation Importance", "Frequency")
   pis<-cbind(topPInms, norm.topPI, PIfreq)
   colnames(pis)<-PIc.names
   rownames(pis)<-c(1:num)
   print.default(pis, quote=FALSE, print.gap=3)
   }
 if(x$PIimp=="Both")
   {
   topPI1<-sort(x$AddRemove.PIimport, decreasing=T)[1:num]
   topPInms1<-names(topPI1)
   norm.topPI1<-round(topPI1/topPI1[1], digits=2)
   PIfreq1<-x$PI.freq[topPInms1]
   topPI2<-sort(x$Perm.PIimport, decreasing=T)[1:num]
   topPInms2<-names(topPI2)
   norm.topPI2<-round(topPI2/topPI2[1], digits=2)
   PIfreq2<-x$PI.freq[topPInms2]
   pis1<-cbind(topPInms1, norm.topPI1, PIfreq1)
   pis2<-cbind(topPInms2, norm.topPI2, PIfreq2)
   colnames(pis1)<-c("Prime Implicant", "Add-in/Leave-out Importance", "Frequency")
   colnames(pis2)<-c("Prime Implicant", "Permutation Importance", "Frequency")
   rownames(pis1)<-rownames(pis2)<-c(1:num)
print(c(1:num))
   cat("Top", num, "prime implicants using the Add-in/Leave-out measure", sep=" ")
   cat("\n")
   print.default(pis1, quote=FALSE, print.gap=3)
   cat("\n")
   cat("\n")
   cat("Top", num, "prime implicants using the permutation measure", sep=" ")
   cat("\n")
   print.default(pis2, quote=FALSE, print.gap=3)
   }
}
