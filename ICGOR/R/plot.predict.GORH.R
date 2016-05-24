plot.predict.GORH <-
function(x,...){
     arg<-list(...)
     plot(x$SurvTime,x$SurvProb,type="l",...)
     class(x)<-"predict.GORH"
}
