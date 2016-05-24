plot.predict.GORMC <-
function(x,...){
     arg<-list(...)
     #if(is.null(arg$add)) arg$add<-FALSE
     plot(x$Survival$Time,x$Survival$SurvProb,type="l",...)
     #if(arg$add) lines(x$Survival$Time,x$Survival$SurvProb,...)
     class(x)<-"predict.GORMC"
}
