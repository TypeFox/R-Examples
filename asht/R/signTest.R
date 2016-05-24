signTest <-
function(x,stat=c("pos-neg","prop pos"),...){
     stat<-match.arg(stat)
     n<-length(x)
     yg<-length(x[x>0])
     yless<-length(x[x<0])
     m<-yg+yless
     b<-binom.exact(yg,m,...)
     b$statistic<- yg/m
     names(b$statistic)<-"proportion positive"
     b$parameter<- yless/m
     names(b$parameter)<-"proportion negative"   
     b$method<-"Exact Sign Test"
     if (stat=="pos-neg"){
         tr<-function(p){ 2*p-1 }
         statName<-"prop pos minus prop neg"
         b$conf.int<-tr(b$conf.int)
         b$null.value<-tr(b$null.value)
         names(b$null.value)<-statName
         b$estimate<-tr(b$estimate)
         names(b$estimate)<-statName
     } else if (stat=="prop pos"){
         statName<-"proportion positive"
         names(b$null.value)<-statName
         names(b$estimate)<-statName
     }
     b
}
