bref.variable0<-function(x,digits=max(3, getOption("digits") - 
    3)){
if(is.numeric(x)){ttt<-signif(mean(na.omit(x)),digits=digits)
names(ttt)<-"M"}
else{

if(is.factor(x)){
MAX<-rev(sort(table(na.omit(x))))[1]
ttt<-round(100*MAX/sum(table(na.omit(x))))
names(ttt)<-paste(names(MAX)," (%)",sep="")}


else{ttt<-NA}

}
ttt
}
