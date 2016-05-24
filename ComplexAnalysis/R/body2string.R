body2string <-
function(func){
Body<-as.character(body(func))
if(Body[1]=="{"){Body<-Body[2:length(Body)]}else{Body<-as.character(as.expression(body(func)))}
emptyspace<-function(subfunc){return(paste(strsplit(subfunc," ")[[1]],collapse="",sep=""))};Body<-tapply(Body,1:length(Body),emptyspace)
Body<-paste(Body,collapse=";")
return(Body)
}
