plotVarUsage <-
function(mogavs,method=c("hist","plot","table")){
  tb<-NULL
  for(i in 2:nrow(mogavs$nonDominatedSet)){
  tb<-rbind(tb,cbind(which(mogavs$nonDominatedSet[i,]==1),rep(mogavs$numOfVariables[i],mogavs$numOfVariables[i])))
}
if(method=="plot"){
plot(tb,xlab="Variable number", ylab="Number of vars in model",pch="*",xlim=c(min(tb[,1]),max(tb[,1])),xaxt="n")
axis(side=1,at=c(min(tb[,1]):max(tb[,1])),labels=c(min(tb[,1]):max(tb[,1])))
}
else if(method=="hist"){
  h<-hist(tb[,1],breaks=c((min(tb[,1])-1):(max(tb[,1]))),freq=TRUE,probability=FALSE,include.lowest=T,right=T,density=NULL,labels=T,axes=F,main="Histogram of variable usage",xlab="Variable number")
  axis(side=1,at=c((min(tb[,1])):max(tb[,1])-0.5),labels=c(min(tb[,1]):max(tb[,1])))
}
else{
  table(tb[,1])
}
}
