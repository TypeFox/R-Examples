SO2Simple<-function (sd) 
{
indices<-as.numeric(as.vector(sd$variables[sd$variables[,4]=="IC",1]))
sd$indivIC[,indices,]
}