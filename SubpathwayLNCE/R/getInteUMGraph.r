getInteUMGraph<-function(LncGenePairs){
if(!exists("envData")) envData<-initialize()
g2<-get("g2",envir=envData)
InterUMGraph<-getInteGraphList(g2,LncGenePairs)
return(InterUMGraph)
}
