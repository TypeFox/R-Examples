
GetPathwayEdgeData<-function(){

if(!exists("envData")) envData<-initializeESEA()

pathwayEdge.db<-get("pathwayEdge.db",envir=envData)
return(pathwayEdge.db)

}