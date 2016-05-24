GetEdgesBackgrandData<-function(){

if(!exists("envData")) envData<-initializeESEA()

edgesbackgrand<-get("edgesbackgrand",envir=envData)
return(edgesbackgrand)

}