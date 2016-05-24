

GetExampleData<-function(exampleData){

if(!exists("envData")) envData<-initializeESEA()




if (exampleData=="dataset")
{
dataset<-get("dataset",envir=envData)

return(dataset)
}



if(exampleData=="class.labels")
{

class.labels<-get("class.labels",envir=envData)

return(class.labels)
}


if (exampleData=="controlcharactor")
{
controlcharactor<-get("controlcharactor",envir=envData)
return(controlcharactor)
}

if (exampleData=="edgesbackgrand")
{
edgesbackgrand<-get("edgesbackgrand",envir=envData)
return(edgesbackgrand)
}

if (exampleData=="pathwayEdge.db")
{
pathwayEdge.db<-get("pathwayEdge.db",envir=envData)
return(pathwayEdge.db)
}



if (exampleData=="PathwayResult")
{
PathwayResult<-get("PathwayResult",envir=envData)
return(PathwayResult)
}

}