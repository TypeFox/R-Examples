

GetExampleData<-function(dataset="prostate"){

if(!exists("MPINetData")) initializeMPINet()




if (dataset=="prostate")
{
prostateTMrisk<-get("prostateTMrisk",envir=MPINetData)

return(prostateTMrisk)
}



if(dataset=="diabetes1")
{

diabetes1risk<-get("diabetes1",envir=MPINetData)

return(diabetes1risk)
}


if (dataset=="diabetes2")
{
diabetes2risk<-get("diabetes2",envir=MPINetData)
return(diabetes2risk)
}

}