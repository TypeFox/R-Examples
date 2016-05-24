GetMiRTargetData<-function(){

if(!exists("envData")) envData<-initializeMiRSEA()

MiRTarget<-get("expMir2Tar",envir=envData)
return(MiRTarget)

}