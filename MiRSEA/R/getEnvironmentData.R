initializeMiRSEA<-function(){
   utils::data("envData",package="MiRSEA")
}

Getenvir<-function(envData){

if(!exists("envData")) initializeMiRSEA()
return(get(envData,envir=envData))

}