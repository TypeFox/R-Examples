initializeESEA<-function(){
   utils::data("envData",package="ESEA")
}

Getenvir<-function(envData){

if(!exists("envData")) initializeESEA()
return(get(envData,envir=envData))

}