initialize<-function(){
   utils::data("envData",package="SubpathwayLNCE")
}

Getenvir<-function(envData){

if(!exists("envData")) initialize()
return(get(envData,envir=envData))

}