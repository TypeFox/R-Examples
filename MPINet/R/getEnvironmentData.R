initializeMPINet<-function(){
   utils::data("MPINetData",package="MPINet")
}

Getenvir<-function(envirdata){

if(!exists("MPINetData")) initializeMPINet()
return(get(envirdata,envir=MPINetData))

}