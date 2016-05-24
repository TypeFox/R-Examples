initializePAGIExample<-function(){
   utils::data("ExampleData",package="PAGI")
   
}

getdataset<-function(){

if(!exists("ExampleData")) initializePAGIExample()

newdataset<-get("dataset",envir=ExampleData)

return(newdataset)
}

getclass.labels<-function(){

if(!exists("ExampleData")) initializePAGIExample()

newclass.labels<-get("class.labels",envir=ExampleData)

return(newclass.labels)
}
