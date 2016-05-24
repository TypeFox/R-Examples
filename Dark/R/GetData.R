GetData<- function(path,.....){
  # where 'path' is file/path/to/data
  if(missing (path)){
    tmp<- TestData(0:20)
    path='~'
  }

  # it might be easier to change the working directory and then perform searches. 
  Home<- getwd() 
  setwd(path)
  # choose data, e.g. subject 4
  # read in the data using read.table() or read.csv()
    # then create list
    obj<- NULL
    # these lines will need editing in your version
  
    obj$time <- tmp$time # the time data   
    obj$thrs <- tmp$thrs # the time data 
    obj$data <- tmp$data # the name of the data, eg 'study_xp_sx4'
  # add other items if known, see the help ?GetData
  
  # then return to the original directory
  setwd(Home)
  
  # set the class of the object this is used for some ancillary functions
  class(obj)='dark' 
  return(obj)
}