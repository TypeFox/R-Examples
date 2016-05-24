addCommas<-function(str){
  #if str includes spaces but isn't surrounded by commas then the commas are added
  
  if( length(grep(" ",str))>0  & length(grep("\"",str))==0 ){
    str <- paste("\"",str,"\"",sep="")
  }
    
  return(str)
  
}