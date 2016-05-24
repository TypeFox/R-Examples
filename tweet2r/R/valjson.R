#' @export

valjson<-function(fileprefix){

  #get number of json files form current directory
  numfiles=length(list.files(pattern=".json",all.files= ))
  #initilitze loop counter
  counter=0
  #initialize number of deleted files
  del_files=0
  for (i in 1:numfiles){
  
    #get actual file name
    file_name=paste(fileprefix,counter,".json",sep="")
    
    #get the file size
    file.info(file_name)$size
    
    #delete the file if it has no tweets (size<1024 kb)
    if(file.info(file_name)$size<1024){
      #update counter
      counter=counter+1
      del_files=del_files+1
      #remove file with no tweets
      file.remove(file_name)
    }
    else{
      file_num=counter-del_files
      #rename_file
      file.rename(from=file_name,to=paste(fileprefix,file_num,".json",sep=""))
      #update counter
      counter=counter+1
    }
  } 
  
  return (paste("Number of files deleted:",del_files, sep=" "))

}

  



