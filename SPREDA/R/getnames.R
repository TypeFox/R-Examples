getnames <-
function(obj){
  tmp=paste(as.character(obj), collapse=" + ")
  
  tt=unique(unlist(strsplit(tmp," + ", fixed=TRUE)))
  if(any(tt=="+") | any(tt=="|")){
    id=as.numeric(tt=="+"|tt=="|")
    return(tt[-which(id!=0)])
  } else{
    return(tt)
  }
}
