partition.ch=function(quantify.ch.fun,t,breaks,include.lowest=T,type=c("list","index"), ... ){
  
  ch.list=unlist(list.historylabels(t))
  names(ch.list)=ch.list
  names(ch.list)[1]="empty"
  ch.quantified=sapply(ch.list,FUN=quantify.ch.fun, ... )
partition.index=cut(ch.quantified,breaks=breaks,
                      include.lowest=include.lowest)
  
  
  if(type[1]=="list"){
    out=split(unlist(ch.list),partition.index)
  }else{
    out=as.numeric(partition.index)
    names(out)=names(ch.list)
  }
  
  return(out)
  
}

