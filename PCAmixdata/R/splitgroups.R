splitgroups<-function(base,groups,name.groups){
  nbr.groups<-length(unique(groups))
  list.groups<-list()
  for (i in 1:nbr.groups){
    list.groups[[i]]<-data.frame(base[,which(groups==i)],check.names=T)
    colnames(list.groups[[i]])<-colnames(base)[(groups==i)]
  }
  names(list.groups)<-name.groups
  return(list.groups)
}
