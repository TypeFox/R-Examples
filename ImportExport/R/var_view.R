var_view <-
function(x){
  varname<-names(x)
  descr<-c(1:length(x));value_labels<-c();formats<-c()
  for(i in 1:length(x)){value_labels[i]<-paste(names(attr(x[,i],"value.labels")),attr(x[,i],"value.labels"),sep="=",collapse="/")
  formats[i]<-paste(class(x[1,i]),collapse=" ")
  if(is.null(attr(x[,i],"label"))){descr[i]<-"No descr"}
  else{descr[i]<-attr(x[,i],"label")}
  }
  return(data.frame(varname,descr,value_labels,formats))
}
