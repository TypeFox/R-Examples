date_to_char<-function(x){
  
  datevar_names<-as.character(subset(var_view(x),formats=="dates times" | formats=="date")$varname)
  if(length(datevar_names) != 0){
  for(i in 1:length(datevar_names)){
    x[[datevar_names[i]]]<-as.character(x[[datevar_names[i]]])
    }
  }
  return(x)
}
