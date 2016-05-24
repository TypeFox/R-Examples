date_to_char <-
function(x){
   var_view<-var_view(x)
  datevar_names<-as.character(var_view[(var_view$formats=="dates times" | var_view$formats=="date"),]$varname)
  if(length(datevar_names) != 0){
  for(i in 1:length(datevar_names)){
    x[[datevar_names[i]]]<-as.character(x[[datevar_names[i]]])
    }
  }
  return(x)
}
