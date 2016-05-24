cat_identifier <-
function(data,cat_index){
  test=data
  for (i in 1:length(cat_index)){
    test[,cat_index[i]]=as.character(test[,cat_index[i]])
  }
  return(test)
}
