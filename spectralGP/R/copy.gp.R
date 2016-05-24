"copy.gp" <-
function(object,object2=NULL,...){
  nullFlag=FALSE
  if(is.null(object2)){
    object2=new.env(parent=globalenv())
    class(object2)=class(object)
    nullFlag=TRUE
  }
  elements=names(object)
  for(index in 1:length(elements)){
    assign(elements[index],get(elements[index],envir=object,inherits=FALSE),envir=object2)
  }
  if(nullFlag){
    return(object2)
  } else{
    return(NULL)
  }
}
