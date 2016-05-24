#' Extracts the model description
#' 
#' Returns a description of the model fitted in the ddf object.
#' 
#' @param model a ddf object
#' @return mod.str a string descibing the fitted model
#' @author Jeff Laake & Laura Marshall
model.description<-function(model){
# builds the description of the model to be printed
# only used internally.
# taken from mrds
  if(!is.null(model$ds$aux$ddfobj$type)){
    key<-switch(model$ds$aux$ddfobj$type,
                hn="Half-normal",
                hr="Hazard-rate")
    mod.str<-paste(key,"key function")
    if(!is.null(model$ds$aux$ddfobj$adjustment)){
      adj.series<-switch(model$ds$aux$ddfobj$adjustment$series,
                         cos="cosine",
                         herm="Hermite polynomial",
                         poly="simple polynomial")
      mod.str<-paste(mod.str,"with",adj.series,"adjustment term")
  
      adj.order<-model$ds$aux$ddfobj$adjustment$order
      if(length(adj.order)>1){
        mod.str<-paste(mod.str,"s",sep="")
      }
      mod.str<-paste(mod.str,"of order",paste(adj.order,collapse=","))
    }
  }else{
    mod.str <- model$mr$call
    mod.str$formula <- model$mr$formula
  }
  return(mod.str)
}
