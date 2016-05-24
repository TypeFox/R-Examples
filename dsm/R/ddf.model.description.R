ddf.model.description<-function(model){
# builds the description of the model to be printed
# only used internally.

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

  return(mod.str)
}
