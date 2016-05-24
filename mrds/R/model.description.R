model.description <- function(model){
# builds the description of the model to be printed
# only used internally.

  key <- switch(model$key,
                hn="Half-normal",
                hr="Hazard-rate",
                th1="Threshold 1",
                th2="Threshold 2")

  if(is.null(key)){
    key <- "Uniform"
  }

  mod.str <- paste(key,"key function")
  if(!is.null(model$adjustment)){
    adj.series <- switch(model$adjustment$series,
                         cos="cosine",
                         herm="Hermite polynomial",
                         poly="simple polynomial")
    mod.str <- paste(mod.str,"with",adj.series,"adjustment term")

    adj.order <- model$adjustment$order
    if(length(adj.order)>1){
      mod.str <- paste(mod.str,"s",sep="")
    }
    mod.str <- paste(mod.str,"of order",paste(adj.order,collapse=","))
  }

  return(mod.str)
}
