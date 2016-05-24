## Save method

setMethod("savedata", "Dataclass", function(object,...){

  if(is.null(filename(object))) 
     stop("This Dataclass object has to be given a filename before it can be saved to harddisk")

  name0 <- as.character(match.call(call=sys.call(-1))$object)            
  eval.parent(parse(text=paste("save(", name0, ", file = \"", 
                                filename(object), "\")" ,sep = "")))
  
  namecomment <- paste(name0,".comment",sep="")
  filenamecomment <- paste(filename(object),".comment",sep = "")
  
  eval.parent(parse(text = paste(namecomment," <- ", name0, sep = ""))) 
  eval.parent(parse(text = paste(namecomment, "@Data <- NULL", sep = ""))) 
  eval.parent(parse(text = paste("save(", namecomment, ", file = \"", 
                                  filenamecomment,"\" )", sep = ""))) 
  eval.parent(parse(text = paste("rm(",namecomment,")", sep = "")))
})

## Load Method for comments
cload <- function(filename){
  eval.parent(parse(text=paste("load(\"",filename,".comment\")", sep = "")))
}

## Simulation

setMethod("savedata", "Simulation", function(object,...){
  if(is.null(filename(object)))
     stop("This simulation has to be given a filename before it can be saved to harddisk")

  name <- deparse(substitute(object))

  eval.parent(parse(text=paste(name,"@Data <- NULL",sep="")))
  eval.parent(substitute(save(object, file = filename(object))))
})

## Contsimulation


setMethod("savedata", "Contsimulation", function(object,...){
  if(is.null(filename(object)))
     stop("This simulation has to be given a filename before it can be saved to harddisk")

  name <- deparse(substitute(object))

  eval.parent(parse(text=paste(name,"@Data <- NULL",sep="")))
  eval.parent(parse(text=paste(name,"@Data.id <- NULL",sep="")))
  eval.parent(parse(text=paste(name,"@Data.c <- NULL",sep="")))
  eval.parent(parse(text=paste(name,"@ind <- NULL",sep="")))
  eval.parent(substitute(save(object, file = filename(object))))
})
