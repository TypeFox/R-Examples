## Slots in Evaluation

## name        - name of Dataclass object, which was called by evaluate
## filename    - filename of this object
## call.ev     - call which created the object, e.g.evalate(contsim, mean)
## result      - result of estimation on data
## estimator   - used estimation function

setMethod("name", "Evaluation", function(object) object@name)
setMethod("Data", "Evaluation", function(object) object@Data)
setMethod("filename", "Evaluation", function(object) object@filename)
setMethod("call.ev", "Evaluation", function(object) object@call.ev)
setMethod("result", "Evaluation", function(object) object@result)
setMethod("estimator", "Evaluation", function(object) object@estimator)

setReplaceMethod("name", "Evaluation", 
   function(object, value){ object@name <- value; object}) ### new 1.8
setReplaceMethod("filename", "Evaluation", 
   function(object, value){ object@filename <- value; object}) ### new 1.8



## Save method
#if(!isGeneric("savedata")) setGeneric("savedata", 
# function(object) standardGeneric("savedata"))

setMethod("savedata", "Evaluation", 
          function(object, estimatorName = NULL, fileN = NULL, ...){
            
            name0 <- as.character(match.call(call=sys.call(-1))$object)
            if(is.null(estimatorName)) 
               (estimatorName <- as.character(call.ev(object)$estimator))
            if(is.null(fileN)) 
               {if(is.null(filename(object))) 
stop("This Dataclass object has to be given a filename before it can be saved to harddisk")
                fileN <- paste(filename(object), ".", estimatorName, sep = "")}
            
#            save(object, file = fileN)
            eval.parent(parse(text = paste("save(", name0,", file = \"",
                                            fileN,"\")", sep = "")))
 
            namecomment <- paste(name0, ".comment", sep = "")  
            commentfile <- paste(fileN, ".comment", sep = "")
            
            eval.parent(parse(text = paste(namecomment," <- ",name0, sep = ""))) 
            eval.parent(parse(text = paste(namecomment,"@result <- NULL", 
                                           sep = ""))) 
            eval.parent(parse(text = paste("save(",namecomment,", file = \"",
                                         commentfile, "\")", sep = ""))) 
            eval.parent(parse(text = paste("rm(",namecomment,")", sep = "")))
            print(fileN)
            print(name0)
            print(commentfile)
            print(namecomment)
          })

