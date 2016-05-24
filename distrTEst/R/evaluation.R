

setMethod("evaluate", signature(object = "Dataclass", estimator = "function"),
              function(object, estimator, ..., resname = NULL,
                   name = NULL, filename = NULL){
  if(is.null(Data(object)))
    stop("No Data found -> simulate first")
  ## is it a Simulation?
  object0 <- object
  if(!is(try(slot(object,"rate"),silent=TRUE),"try-error"))
     {object0@Data <- NULL}
  if(is.null(filename(object)))
    stop("filename of Dataclass object is NULL, hence results couldn't be saved to harddisk")
  if(filename(object) == "")
    stop("filename of Dataclass object is an empty string, hence results couldn't be saved to harddisk")
  Data <- Data(object)
  if (is.null(name)) name <- as.character(substitute(object))
  if (is.null(filename)) filename <- filename(object)
  if (is.null(resname) && !is(object,"Simulation")) resname <- "res"
  if (is.null(resname)) resname <- as.character(match.call(call=sys.call(-1))$estimator)
  ### new 04-10-06
  result <- .convert.result.format(apply(X = Data, MARGIN = 3,
                                         FUN = estimator, ...),
                                   resname = resname)
  # old
  # old  if(!is.vector(Data)) result <- apply(X = Data, MARGIN = 1,
  # old                                       FUN = estimator, ...)
  # old     else result <- apply(X = as.matrix(Data), MARGIN = 2,
  # old                          FUN = estimator, ...)
  new("Evaluation", name = name, call.ev = match.call(call=sys.call(-1)),
       filename = filename, result = result, estimator = estimator,
       Data = object0)
})




### new 04-10-06
setMethod("evaluate", signature(object = "Contsimulation",
           estimator = "function"),
          function(object, estimator, ..., resname = NULL,
                   name = NULL, filename = NULL){
  if(is.null(Data.id(object)))
    stop("No Data found -> simulate first")
  ## is it a Simulation?
  object0 <- object
  if(!is(try(slot(object,"rate"),silent=TRUE),"try-error"))
     {object0@Data.id <- NULL
      object0@Data.c <- NULL
      object0@Data <- NULL
      object0@ind <- NULL}
  if(is.null(filename(object)))
    stop("filename of Dataclass object is NULL, hence results couldn't be saved to harddisk")
  if(filename(object) == "")
    stop("filename of Dataclass object is an empty string, hence results couldn't be saved to harddisk")
  if (is.null(name)) name <- as.character(substitute(object))
  if (is.null(filename)) filename <- filename(object)
  if (is.null(resname))
      resname <- as.character(match.call(call=sys.call(-1))$estimator)

  Data.id <- Data.id(object);
  Data.re <- Data(object)

  result.id <- .convert.result.format(apply(X = Data.id, MARGIN = 3,
                                            FUN = estimator, ...),
                                      resname = resname)
  result.re <- .convert.result.format(apply(X = Data.re, MARGIN = 3,
                                            FUN = estimator, ...),
                                      resname = resname)

  names(result.id) <- paste(names(result.id),"id",sep=".")
  names(result.re) <- paste(names(result.re),"re",sep=".")
  new("Evaluation", name = name, call.ev = match.call(call=sys.call(-1)),
       filename = filename,
       result = as.data.frame(cbind(result.id, result.re)),
       estimator =  estimator, Data = object0)
})
