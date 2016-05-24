## accessor and replacement functions for simObj slots

setGeneric("parms", function(obj, ...) standardGeneric("parms"))
setGeneric("parms<-", function(obj, value) standardGeneric("parms<-"))

setGeneric("init", function(obj, ...) standardGeneric("init"))
setGeneric("init<-", function(obj, value) standardGeneric("init<-"))

setGeneric("times", function(obj, ...) standardGeneric("times"))
setGeneric("times<-", function(obj, value) standardGeneric("times<-"))

setGeneric("inputs", function(obj, ...) standardGeneric("inputs"))
setGeneric("inputs<-", function(obj, value) standardGeneric("inputs<-"))

setGeneric("equations", function(obj, ...) standardGeneric("equations"))
setGeneric("equations<-", function(obj, value) standardGeneric("equations<-"))

setGeneric("solver", function(obj, ...) standardGeneric("solver"))
setGeneric("solver<-", function(obj, value) standardGeneric("solver<-"))

setGeneric("main", function(obj, ...) standardGeneric("main"))
setGeneric("main<-", function(obj, value) standardGeneric("main<-"))

setGeneric("initfunc", function(obj, ...) standardGeneric("initfunc"))
setGeneric("initfunc<-", function(obj, value) standardGeneric("initfunc<-"))

setGeneric("observer", function(obj, ...) standardGeneric("observer"))
setGeneric("observer<-", function(obj, value) standardGeneric("observer<-"))

setGeneric("out", function(obj, ...) standardGeneric("out"))
setGeneric("out<-", function(obj, value) standardGeneric("out<-"))
## the out slot is readonly

setMethod("parms", "simObj",
  function(obj) {
    obj@parms
  }
)

setMethod("parms<-", "simObj", 
    function(obj, value) {
      obj@parms <- value
      obj@out <- NULL
      invisible(obj)
    }
)

setMethod("times", "simObj",
    function(obj) {obj@times}
)

setMethod("times<-", "simObj",
    function(obj, value) {
      if (isfromtoby(obj@times) & hasfromtoby(value)) {
        # modify named vector with from, to and by
        for (i in 1:length(value)) {
          nam <- names(value[i])
          if (nam %in% c("from", "to", "by")) {
            obj@times[nam] <- value[[i]]
          } else {
            print(paste("Warning: vector element ", nam , " ignored"))
          }
        }
      } else {
        # overwrite vector
        if (is.null(names(value)) | isfromtoby(value)) {
          obj@times <- value
        } else {
          print("Warning: Ignored. Invalid (or incomplete) names in right hand side of assignment.")
        }
      }
      obj@out <- NULL
      invisible(obj)
    }
)

setMethod("init", "simObj",
    function(obj) {
      obj@init
    }
)

setMethod("init<-", "simObj",
    function(obj, value) {
      obj@init <- value
      obj@out <- NULL
      invisible(obj)
    }
)

setMethod("init<-", c("gridModel", "matrix"),
    function(obj, value) {
      obj@init <- value
      obj@out <- NULL      
      invisible(obj)
    }
)

setMethod("init<-", c("gridModel", "ANY"),
    function(obj, value) {
      print("error: value must be a matrix")
      invisible(obj)
    }
)

setMethod("inputs", "simObj",
    function(obj) {obj@inputs}
)

setMethod("inputs<-", "simObj",
    function(obj, value) {
      obj@inputs <- value
      obj@out <- NULL      
      invisible(obj)
    }
)

setMethod("main", "simObj",
    function(obj) {obj@main}
)

setMethod("main<-", "simObj",
    function(obj, value) {
      obj@main <- value
      obj@out <- NULL      
      invisible(obj)
    }
)

setMethod("equations", "simObj",
    function(obj) {obj@equations}
)

setMethod("equations<-", "simObj",
    function(obj, value) {
      obj@equations <- value
      obj@out <- NULL      
      invisible(obj)
    }
)

setMethod("initfunc", "simObj",
    function(obj) {obj@initfunc}
)

setMethod("initfunc<-", "simObj",
    function(obj, value) {
      obj@initfunc <- value
      obj@out <- NULL      
      invisible(obj)
    }
)


setMethod("observer", "simObj",
    function(obj) {obj@observer}
)

setMethod("observer<-", "simObj",
    function(obj, value) {
      obj@observer <- value
      invisible(obj)
    }
)


setMethod("solver", "simObj",
    function(obj) {obj@solver}
)

setMethod("solver<-", "simObj",
    function(obj, value) {
      obj@solver <- value
      obj@out <- NULL      
      invisible(obj)
    }
)

setMethod("out", "simObj",
  # returns a list
  function(obj, last=FALSE) {
    o <- obj@out
    if (last) o[length(o)] else o
  }
)

setMethod("out<-", "simObj",
    function(obj, value) {
      if(!is.null(value)) stop("``out'' can only be set to NULL")
      obj@out <- NULL      
      invisible(obj)
    }
)

setMethod("out", "gridModel",
  # if last==TRUE: returns a matrix (identical to init)
  function(obj, last=FALSE) {
    o <- obj@out
    if (last) o[[length(o)]] else o
  }
)

setMethod("out", "odeModel",
  # if last==TRUE: returns a vector (similar to init, time as first value)
  function(obj, last=FALSE) {
    o <- obj@out
    if (is.matrix(o)) o <- as.data.frame(o)
    if (last) unlist(o[nrow(o), ]) else o
  }
)


