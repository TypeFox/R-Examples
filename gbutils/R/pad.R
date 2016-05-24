setClass("objectPad",
         representation("list")
         # , prototype = prototype(list())
         )

setMethod("initialize",
          "objectPad",
          function(.Object,...) {
            # callNextMethod mozhe da e neobhodim za naslednitsi na tozi class.
            # .Object <- callNextMethod()            # dali e na pravilnoto myasto?
            sklad <- list()
            get <- function(x){
              if(missing(x))
                sklad
              else
                sklad[[x]]
            }
            set <- function(x,value){
                sklad[[x]] <<- value
            }
            .Object$get <- get
            .Object$set <- set
            .Object
          }
          )

setGeneric("pad", def = function(x,item){standardGeneric("pad")},
           useAsDefault = function(x,item){ pad(x@pad,item)  } # mozhe s poveche proverki!
           )

setGeneric("pad<-", def = function(x,item,...,value){standardGeneric("pad<-")},
           useAsDefault = function(x,item,...,value){
             if(length(list(...))>0)
               pad(x@pad,item,...) <- value
             else
               pad(x@pad,item) <- value
             x    # the above expects that the required values are changed internally.
                  # otherwise assignments to x and/or x@pad will be necessary.
           }
           )

setMethod("pad", signature(x = "objectPad", item = "missing")
          , function(x){ x$get() }
          )

setMethod("pad", signature(x = "objectPad", item = "ANY")
          , function(x,item){ x$get(item) }
          )

setMethod("pad", signature(x = "ANY", item = "missing")
          , function(x){ x@pad$get() }
          )


setReplaceMethod("pad", signature(x = "objectPad", item = "ANY", value="ANY")
          , function(x,item,value){ x$set(item,value); x }
          )



padcheck <- function(x, item){   # check if item is set, very crude.
  wrk <- pad(x,item)
  !is.null(wrk)
}
