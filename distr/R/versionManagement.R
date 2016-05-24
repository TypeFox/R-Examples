### elementary version management:
##  detects missing slots

setMethod("isOldVersion", "ANY",
# isold <-
function(object) {
   if(!isClass(class(object)))
      stop("'isOldVersion()' only works for formal S4-Classes.")
   slotNames <- slotNames(object)
   ers <- sapply(slotNames, function(x)
                is(try(slot(object, x), silent = TRUE), "try-error"))
   error <- any(ers)
   if(error)
      { warning(gettextf(
        paste("Object '%s' was defined under a deprecated version of class '%s'.",
              "Slot[s] %s is[are] missing.",
              "Try conv2NewVersion().", sep = "\n"),
              deparse(substitute(object)),
              class(object),
              paste("'",slotNames[ers],"'",sep="", collapse=", "))
               )
         return(TRUE)
   }else return(FALSE)     }
)

setMethod("conv2NewVersion", "ANY",
function(object) {
           tryobject <- new(class(object))
           slotNames <- slotNames(object)
           getIfExists <- function(x)
                    if (!is(try(slot(object, x), silent = TRUE),"try-error"))
                         slot(object,x) else slot(tryobject,x)
           lst <- sapply(slotNames, function(x) getIfExists(x))
           names(lst) <- slotNames
           #
           lst <- c(list(Class = class(object)), lst)
           if(is(try (myobj <- do.call("new", args = lst), silent = TRUE), "try-error")){
              myobj <- tryobject
              for(i in 2:length(lst)) 
                  slot(myobj, name=names(lst)[i]) <- lst[[i]]
           }           
           myobj
           })

setMethod("conv2NewVersion", "LatticeDistribution",
function(object) {
           slotNames <- slotNames(param(object))
           slotNames <- slotNames[slotNames!="name"]
           lst <- lapply(slotNames, function(x) slot(param(object),x))
           names(lst) <- slotNames
           lst <- c(list(Class = class(object)), lst)
           myobj <- do.call("new", args = lst)
           myobj
            })

