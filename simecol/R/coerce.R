setAs("simObj", "list",
  function(from) {
    slots <- slotNames(from)
    ll <- list()
    for (s in slots) {
      ll[[s]] <- slot(from, s)
    }
    ll$class <- class(from)
    ll
  }
)

## as.list is already a generic
setMethod("as.list", "simObj",
  function(x, ...) as(x, "list")
)

setAs("list", "simObj",
  function(from) {
    slots <- names(from)
    cl <- from$class[[1]]
    if (is.null(cl))
      stop("list has no 'class' element")
    if (!extends(cl, "simObj"))
      stop(paste("class '", cl, "' is no valid simecol class", sep=""))
    slots <- slots[which(slots!="class")]
    obj <- new(cl)
    for (s in slots) {
      slot(obj, s) <- from[[s]]
    }
    initialize(obj)
  }
)

setGeneric("as.simObj", function(x, ...) standardGeneric("as.simObj"))

setMethod("as.simObj", "list",
  function(x, ...) as(x, "simObj")
)


