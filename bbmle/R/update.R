## setGeneric("update", function(object, formula.,  ..., evaluate=TRUE)
##           standardGeneric("update"))


## FIXME: compare these two
## setMethod("update", "mle2",
## function (object, ..., evaluate = TRUE)
## {
##   call <- object@call
##   extras <- match.call(expand.dots = FALSE)$...
##   if (length(extras) > 0) {
##     existing <- !is.na(match(names(extras), names(call)))
##     for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
##     if (any(!existing)) {
##       call <- c(as.list(call), extras[!existing])
##       call <- as.call(call)
##     }
##   }
##   if (evaluate) eval(call, parent.frame()) else call
## })

 

## update.default, modified with $ turned to @ as appropriate
setMethod("update", "mle2",
          function (object, formula., evaluate = TRUE, ...) 
          {
            call <- object@call
            extras <- match.call(expand.dots = FALSE)$...
            if (!missing(formula.)) 
              call$minuslogl <- update.formula(formula(object), formula.)
            if (length(extras)) {
              existing <- !is.na(match(names(extras), names(call)))
              for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
              if (any(!existing)) {
                call <- c(as.list(call), extras[!existing])
                call <- as.call(call)
              }
            }
            if (evaluate) 
              eval(call, parent.frame())
            else call
          })
