### simple version management

## use "isOldVersion"-"Any" method from distr now...

#setMethod("isOldVersion", "Dataclass", function(object) {
#   if(is(try(slot(object, "version"), silent = TRUE), "try-error"))
#      { warning(gettextf(
#"Object '%s' was defined under a deprecated version of class 'DataClass' (or subclasses). Try conv2NewVersion()...",
#        deparse(substitute(object))
#                         )
#               )
#         return(TRUE)
#   }else return(FALSE)     })


setMethod("conv2NewVersion", "Dataclass", function(object) {
           myobj <- new("Dataclass",
                         Data = aperm(array(object@Data,
                                   c(object@runs, object@samplesize, 1)
                                           ), perm = c(2,3,1)
                                     ),
                         filename = object@filename
                         )
           myobj
           })

