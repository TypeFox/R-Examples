#' Accessor class
#' 
#' This is a virtual class to be contained in other class definitions. It overrides the default accessor \code{$} and is intended to be used with the aoos class system (\code{\link{defineClass}}). Inherit from this class if you want to access public fields in the same way you access lists.
#'
#' @exportClass Accessor
#' @rdname Accessor
setClass("Accessor", contains = "VIRTUAL")

#' @rdname Accessor
#' @export
#' @param x object
#' @param name member name
setMethod("$", signature = c(x = "Accessor"),
          function(x, name) {
            
            privacy <- !any(sapply(envirSearch(list(parent.frame())), 
                                   identical, y = parent.env(x)))
            
            member <- getMember(name, x, privacy)
            
            if(inherits(member, "publicValue")) {
              member()
            } else {
              member
            }
            
          })

#' @rdname Accessor
#' @export
#' @param value value to assign to.
setMethod("$<-", signature = c(x = "Accessor"),
          function(x, name, value) {
            
            privacy <- !any(sapply(envirSearch(list(parent.frame())), 
                                   identical, y = parent.env(x)))
            
            member <- getMember(name, x, privacy)
            
            if(inherits(member, "publicValue")) {
              member(value)
            } else {
              
              if(privacy) {
                stop("If you need to extend object, modify class definition.")
              } else {
                assign(name, value = value, envir = parent.env(x))
              }
            
            }
            
            x
          })