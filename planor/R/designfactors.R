#---------------------------------------------------------------------------
# CLASS designfactors and its METHODS
#---------------------------------------------------------------------------
#  Methods of "designfactors" : names, length, "[", bind
#---------------------------------------------------------------------------
setClass("designfactors",
         representation(fact.info="data.frame",
                        pseudo.info="data.frame",
                        levels="list"))

# --------------------------------------
# "names" method
# --------------------------------------
setMethod("names", signature(x="designfactors"),
          definition=function(x){rownames(x@fact.info)})

# --------------------------------------
# "length" method
# --------------------------------------
setMethod("length", signature(x="designfactors"),
          definition=function(x){nrow(x@fact.info)})

# --------------------------------------
# "[" method
# --------------------------------------

# --------------------------------------
setMethod("[",
signature(x = "designfactors", i = "ANY", j = "ANY", drop = "ANY"),
          definition=function(x,i,j,...,drop){
            x@fact.info <- x@fact.info[i,]
            x@levels <- x@levels[i]
            x <- planor.pseudofactors(x)
            x
          })

#-------------------------------------------------------
# "bind.designfactors"
# --------------------------------------
# ARGUMENTS
# - x: an object of class designfactors
# - y: an object of class designfactors
# RETURN
#   An object of class designfactors,  where the input are binded.
# --------------------------------------
bind.designfactors <- function(x,y){
  # Warning if factors have the same name
  a <- dimnames(x@fact.info)[[1]]
  b <- dimnames(y@fact.info)[[1]]
  if (any(a%in%b) | any(b%in%a))
    warning("Binding of factors with the same name\n")

  z <- new("designfactors",
             fact.info = rbind(x@fact.info,y@fact.info),
             levels = c(x@levels,y@levels))
  # to avoid names duplication
  names(z@levels) <- rownames(z@fact.info)
  # pseudofactors information
  z <- planor.pseudofactors(z)
  z
}

# --------------------------------------
# "bind" method 
# --------------------------------------
setMethod("bind", signature(x="designfactors",y="designfactors"),
          definition=bind.designfactors)
