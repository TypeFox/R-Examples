#' Class aoos
#' 
#' This is an environment with some methods. Every class defined by \code{defineClass} will inherit from aoos. Summary will show a list of public and private members with approximated memory usage.
#' 
#' @exportClass aoos
#' @rdname aoos
setClass("aoos", contains = c("environment", "VIRTUAL"))

#' @rdname aoos
#' @export
setMethod("show", signature = c(object = "aoos"), 
          function(object) {
            cat("Class: ", class(object), "\n", sep = "")
            cat("public member:\n")
            lapply(ls(object), function(n) cat(" ", n, "\n"))
          })

#' @rdname aoos
#' @export
#' @param x object
#' @param name member name
setMethod("$", signature = c(x = "aoos"),
          function(x, name) {
            
            privacy <- !any(sapply(envirSearch(list(parent.frame())), 
                                   identical, y = parent.env(x)))
            
            getMember(name, x, privacy)
            
          })

envirSearch <- function(envList = list(environment())) {
  if(any(sapply(envList, identical, y = emptyenv()))) {
    envList
  } else {
    envirSearch(c(envList, list(parent.env(envList[[length(envList)]]))))
  }
}

getMember <- function(name, object, privacy = FALSE) {
  if(!privacy) {
    getPublicRepresentation(get(name, envir = parent.env(object)))
  } else {
    if(exists(name, envir = object, inherits = FALSE)) {
      getPublicRepresentation(get(name, envir = parent.env(object)))
    } else {
      stop(paste(name, "is not a public member."))
    }
  }
}

#' @rdname aoos
#' @export
#' @param value value to assign to. Will throw an error.
setMethod("$<-", signature = c(x = "aoos"),
          function(x, name, value) {
            
            privacy <- !any(sapply(envirSearch(list(parent.frame())), 
                                   identical, y = parent.env(x)))
            
            if (privacy) {
              stop("If you need to extend object, modify class definition.")
            } else {
              assign(name, value = value, envir = parent.env(x))
            }
            
            x
          })

#' @rdname aoos
#' @param object object
#' @param ... arguments passed to method (not used).
#' @export
summary.aoos <- function(object, ...) {
  out <- envSize(parent.env(object))
  rownames(out) <- NULL
  out
}

envSize <- function(env) {
  
  napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = env)))
  
  names <- ls(env, all.names = TRUE)
  
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)

  obj.subsizes <- napply(names, function(x) {
    if (inherits(x, "publicFunction") && !identical(environment(x), env))
      envSize(environment(x))
  })
  obj.subsizes <- do.call(rbind, obj.subsizes)
  obj.subsizes$Name <- rownames(obj.subsizes)
  
  out <- rbind(data.frame(
    Name = names(obj.type),
    Type = obj.type, 
    "Size.Mib" = round(obj.size / (1024^2), 1),
    stringsAsFactors = FALSE), 
    obj.subsizes)
  
  out[order(out$Name, out$Type), ]
    
}

#' @rdname aoos
#' @export
setMethod("as.environment", "aoos", function(x) parent.env(x))