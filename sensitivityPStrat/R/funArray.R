.Data.list <- NULL

.funVectorFun <- function(...) {
  lapply(.Data.list, FUN=function(.estimatorFunction, ...) {
    if(is.null(.estimatorFunction) || !is.function(.estimatorFunction))
      return(logical())
    .estimatorFunction(...)
  }, ...)
}

makeFunVector <- function(data) {
  x <- .funVectorFun

  if(typeof(data) != 'list') {
    storage.mode(data) <- 'list'
  }
  
  envir <- new.env(parent=environment(fun=.funVectorFun))
  valid <- unlist(lapply(data, FUN=function(x) is.null(x) || is.function(x) || is.na(x)))
  
  if(!all(valid)) {
    stop("some elements from data are not functions")
  }

  envir$.Data.list <- data

  environment(x) <- envir
  if(is.array(envir$.Data.list))
    class(x) <- c('funArray', 'funVector')
  else
    class(x) <- 'funVector'

  x
}

funVector <- function(length = 0) {
  x <- .funVectorFun

  envir <- new.env(parent=environment(fun=.funVectorFun))
  
  envir$.Data.list <- lapply(integer(length),
                             FUN=function(...) function() NULL)
  environment(x) <- envir
  class(x) <- 'funVector'

  x
}

funMatrix <- function(...) {
  fcall <- match.call()

  fcall[[1]] <- as.name('matrix')

  makeFunVector(eval.parent(fcall))
}

funArray <- function(...) {
  fcall <- match.call()

  fcall[[1]] <- as.name('array')

  makeFunVector(eval.parent(fcall))
}

'[.funVector' <- function(x, ..., drop=TRUE) {
  cl <- oldClass(x)
  oldX <- x
  oldEnv <- environment(x)

  x <- oldEnv$.Data.list
  
  x <- NextMethod(.Generic)

  if(length(x) <= 1)
    return(x[[1]])

  newEnvir <- new.env(parent=parent.env(oldEnv))
  
  newEnvir$.Data.list <- x
  environment(oldX) <- newEnvir
  oldX
}

'[.funArray' <- function(x, ..., drop=TRUE) {
  y <- NextMethod('[')

  cl <- oldClass(x)
  
  if(is.array(y))
    class(x) <- cl
  else if(is.function(y))
    return(y)
  else if(length(y) <= 1)
    return(y[[1]])
  else
    class(x) <- cl[!'funArray' %in% cl]

  newEnvir <- new.env(parent=parent.env(environment(x)))
  newEnvir$.Data.list <- y
  environment(x) <- newEnvir

  x
}

## 'str.funVector' <- function(x, ...) {
##   P0 <- function(...) paste(..., sep = "")                                                                            
##   mod <- "function"

##   if (is.array(object)) {
##     rnk <- length(di. <- dim(object))
##     di <- P0(ifelse(di. > 1, "1:", ""), di., ifelse(di. >
##                                                     0, "", " "))
##     pDi <- function(...) paste(c("[", ..., "]"),
##                                collapse = "")
##     le.str <- (if (rnk == 1)
##                pDi(di[1L], "(1d)")
##     else pDi(P0(di[-rnk], ", "), di[rnk]))
##     std.attr <- "dim"
##   }
##   else if (!is.null(names(object))) {
##     mod <- paste("Named", mod)
##     std.attr <- std.attr[std.attr != "names"]
##   }
##   if (has.class && length(cl) == 1) {
##     if (cl != mod && substr(cl, 1, nchar(mod)) !=
##         mod)
##       mod <- P0("'", cl, "' ", mod)
##     std.attr <- c(std.attr, "class")
##   }
##   str1 <- if (le == 1 && !is.array(object))
##     paste(NULL, mod)
##   else P0(" ", mod, if (le > 0)
##           " ", le.str)


## }

'[<-.funVector' <- function(x, ..., value) {
  if (!as.logical(length(value)))
    return(x)

  if(is.list(value)) {
    valid <- unlist(lapply(value, FUN=is.function))
  
    if(!all(valid)) {
      stop("some elements from value are not assignable to funVectors")
    }
  } else if(inherits(value, c('funVector'))) {
    value <- environment(value)$.Data.list
  } else if(is.function(value)) {
    value <- list(value)
  } else {
    stop("type of value cannot be assigned to a funVector")
  }
    
  oldX <- x
  oldEnv <- environment(x)
  
  x <- oldEnv$.Data.list

  x <- NextMethod('[<-')

  oldEnv$.Data.list <- x

  oldX
}

'[<-.funArray' <- function(x, ..., value) {
  NextMethod('[<-')
}

c.funVector <- function(..., recursive=FALSE)
  makeFunVector(c(unlist(lapply(list(...), function(x)
                                if(inherits(x, 'funVector'))
                                   environment(x)$.Data.list
                                else if(is.function(x) || is.list(x))
                                   x
                                else 
                                   function() NULL))))


cbind.funVector <- function (..., deparse.level = 1) {
  dotargs <- list(...)

  newargs <- lapply(dotargs, function(x) if(inherits(x, 'funVector')) environment(x)$.Data.list else x)
                    
  names(newargs) <- names(dotargs)

  makeFunVector(do.call('cbind', newargs))
}


rbind.funVector <- function (..., deparse.level = 1) {
  dotargs <- list(...)

  newargs <- lapply(dotargs, function(x) if(inherits(x, 'funVector')) environment(x)$.Data.list else x)

  names(newargs) <- names(dotargs)
#  str(do.call('rbind', newargs))

  makeFunVector(do.call('rbind', newargs))
}

names.funVector <- function (x) names(environment(x)$.Data.list)

'names<-.funVector' <- function(x, value) {
  names(environment(x)$.Data.list) <- value
  return(x)
}

dimnames.funVector <- function(x) dimnames(environment(x)$.Data.list)

'dimnames<-.funVector' <- function(x, value) {
  dimnames(environment(x)$.Data.list) <- value
  return(x)
}
