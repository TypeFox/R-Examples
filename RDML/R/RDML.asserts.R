# library(assertthat)
is.opt.string <- function(x) {
  if(is.null(x) || is.na(x)) return(TRUE)
  is.string(x)
}
on_failure(is.opt.string) <- function(call, env) {
   paste0(deparse(call$x), " is present but not a string")
}

is.opt.list <- function(x) {
  if(is.null(x)) return(TRUE)
  is.list(x)
}
on_failure(is.opt.list) <- function(call, env) {
  paste0(deparse(call$x), " is present but not a list")
}

is.opt.list.one.el <- function(x) {
  if(is.null(x)) return(TRUE)
  is.list(x) && length(x) == 1  
}
on_failure(is.opt.list.one.el) <- function(call, env) {
  paste0(deparse(call$x), " is present but not a list or length > 1")
}

is.opt.flag <- function(x) {
  if(is.null(x) || is.na(x)) return(TRUE)
  is.flag(x)
}
on_failure(is.opt.flag) <- function(call, env) {
  paste0(deparse(call$x), " is present but not a logical or length > 1")
}

is.opt.double <- function(x) {
  if(is.null(x) || is.na(x)) return(TRUE)
  is.double(x) && length(x) == 1  
}
on_failure(is.opt.double) <- function(call, env) {
  paste0(deparse(call$x), " is present but not a double or length > 1")
}

is.double.matrix <- function(x) {
  is.double(x) && is.matrix(x) 
}
on_failure(is.double.matrix) <- function(call, env) {
  paste0(deparse(call$x), " is not a 'double' matrix")
}

is.opt.double.matrix <- function(x) {
  if(is.null(x) || is.na(x)) return(TRUE)
  is.double(x) && is.matrix(x) 
}
on_failure(is.opt.double.matrix) <- function(call, env) {
  paste0(deparse(call$x), " is present but not a 'double' matrix")
}

is.opt.integer <- function(x) {
  if(is.null(x) || is.na(x)) return(TRUE)
  is.integer(x) && length(x) == 1  
}
on_failure(is.opt.integer) <- function(call, env) {
  paste0(deparse(call$x), " is present but not a integer or length > 1")
}

is.type <- function(x, test.type) {
  class(x)[1] == substitute(test.type)
}
on_failure(is.type) <- function(call, env) {
  sprintf("%s is not a '%s' or length > 1",
          deparse(call$x),
          deparse(call$test.type))
}

is.opt.type <- function(x, test.type) {
  if(is.null(x)) return(TRUE)
  class(x)[1] == substitute(test.type)
}
on_failure(is.opt.type) <- function(call, env) {
  sprintf("%s is present but not a '%s' or length > 1",
          deparse(call$x),
          deparse(call$test.type))
}

is.opt.list.type <- function(x, test.type) {
  if (is.null(x)) return(TRUE)
  if (!is.list(x)) return(FALSE)
  for(el in x) {
    if (class(el)[1] != substitute(test.type))
      return(FALSE)
  }
  TRUE
}
on_failure(is.opt.list.type) <- function(call, env) {
  sprintf("%s is present but not a '%s' list",
          deparse(call$x),
          deparse(call$test.type))
}

is.list.type <- function(x, test.type) {
  if (!is.list(x)) return(FALSE)
  for(el in x) {
    if (class(el)[1] != substitute(test.type))
      return(FALSE)
  }
  TRUE
}
on_failure(is.list.type) <- function(call, env) {
  sprintf("%s is not a '%s' list",
          deparse(call$x),
          deparse(call$test.type))
}

is.opt.count <- function(x) {
  if(is.null(x) || is.na(x)) return(TRUE)
  
}
on_failure(is.opt.integer) <- function(call, env) {
  paste0(deparse(call$x), " is present but not a count or length > 1")
}

is.float <- function(x) {
  is.numeric(x) && length(x) == 1
}
on_failure(is.float) <- function(call, env) {
  paste0(deparse(call$x), " is not a float (numeric) or length > 1")
}

is.enum <- function(x, levels) {
  if (is.null(x)) return(TRUE)
  x %in% c(levels, NA)
}
on_failure(is.enum) <- function(call, env) {
  sprintf("'%s' has level '%s' not in: %s",
          deparse(call$x),
          eval(call$x, envir = env),
          paste0(eval(call$levels, envir = env), collapse = ", "))
}

has.only.names <- function(x, only.names) {
  !(FALSE %in% (names(x) %in% only.names))
}
on_failure(has.only.names) <- function(call, env) {
  paste0(deparse(call$x), " has names not in: ", paste(call$only.names, collapse = ", "))
}

is.to.remove.list <- function(x) {  
  for(el in x[-which(names(x) == "id")]) { 
    if(!is.na(el)) return(FALSE)        
  }
  TRUE
}
on_failure(is.to.remove.list) <- function(call, env) {
  paste0(deparse(call$x), " is not a command to remove from list")
}