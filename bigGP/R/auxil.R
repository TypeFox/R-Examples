remoteRm <- function(..., list = character()){
  dots <- match.call(expand.dots = FALSE)$...
  if (length(dots) && !all(sapply(dots, function(x) is.symbol(x) || 
                                  is.character(x)))) 
    stop("remoteRm: ... must contain names or character strings.")
  names <- sapply(dots, as.character)
  if(length(names) == 0L)
    names <- character()
  list <- c(list, names)
  status <- mpi.remote.exec(localRm, list, ret = TRUE)
  if("try-error" %in% sapply(status, class))
    stop("remoteRm: error on slaves:\n", status)  
  invisible(NULL)
}

localRm <- function(list){
  status <- try( {
    if(!is.character(list))
      stop("localRm: list must be a character string.")
    rm(list = list, pos = .GlobalEnv)
  } )
  invisible(status)
}

remoteLs <- function(all.names = FALSE)
  mpi.remote.exec(ls, .GlobalEnv, all.names = all.names, ret = TRUE)
# hard to do this in generality as need to deal with passing an envir on slave, not master, as the thing to be ls'ed

  
