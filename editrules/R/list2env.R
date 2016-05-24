# legacy function to keep running under R<=2.10
if (!exists("list2env")){
  list2env <- function(x, envir=NULL, parent=parent.frame()){
     if (is.null(envir))
        envir <- new.env(parent=parent)
     for (i in names(x)){
        envir[[i]] <- x[[i]]
     }
     envir
  }
}



