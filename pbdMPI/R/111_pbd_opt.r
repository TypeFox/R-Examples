### Set pbd options.
pbd_opt <- function(..., bytext = "", envir = .GlobalEnv){
  if(!exists(".pbd_env", envir = envir)){
    envir$.pbd_env <- new.env()
  } 

  arg <- list(...)
  if(length(arg) > 0){
    names.arg <- names(arg)
    if(is.null(names.arg) || any(names.arg == "")){
      comm.stop("Options are all named.")
    }

    for(i.arg in 1:length(arg)){
      envir$.pbd_env[[names.arg[i.arg]]] <- arg[[i.arg]]
    }
  }

  if(bytext != ""){
    eval(parse(text = bytext), envir = envir$.pbd_env)
  }

  invisible()
} # End of pbd_opt().

