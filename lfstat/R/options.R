setlfunit <- function(string = ""){
  if(!is.character(string)){
    stop('"string" must be a character string"')
  }
  options("lfstat" = list(unit = string))
}

.onLoad <- function(libname, pkgname) {
  setlfunit()
}


lflabel <- function(x = "Flow"){
  if(getOption("lfstat")$unit == ""){
    return(x)
  } else {
    title <-  paste0('"', x, '"~', getOption("lfstat")$unit)
  }
  parse(text = title)
}
