#' Set Macaulay2 Path
#'
#' This function sets the Macaulay2 path either by (1) passing it a character string or (2) using file.choose. 
#' 
#' @param path a character string, the path to m2
#' @return invisible m2Path
#' @export setM2Path
#' @examples
#' \dontrun{
#'
#' setM2Path() 
#'
#' }
#'
setM2Path <- function(path){

  if(missing(path) && interactive()){
  	
    m2Path <- dirname(file.choose())
    if(.Platform$OS.type == "windows" && str_detect(m2Path,"C:/")){
      m2Path <- str_replace(m2Path, "C:/", "/cygdrive/c/")
    }
    options(m2Path = m2Path)
    return(invisible(m2Path))
    
  } else if(!missing(path)){
  	
    options(m2Path = path)
    return(invisible(path))
    
  } else {
    stop(
      "If the session is not interactive, a path must be specified.", 
      call. = FALSE
    )
  }
}







#' Set Bertini Path
#'
#' This function sets the Bertini path either by (1) passing it a character string or (2) using file.choose. 
#' 
#' @param path a character string, the path to Bertini
#' @return invisible bertiniPath
#' @export setBertiniPath
#' @examples
#' \dontrun{
#'
#' setBertiniPath() 
#'
#' }
#'
setBertiniPath <- function(path){

  if(missing(path) && interactive()){
  	
    bertiniPath <- dirname(file.choose())
    if(.Platform$OS.type == "windows" && str_detect(bertiniPath,"C:/")){
      bertiniPath <- str_replace(bertiniPath, "C:/", "/cygdrive/c/")
    }    
    options(bertiniPath = bertiniPath)
    return(invisible(bertiniPath))
    
  } else if(!missing(path)){
  	
    options(bertiniPath = path)
    return(invisible(path))    
    
  } else {
    stop(
      "If the session is not interactive, a path must be specified.", 
      call. = FALSE
    )
  }
}











#' Set LattE Path
#'
#' This function sets the LattE path either by (1) passing it a character string or (2) using file.choose. 
#' 
#' @param path a character string, the path to LattE (the function count, for example)
#' @return invisible lattePath
#' @export setLattePath
#' @examples
#' \dontrun{
#'
#' setLattePath() 
#'
#' }
#'
setLattePath <- function(path){

  if(missing(path) && interactive()){
  	
    lattePath <- dirname(file.choose())
    if(.Platform$OS.type == "windows" && str_detect(lattePath,"C:/")){
      lattePath <- str_replace(dirname(lattePath), "C:/", "/cygdrive/c/")
    }    
    options(lattePath = lattePath)
    return(invisible(lattePath))
    
  } else if(!missing(path)){
  	
    options(lattePath = path)
    return(invisible(path))    
    
  } else {
    stop(
      "If the session is not interactive, a path must be specified.", 
      call. = FALSE
    )
  }
}












#' Set 4ti2 Path
#'
#' This function sets the 4ti2 path either by (1) passing it a character string or (2) using file.choose. 
#' 
#' @param path a character string, the path to 4ti2 (the function markov, for example)
#' @return invisible markovPath
#' @export setMarkovPath
#' @examples
#' \dontrun{
#'
#' setMarkovPath() 
#'
#' }
#'
setMarkovPath <- function(path){

  if(missing(path) && interactive()){
  	
    markovPath <- dirname(file.choose())
    if(.Platform$OS.type == "windows" && str_detect(markovPath,"C:/")){
      markovPath <- str_replace(markovPath, "C:/", "/cygdrive/c/")
    }        
    options(markovPath = markovPath)
    return(invisible(markovPath))
    
  } else if(!missing(path)){
  	
    options(markovPath = path)
    return(invisible(path))    
    
  } else {
    stop(
      "If the session is not interactive, a path must be specified.", 
      call. = FALSE
    )
  }
}
