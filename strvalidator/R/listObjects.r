################################################################################
# TODO LIST
# TODO: Option to return other info e.g. size. Return dataframe instead of vector.
#       Workaround: unname(sapply(listObjects(), function(x) object.size(get(x, envir = baseenv()))))

################################################################################
# CHANGE LOG (last 20 changes)
# 28.06.2015: Changed parameter names to format: lower.case
# 26.07.2013: 'obj.class' can now be a vector.
# 17.05.2013: New parameters 'obj.class', 'debug'.
# 17.05.2013: Made general. Changed name from listDataFrames -> listObjects.
# <17.05.2013: First version.

#' @title List Objects
#'
#' @description
#' Internal helper function to list objects in an environment.
#'
#' @details
#' Internal helper function to retrieve a list of objects from a workspace.
#' Take an environment as argument and optionally an object class.
#' Returns a list of objects of the specified class in the environment.
#' 
#' @param env environment in wich to search for objects.
#' @param obj.class character string or vector specifying the object class.
#' @param debug logical indicating printing debug information.
#' 
#' @return character vector with the object names.
#' 
#' @export
#'  
#' @examples
#' \dontrun{
#' # List data frames in the workspace.
#' listObjects(obj.class="data.frame")
#' # List functions in the workspace.
#' listObjects(obj.class="function")
#' }

listObjects <- function(env=parent.frame(), obj.class=NULL, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }

  # Result vector.
  res <- character()
  
  # List objects in environment.
  wsObj <- ls(env)

  if(debug){
    print("Objects:")
    print(wsObj)
  }
  
  # Check if specified object class.
  if(!is.null(obj.class)){

    classes <- list()

    # Loop to save all class information.
    for(i in seq(along=wsObj)){
      obj <- get(wsObj[i], envir=env)
      classes[i] <- list(class(obj))
      
    }
    
    # Filter objects with specified classes.
    for(c in seq(along=obj.class)){
      for(i in seq(along=classes)){
        if(obj.class[c] %in% classes[[i]]){
          res <- c(res, wsObj[i])
        }
      }
    }
    
  } else {

    # Return all objects.
    res <- wsObj
    
  }

  if(debug){
    print("Returned objects:")
    print(res)
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(res)
  
}
