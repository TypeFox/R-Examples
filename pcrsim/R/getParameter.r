################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 03.12.2015: getParameter(what="Method") now returns all methods for all kits.
# 10.11.2015: Return only correct 'method' if what=NA.
# 21.10.2015: Implemented 'method'.
# 09.09.2014: Re-written new structure.
# 19.05.2013: Re-written for reading data from text file.

#' @title Get Kit Parameters
#'
#' @description
#' Provides parameters for simulation for different STR kits and methods.
#'
#' @details
#' The function returns various information for kit and parameters specified
#' in parameters.txt.
#' 
#' @param kit string or integer specifying the kit.
#' @param what string to specify which information to return. Default is 'NA' which return all info.
#'  Not case sensitive.
#' @param method string to specify which method to return. Default is 'NA' which return all info.
#'  Not case sensitive.
#' @param show.messages logical, default TRUE for printing messages to the R promt.
#' @param .kit.param data frame, run function on a data frame instead of the parameters.txt file.
#' @param debug logical indicating printing debug information.
#' 
#' @return vector of data frame with kit information.
#' 
#' @importFrom data.table data.table
#' @importFrom utils head read.delim
#' 
#' @export 
#' @examples
#' 
#' # Returns vector of available kits.
#' getParameter()
#' 
#' # Returns vector of all methods.
#' getParameter(what="methods")
#' 
#' # Returns methods for specified kit.
#' getParameter(kit="SGMPlus", what="methods")
#' 
#' # Returns vector of available options.
#' getParameter(what="options")
#' 
#' # Returns vector of markers for specified kit.
#' getParameter(kit="SGMPlus", what="Marker")
#' 
#' # Returns data frame of all information for specified kit and method.
#' getParameter(kit="SGMPlus", method = "Default")

getParameter<-function(kit=NULL, what=NA, method=NA, show.messages=FALSE,
                       .kit.param=NULL, debug=FALSE) {
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("kit:")
    print(kit)
    print("what:")
    print(what)
    print("method:")
    print(method)
    print("show.messages:")
    print(show.messages)
    print(".kit.param:")
    print(head(.kit.param))
  }

  # CHECK PARAMETERS  #########################################################
  
  if(!is.logical(show.messages)){
    stop(paste("'show.messages' must be logical!"))
  }
  
  if(length(method) == 0){
    method <- NA
  }
  if(length(what) == 0){
    what <- NA
  }
  
  # PREPARE  ##################################################################
  
  .separator <- .Platform$file.sep # Platform dependent path separator.

  # Define options for 'what'.
  options <- c("Method", "Marker", "Efficiency", "Stutter",
               "Threshold", "Scaling", "Options")
  
  # Convert to upper case.  
  what <- toupper(what)
  method <- toupper(method)
  
  # LOAD KIT INFO  ############################################################ 
  
  if(is.null(.kit.param)){
    # Get package path.
    packagePath <- path.package("pcrsim", quiet = FALSE)
    subFolder <- "extdata"
    fileName <- "parameters.txt"
    
    filePath <- paste(packagePath, subFolder, fileName, sep=.separator)
    
    .kit.param <- read.delim(file=filePath, header = TRUE, sep = "\t", quote = "\"",
                           dec = ".", fill = TRUE, stringsAsFactors=FALSE)
    
  }
  
  # Available kits. Must match else if construct.
  kitShorNames <- unique(.kit.param$Short.Name)
  
  # Check if NULL
  if (is.null(kit)) {
    
    if(!is.na(what)){

      if(grepl("METHOD", what)){
        # Use grepl to also match 'Methods'.
        
        res <- unique(.kit.param$Method)
        
        # Print available kits
        if (show.messages){
          print("All methods:")
          print(res)
        }
        
      } else if(grepl("OPTION", what)){
        # Use grepl to also match 'Options'.
        
        # Return options.
        res <- options
        
      } else if (what == "ID"){
        
        # Return id's.
        res <- unique(paste(.kit.param$Short.Name, ":", .kit.param$Method, sep=""))
        
      }
      
    } else {

      res <- kitShorNames
      
      # Print available kits
      if (show.messages){
        print("Available kits:")
        print(kitShorNames)
      }
      
    }  
      
  } else { # String provided.
    
    # Check if number or string.
    if (is.numeric(kit)) {
      
      # Set index to number.
      index <- kit
      
    } else {
      
      # Find matching kit index (case insensitive)
      index <- match(toupper(kit), toupper(kitShorNames))
      
    }
    
    # No matching kit.
    if (is.na(index)) {
      
      # Print available kits
      if (show.messages){
        print("No matching kit!")
        print("Available kits:")
        print(kitShorNames)
      }
      res <- NA
      
    } else { # Assign matching kit information.
      
      # Get matching kit.
      res <- .kit.param[.kit.param$Short.Name == kitShorNames[index], ]
      
    } 
    
  }
  
  if(debug){
    print(str(res))
    print(head(res))
    print(tail(res))
  }

  # If kit is not null.
  if (!is.null(kit)) {
    
    if (is.na(what)) {

      if(!is.na(method)){

        # Filter out method specific parameters.
        res <- res[toupper(res$Method) == method, ]

      }
      
      #return(res)

    } else if (grepl("METHOD", what)) {
      # Use grepl to also match 'Methods'.
        
      #return(unique(res$Method))
      res <- unique(res$Method)
      
    } else if (what == "ID") {
      
      res <- unique(paste(res$Short.Name, ":", res$Method, sep=""))
      
      #return(id)
      
    } else if (what == "MARKER") {
      
      res <- as.vector(unique(res$Marker))
      #return(as.vector(unique(res$Marker)))
      
    } else if (what == "EFFICIENCY") {
      
      if(length(unique(res$Method)) > 1){
        
        if(!is.na(method)){
          
          # Extract method parameters.
          res <- res[toupper(res$Method) == method, ]
          
        } else {
          
          stop(paste("There are multiple methods for the ", kitShorNames[index],
                     "kit, but no method is provided!"))
        }
        
      }
      
      # Extract information.
      res <- data.frame(Marker = res$Marker, PCR.Efficiency = res$PCR.Efficiency)

      #return(res)
      
    } else if (what == "STUTTER") {
      
      if(length(unique(res$Method)) > 1){
        
        if(!is.na(method)){
          
          # Extract method parameters.
          res <- res[toupper(res$Method) == method, ]
          
        } else {
          
          stop(paste("There are multiple methods for the ", kitShorNames[index],
                     "kit, but no method is provided!"))
        }
        
      }
      
      # Extract information.
      res <- data.frame(Marker = res$Marker, Stutter.Probability = res$Stutter.Probability, 
                        Stutter.Max = res$Stutter.Max, Stutter.Type.Repeat = res$Stutter.Type.Repeat,
                        Stutter.Type.Bp = res$Stutter.Type.Bp, Stutter.Intercept = res$Stutter.Intercept,
                        Stutter.Slope = res$Stutter.Slope, Stutter.Sigma = res$Stutter.Sigma,
                        stringsAsFactors = FALSE)
      
      #return(res)
      
    } else if (what == "THRESHOLD") {

      # Select columns and get unique rows by key.
      res <- data.table(res) # Convert to data.table...
      selected <- c("Method", "T.Intercept", "T.Slope", "T.Sigma")
      res <- res[, selected, with=FALSE]
      res <- unique(res, by=selected)

      if(length(unique(res$Method)) > 1){
        
        if(!is.na(method)){
          
          # Extract method parameters.
          res <- res[toupper(res$Method) == method, ]
          
        } else {
          
          stop(paste("There are multiple methods for the ", kitShorNames[index],
                     "kit, but no method is provided!"))
        }
        
      }
      
      res <- as.data.frame(res) # ...and back to avoid possible negative effects.
      names(res) <- c("Method", "Intercept", "Slope", "Sigma")

      #return(res)
      
    } else if (what == "SCALING") {
      
      # Select columns and get unique rows by key.
      selected <- c("Method", "Scaling.Intercept", "Scaling.Slope", "Scaling.Sigma")
      res <- data.table(res) # Convert to data.table...
      res <- res[, selected, with=FALSE]
      res <- unique(res, by=selected)
      
      if(length(unique(res$Method)) > 1){

        if(!is.na(method)){
          
          # Extract method parameters.
          res <- res[toupper(res$Method) == method, ]

        } else {
          
          stop(paste("There are multiple methods for the ", kitShorNames[index],
                     "kit, but no method is provided!"))
        }
        
      }

      res <- as.data.frame(res) # ...and back to avoid possible negative effects.
      names(res) <- c("Method", "Intercept", "Slope", "Sigma")
      
      #return(res)
      
    } else if (what == "OPTIONS") {
      
      # Return options.
      res <- options
      
    } else {
     
      warning(paste(what, "not supported! \nwhat = {", 
                    options, "}"))
      #return(NA)
      res <- NA
     
    }
    
  } # End if is.null(kit).
  
  # Return kit information (or NA)
  return(res)
  
}