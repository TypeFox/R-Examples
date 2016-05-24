################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 15.12.2014: Changed parameter names to format: lower.case
# 24.07.2014: Fixed 'NA' bug when recycling names.
# 26.07.2013: Added 'any' to class if-case.
# 11.07.2013: Fixed scope bug.
# 10.07.2013: First version.


#' @title Export
#'
#' @description
#' Exports or saves various objects.
#'
#' @details
#' Export objects to a directory on the file system. Currently only objects
#' of class data.frames or ggplot are supported. data.frame objects will
#' be exported as '.txt' and ggplot objects as '.png'.
#' .RData applies to all supported object types.
#' 
#' @param object string, list or vector containing object names to be exported. 
#' @param name string, list or vector containing file names. 
#' Multiple names as string must be separated by pipe '|'.
#' If not equal number of names as objects, first name will be used to construct names. 
#' @param use.object.name logical, if TRUE file name will be the same as object name.
#' @param env environment where the objects exists.
#' @param path string specifying the destination folder exported objects.
#' @param ext string specifying file extension.
#' Default is 'auto' for automatic .txt or .png based on object class.
#' If .RData all objects will be exported as .RData files.
#' @param delim string specifying the delimiter used as separator.
#' @param width integer specifying the width of the image.
#' @param height integer specifying the height of the image.
#' @param res integer specifying the resolution of the image.
#' @param overwrite logical, TRUE if existing files should be overwritten.
#' @param debug logical indicating printing debug information.
#' 
#' @return NA if all objects were exported OR,
#' data.frame with columns 'Object', 'Name', and 'New.Name' with objects
#' that were not exported.
#' 
#' @importFrom utils write.table
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot
#' 
#' @keywords internal
#' 

export <- function(object, name=NA, use.object.name=is.na(name),
                   env=parent.frame(), path=NA, 
                   ext="auto", delim="\t", 
                   width=3000, height=2000, res=250,
                   overwrite=FALSE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("object")
    print(object)
    print("name")
    print(name)
    print("use.object.name")
    print(use.object.name)
    print("env")
    print(environmentName(env))
    print("path")
    print(path)
    print("ext")
    print(ext)
    print("delim")
    print(delim)
    print("width")
    print(width)
    print("height")
    print(height)
    print("res")
    print(res)
    print("overwrite")
    print(overwrite)
  }
  
  .uniqueFile <- function(object, name, path){
    # Returns TRUE if filname is unique and FALSE if file exist.

    if(debug){
      print("failed")
      print(failed)
    }
    
    # Check if file exist.
    if(file.exists(path)){
      
      # Make a data frame of objects that could not be exported.
      failed <<- rbind(failed, 
                       data.frame(Object=object, Name=name, New.Name=as.character(NA),
                                  stringsAsFactors=FALSE))
      
      if(debug){
        print(paste("file '", name, "' already exist!", sep=""))
      }
      
      # Set flag.
      available <- FALSE
      
    } else {
      
      # Set flag.
      available <- TRUE

      if(debug){
        print(paste("file '", name, "' is available!", sep=""))
      }
    }
    
    return(available)
    
  }
  
  # Constants.
  .separator <- .Platform$file.sep # Platform dependent path separator.
  
  # Variables.
  # Create empty result data frame with NAs.
  failed <- data.frame(t(rep(NA,3)), stringsAsFactors=FALSE)
  # Add column names.
  names(failed) <- c("Object","Name", "New.Name")
  # Remove all NAs
  failed <- failed[-1,]
  
  # CHECK DATA ----------------------------------------------------------------
  
  if(debug){
    print("CHECK DATA")
  }
  
  if(!all(is.na(name))){
    if(is.list(name)){
      # Convert list to vector.
      name <- unlist(name)
    }
    if(all(name == "")){
      # Replace empty string with 'NA'.
      name <- NA
    }
  } else if (!use.object.name){
    stop("'name' is required",
         call. = TRUE)
  }
  
  if(is.na(path)){
    stop("'path' is required",
         call. = TRUE)
  }
  
  if(!is.logical(use.object.name)){
    stop("'use.object.name' must be logical",
         call. = TRUE)
  }
  
  if(!is.logical(overwrite)){
    stop("'overwrite' must be logical",
         call. = TRUE)
  }
  
  if(width < 20 || !is.numeric(width) || is.na(width)){
    stop("'width' must a numeric <20",
         call. = TRUE)
  }
  
  if(height < 20 || !is.numeric(height) || is.na(height)){
    stop("'height' must a numeric <20",
         call. = TRUE)
  }
  
  if(res < 20 || !is.numeric(res) || is.na(res)){
    stop("'res' must a numeric <20",
         call. = TRUE)
  }
  
  # PREPARE ----------------------------------------------------------------
  
  if(debug){
    print("PREPARE")
  }
  
  # Create file names.
  nbObj <- length(object)
  if(use.object.name){
    # Copy object names to name variable.
    name <- object
    
    # Check if not vector.
  } else if (length(name) == 1){
    
    # Split using pipe separator, if present.
    name <-  unlist(strsplit(name, "\\|"))
    
    if(length(name) == nbObj){
      # If equal length, make names.
      name <- make.names(name, unique=TRUE)
    } else {
      # If not equal length, use first and make names.
      name <- make.names(rep(name[1], nbObj), unique=TRUE)
    }
  }

  if(debug){
    print(paste("Number of objects:", nbObj))
    print("name:")
    print(name)
  }
  
  # Add trailing path separator if not present.
  if(substr(path, nchar(path), nchar(path)+1) != .separator){
    path <- paste(path, .separator, sep="")
  }
  
  # EXPORT -------------------------------------------------------------------

  if(debug){
    print("EXPORT")
  }
  
  for(o in seq(along=object)){
    
    # Check file format.
    if(ext == ".RData"){

      # Construct complete file name.
      fullFileName <- paste(path, name[o], ".RData", sep="")
      
      # Overwrite.
      if(overwrite){
        ok <- TRUE
      } else {
        # Check if file exist.
        ok <- .uniqueFile(object=object[o], name=name[o], path=fullFileName)
      }
      
      # Check if ok to export.
      if(ok){
        # Save as R-object.
        save(list=object[o], file=fullFileName, envir=env)
      } else {
        # Write warning.
        warning(text=paste("Object '", object[o], "' was not exported.\n",
                           "The file '", fullFileName,"' already exist!",sep=""))
      }         
      
    } else if (ext == "auto") { # Automatic detection based on class type.

      if(any(class(get(object[o], envir=env)) == "data.frame")){
        
          # Construct complete file name.
          fullFileName <- paste(path, name[o], ".txt", sep="")
          
          # Overwrite.
          if(overwrite){
            ok <- TRUE
          } else {
            # Check if file exist.
            ok <- .uniqueFile(object=object[o], name=name[o], path=fullFileName)
          }
          
          # Check if ok to export.
          if(ok){
            # Save as text file.
            write.table(x=get(object[o], envir=env),
                        file = fullFileName,
                        append = FALSE, quote = FALSE, sep = delim,
                        dec = ".", row.names = FALSE,
                        col.names = TRUE)
          } else {
            # Write warning.
            warning(text=paste("Object '", object[o], "' was not exported.\n",
                               "The file '", fullFileName,"' already exist!",sep=""))
          }         
          
      } else if(any(class(get(object[o], envir=env)) == "ggplot")){
        
        # Construct complete file name.
        fullFileName <- paste(path, name[o], ".png", sep="")
        
        # Overwrite.
        if(overwrite){
          ok <- TRUE
        } else {
          # Check if file exist.
          ok <- .uniqueFile(object=object[o], name=name[o], path=fullFileName)
        }
        
        # Check if ok to export.
        if(ok){
          # Save as image.
          png(filename=fullFileName, width=width, height=height, res=res)
          plot(get(object[o], envir=env))
          dev.off()
        } else {
          # Write warning.
          warning(text=paste("Object '", object[o], "' was not exported.\n",
                             "The file '", fullFileName,"' already exist!",sep=""))
        }         
        
      } else {
        
        warning(paste("Object",object[o],
                      "could not be exported. Object class '",
                      class(get(object[o], envir=env)), "' not supported!",
                      sep=""),
                call. = TRUE, immediate. = FALSE, domain = NULL)
        
      }
    } else {

      stop(paste("Object",object[o],
                    "could not be exported. File extension '",
                    ext, "' not supported!",
                    sep=""),
              call. = TRUE)
      
    }
      
  }
 
  if(debug){
    print("failed")
    print(failed)
  }
  # If 0 rows all objects were sucessfully exported.
  if(nrow(failed) == 0){
    failed <- NA
  }
  
  return(failed)
  
}
