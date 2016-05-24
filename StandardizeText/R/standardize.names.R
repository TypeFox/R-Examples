standardize.names <-
function(original.input,col.original=NULL,standard.input,col.standard=NULL,regex=NULL,codes=NULL,match=FALSE,only.names=FALSE,na.rm=FALSE,suggest="prompt",print.changes=TRUE,verbose=FALSE) {
  
  #Check arguments
  arguments.valid <- check.arguments(original.input,col.original,standard.input,col.standard,regex,codes,match,only.names,na.rm,suggest,print.changes,verbose)
  if(is.null(arguments.valid)){
    return()
  }
  col.location <- arguments.valid[[1]]
  col.stdnames <- arguments.valid[[2]]
  
  #Original input
  #If a vector 
  if(col.location == 0) {
    original.names <- levels(factor(as.character(original.input)))
  #If a dataframe
  } else if(col.location > 0) {
    original.names <- levels(factor(as.character(original.input[,col.location])))
  } else {
    return()
  }
 
  #Standard input
  #If a vector 
  if(col.stdnames == 0) {
    if (match) {
      standard.names <- standard.input
    } else {
      standard.names <- levels(factor(standard.input))
    }
  #If a dataframe
  } else if(col.stdnames > 0) {
    standard.names <- levels(factor(standard.input[,col.stdnames]))
  } else {
    return()
  }

  #Regex
  if(!is.vector(regex) && is.null(regex)) {
    regex <- preformat.names(standard.names)
  } 
  #Codes
  if(!is.vector(codes) && is.null(codes)) {
    codes <- autogen.codes(length(regex))
  }
  #Regex/codes
  regex.codes <- data.frame(regex,codes)
  names(regex.codes) <- c("regex","code")
  regex.codes$regex <- as.character(regex.codes$regex)
  regex.codes$code <- as.character(regex.codes$code)
  
  #Generate standard names
  if (match) {
    standard.values <- data.frame(codes,standard.names)
    names(standard.values) <- c("code","name")
  } else {
    standard.values <- generate.names(standard.names,regex.codes,set.na=na.rm,suggest=FALSE,verbose=verbose)
  }
  
  #Generate standardized names
  use.suggestions <- (suggest!="none" && suggest!="n")
  standardized.values <- generate.names(original.names,regex.codes,standard.values,set.na=na.rm,suggest=use.suggestions,verbose=verbose)
  suggested.values <- data.frame(original=character(0),suggestion=character(0))
  if(is.list(standardized.values)) {
    suggested.values <- standardized.values[[2]]
    suggested.values[,1] <- as.character(suggested.values[,1])
    suggested.values[,2] <- as.character(suggested.values[,2])
    standardized.values <- standardized.values[[1]]
  }
  
  #Obtain original names
  original.values <- original.names
  modified.values <- standardized.values
  
  #Determine changed names
  for(i in 1:length(original.values)) {
    if(is.na(original.values[i]) || is.na(modified.values[i]) || original.values[i] == modified.values[i]) {
      original.values[i] <- NA
      modified.values[i] <- NA
    }
  }
  original.values.nona <- na.omit(original.values)
  modified.values.nona <- na.omit(modified.values)
  if(length(original.values.nona) > 0) {
    if(print.changes){
      cat("\nThe following names were changed:\n")
      print(data.frame(Original=na.omit(original.values),Modified=na.omit(modified.values)))
    } else {
      cat(sprintf("\n%d names were changed.\n",length(original.values.nona)))
    }
  } else{
    cat("\nNo names were changed.\n")
  }

  #Suggestions
  if(use.suggestions && dim(suggested.values)[1] > 0) {
    standardized.values <- apply.suggestions(suggest,suggested.values,original.names,standardized.values,na.rm,print.changes)
  }
  
  #If only returning vector of names
  if(only.names) {
    if(na.rm) {
      return(as.vector(na.omit(standardized.values)))
    }else{
      return(standardized.values)
    }
  }
  
  #Copy of input
  output.values <- original.input  
  #If input a vector
  if (col.location == 0) {
    #Factor
    output.values <- factor(output.values)
    #Change name levels
    levels(output.values) <- standardized.values
    #If removing NA names
    if (na.rm) {
      output.values <- output.values[!is.na(output.values)]
      output.values <- factor(output.values)
    }
    #Convert to vector
    output.values <- as.vector(output.values)
  #If input a dataframe
  } else {
    #Factor
    output.values[,col.location] <- factor(output.values[,col.location])
    #Change name levels
    levels(output.values[,col.location]) <- standardized.values
    #Fix structure
    output.values<-data.frame(output.values)
    #Fix modified column name
    names(output.values)[col.location] <- names(original.input)[col.location]
    #If removing NA names
    if (na.rm) {
      output.values <- output.values[!is.na(output.values[,col.location]),]
      output.values[,col.location] <- factor(output.values[,col.location])
    }
  }

  #Return vector or data frame
  return(output.values)
}
