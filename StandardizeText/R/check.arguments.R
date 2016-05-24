check.arguments <-
function(original.input,col.original,standard.input,col.standard,regex,codes,match,only.names,na.rm,suggest,print.changes,verbose) {
  
  #Check simple arguments
  if(!check.value(suggest,"string")) {
    cat("Error: 'suggest' must be a character string'\n")
    return()
  }else if(!check.value(match,"boolean")) {
    cat("Error: 'match' must be a boolean value'\n")
    return()
  }else if(!check.value(only.names,"boolean")) {
    cat("Error: 'only.names' must be a boolean value'\n")
    return()
  }else if(!check.value(na.rm,"boolean")) {
    cat("Error: 'na.rm' must be a boolean value'\n")
    return()
  }else if(!check.value(print.changes,"boolean")) {
    cat("Error: 'print.changes' must be a boolean value'\n")
    return()
  }else if(!check.value(verbose,"boolean")) {
    cat("Error: 'verbose' must be a boolean value'\n")
    return()
  }
  
  #Check if original input valid
  col.location <- check.input(original.input,col.original)
  if (is.null(col.location)) {
    return()
  }
  
  #Check if standard input valid
  col.stdnames <- check.input(standard.input,col.standard,"'standard'") 
  if (is.null(col.stdnames)) {
    return()
  }
  
  #Check if original names match standard if requested
  if(match) {
    if(col.stdnames != 0) {
      cat("Error: matched 'standard' must be a vector\n")
      return()
    }
    else if(length(standard.input) != length(unique(standard.input))) {
      cat("Error: matched 'standard' vector contains nonunique elements\n")
      return()
    }
  }
  
  #Check if regex valid
  if(!check.value(regex,"vector")) {
    cat("Error: 'regex' must be a vector\n")
    return()
  } else if (length(regex) != length(unique(regex))) {
    cat("Error: 'regex' contains nonunique entries\n")
    return()
  }
             
  #Check if codes valid
  if(!check.value(codes,"vector")) {
    cat("Error: 'codes' must be a vector\n")
    return()
  } else if (length(codes) != length(unique(codes))) {
    cat("Error: 'codes' contains nonunique entries\n")
    return()
  }
  
  #If regex and codes provided
  if(!is.null(regex) && !is.null(codes)) {
    #Check if codes match regex
    if (length(regex) != length(codes)) {
      cat("Error: 'regex' and 'codes' different lengths\n")
      return()
    }
    #Check if standard names match regex if requested
    if (match && (length(regex) != length(standard.input))) {
      cat("Error: 'standard' and 'regex' different lengths\n")
      return()
    }
  }
  return(list(col.location,col.stdnames))
}
