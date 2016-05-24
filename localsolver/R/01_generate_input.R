#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------


lsp.read.command <- function(param.name, value.class){
  if(value.class == "integer"){
    return("readInt")
  }
  if (value.class == "numeric") {
    return("readDouble")
  }
  stop(sprintf("%s: unsupported parameter value class(%s)", param.name, value.class))
}

#returnFunctionContentClass <- 

generate.input <- function(data, inp.append, dat.append, indexFromZero) {
  range_str <- function(length) {
    return(sprintf("%d..%d", ifelse(indexFromZero, 0, 1), ifelse(indexFromZero, length - 1, length)))
  }  
  
  inp.append("function input() {\n")
  
  for (param.name in names(data)){
    value <- data[[param.name]]
    
    if(class(value) == "list"){
      listLength <- length(value)
      if(listLength <= 0){
        err <- sprintf("The length of the list %s must be bigger than 0!", param.name)
        stop(err)
      }
      
      # we find the first index of the list which is not null to check the types of the list contents
      
      i <- 1
      nonEmptyIndexFound <- TRUE
      while(i <= listLength & nonEmptyIndexFound){
        if(!is.null(value[[i]])){
          nonEmptyIndexFound <- FALSE
        }else{
          i <- i + 1
        }
      }
      
    
      for(i in 1:listLength){
        record <- value[[i]]
        recordLength <- length(record)
        if(recordLength > 0){
          dat.append(paste(record, collapse = "\n"))
          dat.append('\n')
          lsp.read.cmd <- lsp.read.command(param.name, class(record))          
          inp.append(sprintf("\t%s[%d][%s] = %s();\n", param.name, i - ifelse(indexFromZero, 1, 0), range_str(recordLength), lsp.read.cmd))
        }
      }
    }else{
  
      if(is.null(dim(value))) { # 1-dim object or value
        lsp.read.cmd <- lsp.read.command(param.name, class(value))
  
        dat.append(paste(value, collapse = "\n"))
        dat.append('\n')
        
        len <- length(value)
        if (len == 1) {        
          inp.append(sprintf("\t%s = %s();\n", param.name, lsp.read.cmd))
        } else { # 1-dim vector len >= 2        
          inp.append(sprintf("\t%s[i in %s] = %s();\n", param.name, range_str(len), lsp.read.cmd))
        }      
        next
      }
      
      # array of 2 or 3 dimensions
      dimLength <- length(dim(value))
      
      # array of dim = 2
      
      if(dimLength == 2){
  
        # writing data into text file
        for(i in 1:nrow(value)){   # 1:dim(value)[1]
          dat.append(paste(value[i,], collapse = "\n"))
          dat.append('\n')
        }
        
        # writing commands into lsp file
        lsp.read.cmd <- lsp.read.command(param.name, class(value[1,1]))
        
        inp.append(sprintf("\tfor[i in %s]{\n", range_str(nrow(value))))
        inp.append(sprintf("\t\t%s[i][%s] = %s();\n", param.name, range_str(ncol(value)), lsp.read.cmd))      
        inp.append("\t}\n\n")
      }else if(dimLength == 3){
        
        dimensions <- dim(value)
        for (i in 1:dimensions[1]){
          for (j in 1:dimensions[2]){
            dat.append(paste(value[i,j], collapse = "\n"))
            dat.append('\n')
  
          }
        }
        
        lsp.read.cmd <- lsp.read.command(param.name, class(value[1,1,1]))
        inp.append(sprintf("\tfor [i in %s]{\n", range_str(dimensions[1])))
        inp.append(sprintf("\t\tfor [j in %s]{\n", range_str(dimensions[2])))      
        inp.append(sprintf("\t\t\t%s[i][j][%s] = %s();\n", param.name, range_str(dimensions[3]), lsp.read.cmd))
        inp.append("\t}\n\n")
      }
    } # end of iteration with paramName

  }
  inp.append("}")
}

