#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------

parse.ls.output <- function(output.exprs, out.file){
  lines <- readLines(con=out.file)
  lineIx <- 1
  
  out.values <- list()
  for(exprName in names(output.exprs)) {
    desc <- output.exprs[[exprName]]
    
    while(lineIx <= length(lines) && substring(lines[[lineIx]], 1, 2) != '$:') {
      lineIx <- lineIx + 1
    }
    
    if (lineIx > length(lines)) {
      stop("Premature LocalSolver output end detected.")
    }
    
    exprValueTxt <- substring(lines[[lineIx]], 3)
    lineIx <- lineIx + 1
    
    out.values[[desc$expr]] <- parse.var.out.values(desc$dimensions, exprValueTxt)    
  }
    
  return(out.values)
}


parse.var.out.values <- function(var.dimensions, outputValues){
  
  outputValues <- as.numeric(unlist(strsplit(outputValues, split = " ")))
  
  if(length(var.dimensions) == 1){
    return(outputValues) # A number or a vector    
  }
  
  # array
  if(length(var.dimensions) == 2){
    iter <- 1
    newArray <- array(dim = var.dimensions)
    for(i in 1:var.dimensions[1]){
      for(j in 1:var.dimensions[2]){
        newArray[i, j] <- outputValues[iter]
        iter <- iter + 1
      }
    }
    
    return(newArray)      
  }
  
  # 3-dimensional array 
  iter <- 1
  newArray <- array(dim = var.dimensions)
  for(i in 1:var.dimensions[1]){
    for(j in 1:var.dimensions[2]){
      for(k in 1:var.dimensions[3]){
        newArray[i, j, k] <- outputValues[iter]
        iter <- iter + 1            
      }
    }
  }
  
  return(newArray)      
}

