# Utilities ############################################

# Attribute used while printing the model expression
expr <- function(obj){
  return(attr(obj,'expr'))
}

'expr<-' <- function(obj, value){
  attr(obj,'expr') <- value
  return(obj)
}

latex.tab.content <- function(d){
  nrows <- dim(d)[[1]]
  ncols <- dim(d)[[2]]
  
  
  printRow <- function(rd){
    if(is.na(rd[[1]])){
      row <- ''; 
    }else{        
      row <- paste(rd[[1]])
    }
    for(j in 2:ncols){
      if(is.na(rd[[j]])){
        row <- paste(row,'&')
      }else if(is.numeric(rd[[j]])){
        row <- paste(row,sprintf(' & %.2e',rd[[j]]))
      }else{
        
        row <- paste(row,sprintf(' & %s',rd[[j]]))
      }
    }
    cat(paste(row,' \\\\\n'))
  }
  
  printRow(colnames(d))
  for(i in 1:nrows){
    cat(rownames(d)[[i]],' & ')
    printRow(d[i,])
  }
}




group.by.d_start <- function(d){
  d$d_start <- round(d$d_start,2)
  
  evals <- levels(d$eval_name)  
  
  df <- data.frame(d_start <- unique(d$d_start))
  
  nrows <- dim(df)[[1]]
  for(ev in evals){    
    col <- matrix(0,nrows,1)
    colnames(col) <- ev
    for(row in 1:nrows){
      mask <- (abs(d$d_start - df$d_start[[row]])<0.01 & d$eval_name == ev)      
      col[[row]] <- mean(d$d_stop[mask])      
    }
    
    df <- cbind(df,col)
  }
  
  return(df)
}


group.by.samplesPerDim <- function(d){
  
  evals <- levels(d$eval_name)  
  
  df <- data.frame(SamplesPerDim = unique(d$SamplesPerDim))
  
  nrows <- dim(df)[[1]]
  
  for(ev in evals){    
    col <- matrix(0,nrows,1)
    sdCol <- matrix(0,nrows,1)
    colnames(col) <- ev
    colnames(sdCol) <- paste(ev,'_sd')
    for(row in 1:nrows){
      mask <- (d$SamplesPerDim == df$SamplesPerDim[[row]] & d$eval_name == ev)      
      col[[row]] <- mean(d$d_stop[mask])      
      sdCol[[row]] <- sd(d$d_stop[mask])
    }
    
    df <- cbind(df,col,sdCol)    
  }
  
  return(df)
}
