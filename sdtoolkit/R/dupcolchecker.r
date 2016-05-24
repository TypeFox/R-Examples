
#Code for checking if input columns are duplicates of one another
#x is the input matrix
#exact determines whether to use identical() vs all.equal, which allows slightly more leeway
#depending on user input, returns either a new matrix with duplicates removed, 
#or a vector indices that are duplicated
#or the same matrix, with no

dupcolchecker <- function(x, exact=FALSE){
  
  d <- ncol(x)
  
  dupmat <- matrix(data=FALSE, ncol=d, nrow=d-1)
  
  if (exact){
    
    for (i in 1:(d-1)){
    
      for (j in (i+1):d){
      
        dupmat[i,j] <- identical(x[,i], x[,j])
        
      }
    
    }

  } else{  
    
    for (i in 1:(d-1)){
    
      for (j in (i+1):d){
      
        dupmat[i,j] <- isTRUE(all.equal(x[,i],x[,j]))
        
      }
    
    }

  }
 
  #identify rows which have dups:
  hasdups <- apply(dupmat,1,any, na.rm=TRUE)
 
  if(any(hasdups)){   #only bother with stuff if there actually are duplicates
 
    #identify if there are more than pairwise duplications and remove the latter:
    numdups <- apply(dupmat,1,sum, na.rm=TRUE)
    
    dupdups <- c()  
      
    for (i in (1:(d-1))[hasdups]){
    
      if(numdups[i]>1){ #skip display of the extra columns
      
        dupdups <- c(dupdups,(1:d)[dupmat[i,]])
        
      }
      
    }
        
    dupdups <- unique(dupdups)
    
    hasdups[dupdups] <- FALSE
   
   
    #print warnings 
   
    cat("\n")
    cat("==== WARNING ====") 
    cat("The following sets of columns appear to be duplicates of each other. \n")
    nicecat("(Note that the column numbers may be different than those in the original dataset
    if there were also also columns with identical variables that were removed.)")
    cat("\n")
   
    for (i in (1:(d-1))[hasdups]){
    
      dupcols <- c(i,(1:d)[dupmat[i,]])
    
      cat("Column numbers", dupcols," with column names:", colnames(x)[dupcols], "\n")
                                                      
    }
    
    cat("\n")
    cat("Duplicate columns are likely to either break PRIM and/or give poorly","\n","interpretable boxes.","\n")
    cat("For each set of duplicates, please choose the column you would like to keep, \n")
    cat("and delete the others from your input matrix.","\n", "\n")
    cat("OR, you can have this program automatically preserve the lowest index columns","\n")
    cat("and automatically delete the repeats at higher indices.","\n")
    
    nonorig <- apply(dupmat,2,any,na.rm=TRUE)
    
    cat("Specifically, it would remove the following dimension indices:","\n")
    cat("\n")
    cat((1:d)[nonorig],fill=TRUE)
    cat("\n") 
    
    cat("Otherwise, it will return a vector containing the indices just printed above.")
    cat("\n")
    
    
    yayornay <- readline("Would you like to have columns removed?  Enter \"y\" or \"n\" please. \n")
    if (yayornay=="y"){
      xnew <- x[,!nonorig]
      return(xnew)
    } else { 
      return((1:d)[nonorig])
    }

  }else {
    cat("\n","Checked for duplicate columns - none found.","\n","\n")
    return("nodups")
  }
  

}
  
