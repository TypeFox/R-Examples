# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     



# "Passion makes the world go round. Love just makes it a safer place."
# - Ice T, The Ice Opinion
	

setMethod(f="show",signature("DependencyModel"),
  function(object){
    cat("***", object@method, "dependency model for window size: ")
    windowSize <- getWindowSize(object)
    if (length(windowSize) == 1) {
      cat(windowSize)
    }
    else {
      cat(windowSize[1],", ",windowSize[2], sep="")
    }
    cat(" *** \n")
  
    #Score
    cat("Dependency score:", object@score,"\n")

    if (is.null(object@W$X)){
      #W print
      cat("- ")
      cat(matrix.print(object@W$total,"W"),"\n", sep="")

      #Phi
      cat("- ")
      cat(matrix.print(object@phi$total,"phi"),"\n", sep="")
    }
    else {
      #WX print
      cat("- ")
      cat(matrix.print(object@W$X,"WX"),"\n", sep="")

      #WY print
      cat("- ")
      cat(matrix.print(object@W$Y,"WY"),"\n", sep="")
	
      #Phi X
      cat("- ")
      cat(matrix.print(object@phi$X,"phiX"),"\n", sep="")
    
      #Phi Y
      cat("- ")
      cat(matrix.print(object@phi$Y,"phiY"),"\n", sep="")
	}
    cat("************************************************\n")
  }
)


matrix.print <- function(matrix,name){
	# Prints size and 4 first values of matrix for show-methods

	if(any(is.na(matrix))){
		return(cat(name,": NA",rep=""))
	}
	else { 
		string1 <- paste(name, ": [1:", dim(matrix)[1], ", 1:",  dim(matrix)[2] ,"] ", sep = "")
		values <- format(matrix[1:min(4,length(matrix))],digits=3)
		string2 <- ""
		if(length(matrix) > 4) string2 <- "..."
		return(cat(string1,values,string2))
	}
}

