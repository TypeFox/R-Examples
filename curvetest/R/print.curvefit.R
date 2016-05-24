print.curvefit <-
structure(function(x,...){
      fits<-x
      cat("Class of 'curvefit', curve fitting results.....\n")
      cat(paste("formula="));   
      print(fits$formula);
      x<-as.character(deparse(substitute(fits)))
      cat(paste(nrow(fits$data.model), "observations used...\n"))
      cat(paste("Other information can be accessed via functions str() or names().\n",sep=""))
 }, modifiers = "public")
