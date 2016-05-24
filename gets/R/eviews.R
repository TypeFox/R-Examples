eviews <-
function(object, file=NULL, print=TRUE,
  return=FALSE)
{
  out <- list()
  out$object.name <- deparse(substitute(object))

  ##index, data, names:
  out$index <- object$aux$y.index
  out$data <- cbind(object$aux$y, object$aux$mX)
  out$data <- as.data.frame(out$data)
  out$data <- cbind(as.character(out$index), out$data)
  out$names <- c("index", object$aux$y.name, object$aux$mXnames)
  where.mconst <- which(out$names=="mconst")
  if(length(where.mconst) > 0){ out$names[where.mconst] <- "c" }
  colnames(out$data) <- out$names

  ##equation command:
  tmp <- paste(out$names[-1], collapse=" ")
  vcov.type <- NULL
  if(object$aux$vcov.type=="white"){
    vcov.type <- "(cov=white)"
  }
  if(object$aux$vcov.type=="newey-west"){
    vcov.type <- "(cov=hac)"
  }
  out$equation <- paste("equation ", out$object.name,
    ".ls", vcov.type, " ", tmp, sep="")

  ##if print=TRUE and is.null(file):
  if(print && is.null(file)){

    ##EViews code to estimate the model:
    cat("EViews code to estimate the model:\n")
    cat("\n")
    cat(" ", out$equation, "\n")
    cat("\n")

    ##R code to export the data:
    cat("R code (example) to export the data of the model:\n")
    cat("\n")
    cat(paste("  eviews(", out$object.name, ", file='C:/Users/myname/Documents/getsdata.csv')\n", sep=""))
    cat("\n")

  } #close if(print)

  ##if save data:
  if(!is.null(file)){
    write.csv(out$data, file, row.names=FALSE)
    ##if print=TRUE:
    if(print){
      cat("Data saved in:\n")
      cat("\n")
      cat("  ", file, "\n", sep="")
      cat("\n")
      cat("EViews code to estimate the model:\n")
      cat("\n")
      cat(" ", out$equation, "\n")
      cat("\n")
    }
  } #end if(!is.null(file))

  ##out:
  if(return){ return(out) }

}
