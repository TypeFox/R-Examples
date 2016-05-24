label <- 
function(x, value=NULL, data=mydata) {

  x.name <- deparse(substitute(x)) 
  options(xname = x.name)

  # get data frame name
  dname <- deparse(substitute(data))
  options(dname = dname)

  if (nchar(x.name)>0) if (!exists(x.name, where=data)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
    "For data frame: ", dname, "\n",
    "This variable does not exist: ", x.name, "\n")
  }

  # get conditions and check for data existing
  xs <- .xstatus(x.name, dname)
  in.global <- xs$ig 

  # see if the data frame exists, if x not in Global Env or function call
  if (!in.global) {
    if (!exists(dname)) {
      if (dname == "mydata") 
        txtA <- ", the default data frame name, " else txtA <- " "
      txtB1 <- "Create the labels by reading with Read and labels option\n"
      txtB <- paste(txtB1, sep="")
      cat("\n"); stop(call.=FALSE, "\n","------\n",
          "Data frame ", dname, txtA, "does not exist\n\n", txtB, "\n")
    }
  }

  if (is.null(value)) {  # display an existing label
    if (!is.null(x.name)) {
      gl <- .getlabels()
      lbl <- gl$xl
      if (is.null(lbl)) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "The variable label does not exist for variable: ", x.name, "\n\n")
      }
      cat(x.name, ": ", lbl, "\n", sep="")
     }
    else {
      mylabels <- attr(data, which="variable.labels")
      for (i in 1:length(mylabels))
        cat(names(mylabels)[i], ": ", mylabels[i], "\n", sep="")
    }
  }

  else {  # assign a label to a var in a data frame and return data frame
    mylabels <- attr(data, which="variable.labels")
    lbl.len <- length(mylabels)
    if (x.name %in% names(mylabels)) { #cat("IS IN\n")
      lbl.index <- which(names(mylabels) == x.name)
      indx <- lbl.index
    }
    else
      indx <- length(mylabels) + 1
    mylabels[indx] <- value
    names(mylabels)[indx] <- x.name
    cat("\n")
    cat("Variable Name:",  names(mylabels)[indx], "\n")
    cat("Variable Label:", mylabels[indx], "\n")
    cat("\n")
    attr(data, which="variable.labels") <- mylabels
    return(data) 
  }
  
  cat("\n")
}
