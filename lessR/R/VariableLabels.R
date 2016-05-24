VariableLabels <- 
function(x, value=NULL, data=mydata, quiet=getOption("quiet")) {

  x.name <- deparse(substitute(x)) 
  options(xname = x.name)

  # get data frame name
  dname <- deparse(substitute(data))
  options(dname = dname)

  fmt <- "none"  # see if x is a file reference
  if (grepl(".csv", x.name)) fmt <- "csv" 
  if (grepl(".xlsx", x.name)) fmt <- "Excel" 

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
    # see if variable exists in the data frame
    x.is.var <- FALSE
    if (nzchar(x.name)) if (exists(x.name, where=data)) 
      x.is.var <- TRUE


  # ------------------------------------------

  # display existing label or all labels
  if (is.null(value)  &&  !in.global  &&  fmt=="none") {  
    if (nzchar(x.name)) {
      gl <- .getlabels()
      lbl <- gl$xl
      if (is.null(lbl)) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "The variable label does not exist for variable: ", x.name, "\n\n")
      }
      cat("\n")
      cat(x.name, ": ", lbl, "\n", sep="")
     }
    else {  # display all labels
      mylabels <- attr(data, which="variable.labels")
      cat("\n")
      for (i in 1:length(mylabels))
        cat(names(mylabels)[i], ": ", mylabels[i], "\n", sep="")
    }
  }

  else {  # assign labels to vars in a data frame and return data frame

    if (fmt == "none")  # no external file
      if (!x.is.var)  # x is a file var names and labels from the console
        mylabels <- read.csv(text=x, col.names=c("Var", "label"), row.names=1,
                             header=FALSE)
      else {  # x is a single variable
        mylabels <- c(x.name, value)
      }
    else if (fmt == "csv")  # x is a file name
      mylabels <- read.csv(x, col.names=c("Var", "label"), row.names=1,
                           header=FALSE)
    else { 
      if (fmt=="Excel"  &&  grepl("http", x, fixed=TRUE)) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
            "At this time the underlying read_excel function\n",
            "does not support reading Excel files from the web.\n",
            "Can use the  download.file  function to download to the\n",
            "local file system.\n\n",
            "  download.file(\"", x, "\", \"MYFILE.xlsx\")\n\n",
            "Replace MYFILE with the desired file name.\n",
            "Enter  getwd()  to see where the file was saved. \n\n")
      }
      mylabels <- read_excel(x, col_names=c("Var", "label"))
      mylabels <- data.frame(mylabels, row.names=1)
    }

    if (in.global  ||  fmt!="none") { # transfer table of labels and maybe units to data
      
      attr(data, which="variable.labels") <- as.character(mylabels$label)
      names(attr(data, which="variable.labels")) <- as.character(row.names(mylabels))
      if (ncol(mylabels) == 2) {
        attr(data, which="variable.units") <- as.character(mylabels$unit)
        names(attr(data, which="variable.units")) <- as.character(row.names(mylabels))
      }
      if (!quiet) details.brief(data)
    }

    else {  # single assignment of a variable label
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
    }

    return(data) 
  }
  
  cat("\n")
}
