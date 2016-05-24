Write <- 
function(ref=NULL, data=mydata, format=c("csv", "R"), ...) {

  format <- match.arg(format)

  dname <- deparse(substitute(data))

  if (!exists(dname, where=.GlobalEnv)) {
    cat("\n");
    if (grepl('"', dname))
      cat(">>> Do NOT have quotes around the data frame name.\n\n")
    stop(call.=FALSE, "\n","------\n",
         "Data frame ", dname, " does not exist.\n")
    cat("\n")
  }

  if (format == "csv") {
    if (is.null(ref))
      file.data <- paste(dname, ".csv", sep="")
    else {
       txt <- ifelse (grepl(".csv", ref), "", ".csv")
       file.data <- paste(ref, txt, sep="")
    }
    write.csv(data, file=file.data, ...)
    .showfile(file.data, c(dname, "data values"))

    mylabels <- attr(data, which="variable.labels") # save variable labels
    if (!is.null(mylabels)) {
      mylabels <- data.frame(mylabels)
      file.lbl <- substr(file.data,1,nchar(file.data)-4)
      file.lbl <- paste(paste(file.lbl,"_lbl",sep=""), ".csv" ,sep="")
      write.table(mylabels, file=file.lbl, col.names=FALSE, dec=".", sep=",")
      .showfile(file.lbl, c(dname, "variable labels"))
    }
  }
  

  else if (format == "R") {
    if (is.null(ref))
      file.data <- paste(dname, ".rda", sep="")
    else {
      txt <- ifelse (grepl(".rda", ref), "", ".rda")
      file.data <- paste(ref, txt, sep="")
    }
    save(list=dname, file=file.data, ...)
    .showfile(file.data, c(dname, "data frame contents"))
  }
  

}
