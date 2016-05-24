createJobCSV <- function(existing.csv, csv.location, overwrite=FALSE, append=FALSE, sep=",",row.names=FALSE, col.names=FALSE) {
  if (is.null(existing.csv)) {
    ex.df = data.frame()
    ex.v = c(1:6)
    ex.v[1] <- "DATE"
    ex.v[2] <- "APIKEY"
    ex.v[3] <- "FILENAME"
    ex.v[4] <- "LANGUAGE"
    ex.v[5] <- "JOBID"
    ex.v[6] <- "TRANSCRIPT"
    ex.df <- rbind(ex.df,ex.v)
    names(ex.df) <- ex.v
    row.names(ex.df) <- NULL

    write.table(ex.df,
                file = csv.location, 
                append=append, sep=sep, 
                row.names=FALSE, 
                col.names = FALSE)
    return(TRUE)
  }
  else {
    print("An existing file is supplied, so a new one will not be created")
    return(FALSE)
  }
}