#' @rdname replacemetvar
#' @export
replacemetdata <- function (metdfr, 
                            oldmetfile = "met.dat", 
                            columns=NULL,
                            newmetfile = "met.dat", 
                            khrs=NA,
                            setdates=TRUE){
    
    metlines <- readLines(oldmetfile)
    datastart <- grep("DATA START", metlines, ignore.case=TRUE)
    
    preamble <- readLines(oldmetfile)[1:datastart]
    
    if(is.na(khrs))
      khrs <- readPAR(oldmetfile,"khrsperday","metformat")
    
    N <- nrow(metdfr)
    if(N %% khrs != 0){
      extralines <- N %% khrs
      metdfr <- metdfr[1:(N-extralines),]
    }
  
    writeLines(preamble, newmetfile)
    
    # set end date based on how many rows of data available.
    if(setdates){
      startdate <- readPAR(oldmetfile,"startdate","metformat")
      startDate <- as.Date(startdate[1], "'%d/%m/%y'")
      if(is.na(startDate))
        startDate <- as.Date(startdate[1], "%d/%m/%y")
      enddate <- startDate + N/khrs
      replacePAR(newmetfile, "enddate","metformat", format(enddate, "%d/%m/%y"))  
    }
    
    replacePAR(newmetfile, "nocolumns","metformat", ncol(metdfr))	
    replacePAR(newmetfile, "khrsperday","metformat", khrs)  
    
    if(!is.null(columns))
      replacePAR(newmetfile,"columns","metformat",columns,noquotes=TRUE)	
    
    g <- readLines(newmetfile,100)
    g[grep("data start",g,ignore.case=TRUE)] <- "DATA STARTS"
    writeLines(g,newmetfile)
    
    write.table(metdfr, newmetfile, sep = " ", row.names = FALSE, 
        col.names = FALSE, append = TRUE)
}

