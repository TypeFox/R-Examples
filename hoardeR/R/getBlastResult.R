getBlastResult <- function(RID){
    timeElapsed <- 0
    tries <- 0
    newError <- TRUE
    while(newError & tries < 6){
      newError <- FALSE
      tryCatch(blastRes <- scan(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?RID=",RID, "&SHOW_OVERVIEW=no&FORMAT_TYPE=XML&ALIGNMENTS=0&NCBI_GI=yes&CMD=Get", sep = ""), what = "", sep = "\n", quiet = TRUE)
                  
                , error = function(e){
                          cat("An error occured, try",tries,"\n")
                          Sys.sleep(10)
                          newError <<- TRUE
                          tries <<- tries + 1
                          }
              )
    }
  ready <- FALSE
  if(blastRes[length(blastRes)]=="</html>"){
    # Case that blast run is not ready yet:
    timeElapsed <- strsplit(strsplit(blastRes[grepl("Time since submission",blastRes)],"<td>")[[1]][3],"</td>")[[1]][1]
    cat("Run",RID,":",timeElapsed,"\n")
  } else {
    ready <- TRUE
  }
  res <- list(blastRes=blastRes, ready=ready, time=timeElapsed, RID=RID)
  res
}