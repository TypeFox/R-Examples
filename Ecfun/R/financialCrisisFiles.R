financialCrisisFiles <- function(files=c("22_data.xls", "23_data.xls",
       "Varieties_Part_III.xls", "25_data.xls"), ...){
##
## 0.  checks
##
  fi <- file.info(files)
  oops <- which(is.na(fi[, 1]))
  if(length(oops)>0)
      stop('file not found: ', files[oops[1]])
#  dots
  dots <- list(...)
##
## primary output list
##
  nFiles <- length(files)
  out <- vector(mode='list', length=nFiles)
##
## populate it
##
#  library(gdata)
  for(i in 1:nFiles){
##
## 1.  Read the first sheet
##
#    Con0 <- xls2csv(files[i])
    cat('Read ', files[i], '; ', sep='')
    Contents <- gdata::read.xls(files[i], pattern='Country',
                         stringsAsFactors=FALSE)
##
## 2.  Extract the names of the Countries
##
#    headers <- which(Contents[[2]]=='Country')
#    n. <- nrow(Contents)
#    ctries <- Contents[-(1:headers), 2]
    ctr <- Contents[[2]]
    oops <- (is.na(ctr) | (ctr==''))
    ctries <- ctr[!oops]
##
## 3.  Eliminate blank spaces in names
##
    Ctries <- gsub(' ', '', ctries)
##
## 4.  Find the sheets corresponding to each compressed name
##
    sheets <- gdata::sheetNames(files[i])
    ns <- length(sheets)
    found <- numeric(ns)
    for(js in 1:ns){
        fj <- pmatch(sheets[js], Ctries)
        if(is.na(fj)){
            fjs <- switch(sheets[js],
                         UK=pmatch('UnitedKingdom', Ctries),
                         US=pmatch('UnitedStates', Ctries)
                         )
            if(!is.null(fjs))fj <- fjs
        }
        found[js] <- fj
    }
    nc <- length(Ctries)
    if(nc != sum(!is.na(found))){
        cat('Country name from Contents not found in sheet names\n')
        print(Ctries)
        print(sheets)
        stop()
    }
    Ctri. <- sort(sheets[!is.na(found)])
##
## 5.  Construct the output list
##
    fileList <- vector(mode='list', length=nc)
    for(iCtri in 1:nc){
        fileList[[iCtri]] <- dots
    }
    names(fileList) <- Ctri.
    out[[i]] <- fileList
  }
  names(out) <- files
##
## done
##
  class(out) <- 'financialCrisisFiles'
  out
}

