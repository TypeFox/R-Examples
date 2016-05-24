generateName <- function(base="tmp",
                            oldnames, ## typically for creating name of new variable
                            pos=".GlobalEnv" ## typically not used: for creating name of object for debugging
                            ) {
  if ( ! missing(oldnames)) { ## because presumably the test should be TRUE if preexisting was provided yet is NULL
    allmatches <- pmatch(x=base,oldnames)
  } else {
    pattern <- paste(base,"*",sep="")
    allmatches <- ls(pattern=pattern,pos=pos)
  }
  allremainders <- substring(allmatches,nchar(base)+1) 
  allremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  if (length(allremainders) == 0L) {
    validname <- base 
  } else {
    num <- max(allremainders)+1
    validname <-paste ( base , num , sep="") 
  }
  return(validname)
}

generateFileName <- function(base="tmp",ext="") { ## for a file
   pattern <- paste(base,"*",ext,sep="")
   allmatches <- dir(pattern=pattern)
   allremainders <- substring(allmatches,nchar(base)+1)
   allremainders <- unlist(strsplit(allremainders,ext)) ## removes the extension from the remainder 
   allremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
   if (length(allremainders) == 0) {
            num <- 0
   } else num <- max(allremainders)+1
   validFileName <-paste ( base , num , ext,sep="") 
   return(validFileName)
}
