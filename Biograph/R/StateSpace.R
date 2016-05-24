StateSpace <-
function (d,newnamstates) 
{   # d can be an entire Biograph object (e.g. Bdata) of a vector of state sequences (e.g. Bdata$path)
	if (missing(newnamstates)) newnamstates <- NULL
	 # Get the state space from the data (by inspecting Bdata$path)
      if (class(d)=="data.frame") d <- unique(as.character(d$path)) else {if (class(d)=="character") d = d else stop ("Error in StateSpace" ) }
    namstates <- array(NA, 50)
    numstates <- 0
    for (i in 1:length(d)) {
        str_char <- stringf(d[i])
        for (j in 1:length(str_char)) {
            if (str_char[j] %in% namstates) 
                test22 <- 1
            else {
                numstates <- numstates + 1
                namstates[numstates] <- str_char[j]
            }
        }
    }
    namstates <- subset(namstates, !is.na(namstates))
    absorb <- namstates
    for (j in 1:numstates) {
    	# State j is in the path and is not the last element
    	# State j i aborbing if it always is last element of path
        for (i in 1:length(d)) {
            if (namstates[j] %in% stringf(d[i]) & namstates[j] != 
                substr(d[i], nchar(d[i]), nchar(d[i]))) {
                absorb[j] <- NA
                break
            }
        }
    }
    absorbstates <- subset(absorb, !is.na(absorb))
    if (length(absorbstates)==0) absorbstates <- NULL

# Check whether new state sequence is provided
   if (is.null(newnamstates)) newnamstates <- namstates
# Consistency check: Check whether labels in newnamstates are all in namstates 
   z<- match (newnamstates, namstates)

   if (NA %in% z) 
      stop ("No match between namstates and newnamstates. Enter a new global variable 'namstates'. ") else  
      {namstates <- newnamstates }  #assign("namstates",newnamstates,envir=.GlobalEnv)}

   # assign("numstates",numstates,envir=.GlobalEnv)
    return(list(namstates = namstates,
                absorbstates = absorbstates))
}
