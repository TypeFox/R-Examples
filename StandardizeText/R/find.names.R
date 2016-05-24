find.names <-
function(original.names,regex.codes,standard.names=NULL,set.na=FALSE,suggest=FALSE,verbose=FALSE) {
  
  #Format the input names
  formatted.names <- preformat.names(original.names)
  #Grep function
  #With anchor tags ("^" and "$"); ensures only complete matches
  grep.regex <- function(x) {
    y <- grep(paste("^",x,"$",sep=""),formatted.names)
    if(length(y) > 1) {y <- integer(0)}
    return(y)
  }
  name.locations <- as.numeric(sapply(regex.codes[,1],grep.regex))
  
  #Find location of name matching a particular regex
  name.matches <- rep(0,length(original.names))
  name.matches[sort(na.omit(name.locations))] <- 1
  notfound.locations <- which(name.matches==0)
  notfound.values <- data.frame(code=numeric(0),name=character(0))
  notfound.names <- original.names[notfound.locations]
  notfound.length <- length(notfound.locations)
  
  initial <- is.null(standard.names)
  #If unrecognized names
  if (notfound.length > 0) {
    #Generate codes for unclassified names
    unclassified.codes <- autogen.codes(notfound.length)
    if (initial) {
      unclassified.codes <- paste("X0",unclassified.codes,sep="")
    } else {
      unclassified.codes <- paste("X1",unclassified.codes,sep="")
    }
    string.xchar <- function(y,n) {sum(grepl(paste("^(X{",as.character(n),"}[0-9])",sep=""),y))>0}
    n = 1
    while (string.xchar(regex.codes[,2],n)) {
      unclassified.codes <- paste("X",unclassified.codes,sep="")
      n = n + 1
    }
    notfound.codes <- unclassified.codes[seq(notfound.length)]
    notfound.values <- data.frame(code=notfound.codes,name=notfound.names)    
    name.locations <- c(name.locations,notfound.locations)

    #Suggestions
    if(suggest) {
	  exclude.codes <- NULL
      if(!initial) {
        exclude.codes <- standard.names[which(is.na(standard.names[,2])),1]
      }
      notfound.formattednames <- formatted.names[notfound.locations]
      grep.suggestions <- function(x) {
        y <- grep(x,notfound.formattednames)
        if(length(y) > 1) {y <- integer(0)}
        return(y)
      }
      unused.regex <- regex.codes[which(is.na(name.locations)),]
	  if(length(exclude.codes)>0){
        exclude.unused <- which(unused.regex[,2] %in% exclude.codes)
        unused.regex <- unused.regex[-exclude.unused,]
      }
      suggestion.locations <- as.numeric(sapply(unused.regex[,1],grep.suggestions))
      suggestion.originals <- notfound.names[suggestion.locations]
      suggestion.codes <- as.character(unused.regex[,2])
      suggestion.norepeatednames <- !(duplicated(suggestion.originals)|duplicated(suggestion.originals,fromLast=TRUE))
      notfound.suggestions <- na.omit(data.frame(name=suggestion.originals[suggestion.norepeatednames],suggestion=suggestion.codes[suggestion.norepeatednames]))
      notfound.values <- merge(notfound.values,notfound.suggestions,all=T)
      notfound.values <- notfound.values[,c(2,1,3)]
    } else {
      notfound.values <- cbind(notfound.values,suggestion=rep(NA,dim(notfound.values)[1]))
    }
    nosuggestions.locations <- which(is.na(notfound.values[,3]))
    
    #Print actions taken
    if(length(nosuggestions.locations) > 0) {
      if (!initial) {
        if (set.na) {
          action.taken <- "removed"
        } else {
          action.taken <- "left unchanged"
        }
        if(verbose || set.na) {
          cat(sprintf("\nThe following names were not recoginized and %s:\n",action.taken))
          print(as.character(notfound.values[nosuggestions.locations,2]))
        } else {
          cat(sprintf("\nNote: %d names were not recoginized and %s.\n",length(notfound.values[nosuggestions.locations,2]),action.taken))
        }
      }
      if (set.na) {
        notfound.values[nosuggestions.locations,2] <- rep(NA,length(notfound.values[nosuggestions.locations,2]))
      }
    }
  }
  return(list(name.locations,notfound.values))
}
