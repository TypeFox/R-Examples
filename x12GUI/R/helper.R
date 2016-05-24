#tests if a given value is a valid periodvalue for x12 i.e. if it is a integer number
#or the abbrevation of a month name 
isPeriod <- function(s){
  (grepl("^(\\s)*[0-9]+(\\s)*$",s) || 
    tolower(s) %in% c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'))
}

#transforms a given periodvalue in a month name i.e. if it is already a itneger value it is passed thru
#but in case of a monthly abbreveation it is changed to a number
as.Period <- function(s){
  if(isNumber(s))return(as.numeric(s))
  else{
    s = trim(s)
    if(tolower(s) %in% c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')){
      return(which(tolower(s) == c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')))
    }
  }
  NULL
}

#checks per regular expression if string contains a integer
#used to make code more readable
isNumber <- function(s, integer=TRUE){
  if(integer==TRUE){
    grepl("^(\\s)*[0-9]+(\\s)*$",s)
  }
  else {
    grepl("^(\\s)*[0-9]+(\\.[0-9]+)?(\\s)*$",s)
  }
}

#checks if list contains NULL-element
containsNULL <- function(s){
  nulls <- sapply(s, FUN=function(s){is.null(s)})
  if(length(unique(nulls))==2||nulls[1]==TRUE)return(TRUE)
  else return(FALSE)
}

#calculates the amount of seasonalperiods of a specific timeseries in 
#interval c(startyear, startperiod, endyear, endperiod)
calcPeriods <- function(interval, tss){
  f <- frequency(tss)
  points <- (interval[3]-(interval[1]+1))*f + (f-interval[2]) + interval[4]
  return(points)
}

#splits the regvariables from getP() into two seperate lists, one containing the regvariables for the regvariables entry
#the other the manual outliners for the corresponding gui
splitRegvariables <- function(rv){
  outliers <- vector()
  regvariables <- vector()
  
  for(s in rv){
    if(grepl("^(\\s)*(AO|LS|TC)[0-9]*\\.([0-9]*|[A-Za-z]{3})(\\s)*$",s)==TRUE){
      outliers <- append(outliers, trim(s))
    }
    else{
      regvariables <- append(regvariables, trim(s))
    }
  }	
  return(list(outliers=outliers, regvariables=regvariables))
}

#splits a list like list("AO1950.5","LS94234.34") into list(c("AO",1950,5),c("LS",94234,34))
splitOulierstring <- function(os){
  lapply(os, FUN=function(s){
    k <- strsplit(s,"\\.")
    return(c(substr(k[[1]][1],1,2),gsub("(AO|LS|TC)","",k[[1]][1]),k[[1]][2]))
  })
}

#checks per regular expression if string is empty(only consists of whitespaces)
#used to make code more readable
isEmpty <- function(s){
  grepl("^(\\s)*$",s)
}

#splits a string of form (param1 param2 param3) into a vector(c(param1,param2,param3))
cutParam <- function(c){
  out <- strsplit(c,"((\\s)+)|((\\s)*[,](\\s)*)")
  ret <- as.vector(sapply(out, FUN = function(s)str_trim(s)))
}

#concatents a vector of params into csv form
concParam <- function(c){
  paste(c, collapse=" ")
}

#removes leading/tailing whitespaces
trim <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

#helper methode to simplyfie/more readable set widgets inactive depending on checkbox
toggle <- function(a, b, invert=FALSE){
  if(invert == FALSE)lapply(a, function(s){s$SetSensitive(b$GetActive())})
  else lapply(a, function(s){s$SetSensitive(!b$GetActive())})
}

#exchanges middleparts of a filepath against dots to shorten the text. 
#/folder1/subfolder1/morefolders/muchmorefolders/filename
#/folder1/.../muchmorefolders/filename
capPath <- function(s, length=18){
  if(str_length(s)<=length)return(s)
  paste(str_sub(s,0,floor(length/2)),"...",str_sub(s,str_length(s)-floor(length/2)))
}