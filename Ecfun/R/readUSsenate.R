readUSsenate <- function(url.=
"https://en.wikipedia.org/wiki/List_of_current_United_States_Senators",
      stateAbbreviations=Ecdat::USstateAbbreviations,
      fixNonStandard=subNonStandardNames, ...){
##
## 1.  download content
##
#  library(RCurl)
  Start <- paste(date(), ': readUSsenate(', url., ')', sep='')
  cat(Start)
  startTime <- proc.time()
  senate.gov <- try(RCurl::getURL(url.))
  et <- max(proc.time()-startTime, na.rm=TRUE)
  Read <- paste('|', nchar(senate.gov), 'bytes read in',
                round(et, 2), 'seconds\n')
  cat(Read)
  if(class(senate.gov)=='try-error'){
      cat(senate.gov)
      stop()
  }
##
## 2.  readHTMLTable
##
#  library(XML)
  senate <- XML::readHTMLTable(senate.gov, stringsAsFactors=FALSE)
##
## 3.  Find table with 100 rows
##
  len <- sapply(senate, length)
  Senate <- senate[len>0]
  ns <- lapply(Senate, dim)
  NS <- sapply(ns, max)
  ns. <- which.max(NS)
  if(NS[ns.]!=100) {
      warning('Problem with number of senators read at URL = ', url.,
              ';  max dimension = ', NS[ns.], ' != 100')
  }
#  ns. <- which(ns[1, ]>=100)
#  if(length(ns.)!=1){
#      stop('There is not exactly one table with 100 rows in url = ',
#           url.)
#  }
  sen <- Senate[[ns.]]
##
## 4.  Delete constant colulms
##
  dataCols <- sapply(sen, function(x){
      length(unique(x))>1
  } )
  Sen <- sen[dataCols]
##
## 5.  parse State
##
  USstates <- stateAbbreviations$Name
  rownames(stateAbbreviations) <- USstates
  USPS <- grep('USPS', names(stateAbbreviations))
  if(length(USPS)<1)
      stop('USPS codes not in stateAbbreviations()')
  senState <- stateAbbreviations[Sen$State, USPS[1]]
##
## 6.  parse more
##
  Sen$State <- factor(Sen$State)
  Sen$state <- factor(senState)
  Sen$Party <- factor(Sen$Party)
#
  Nm <- camelParse(Sen$Name)
  Nm. <- sapply(Nm, "[", 1)
  Sen$Name <- Nm.
  nms <- names(Sen)
  Nms <- sub("Prior Experience", "Experience", nms)
  Nms2 <- sub("Assumed Office", "assumedOffice", Nms)
  Nms3 <- sub("Born In", "Born", Nms2)
  Nms4 <- sub("End Office", "endOffice", Nms3)
  names(Sen) <- Nms4
##
## 7.  Rearrange to put 'state' after 'State';
##     assume 'state' is the last column
##
  nc <- ncol(Sen)
  o <- c(1, nc, 2:(nc-1)) # assume '
  Sen. <- Sen[o]
##
## 8.  parseName
##
  pN <- parseName(Sen.$Name, fixNonStandard=fixNonStandard, ...)
  SEN <- cbind(Sen., as.data.frame(pN, stringsAsFactors=FALSE))
  SEN
}

