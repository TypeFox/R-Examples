readUSstateAbbreviations <- function(url.=
"https://en.wikipedia.org/wiki/List_of_U.S._state_abbreviations",
      clean=TRUE, Names=c('Name', 'Status', 'ISO', 'ANSI.letters',
          'ANSI.digits', 'USPS', 'USCG', 'Old.GPO', 'AP', 'Other') ){
##
## 1.  download content
##
#  library(RCurl)
  Start <- paste(date(), ': readUSstateAbbreviations(',
                 url., ')', sep='')
  cat(Start)
  startTime <- proc.time()
  abbrev <- try(RCurl::getURL(url.))
  et <- max(proc.time()-startTime, na.rm=TRUE)
  Read <- paste('|', nchar(abbrev), 'bytes read in',
                round(et, 2), 'seconds\n')
  cat(Read)
  if(class(abbrev)=='try-error'){
      cat(abbrev)
      stop()
  }
##
## 2.  Find the primary table
##
#  library(XML)
  Abbrev <- XML::readHTMLTable(abbrev, stringsAsFactors=FALSE)
#
  len <- sapply(Abbrev, length)
  abbr <- Abbrev[len>0]
  ns <- sapply(abbr, dim)
  ns. <- which(ns[1, ]>=50)
  if(length(ns.)!=1){
      stop('There is not exactly one table with more than 50 rows in url = ',
           url)
  }
  Abbr <- abbr[[ns.]]
##
## 3.  The Wikipedia table has problem headers, so fix ...
##
  colNms <- names(Abbr)
  if(all(substring(colNms, 1, 1) == 'V')){
      colNms. <- make.names(Abbr[1,])
      Abbr <- Abbr[-1,]
      names(Abbr) <- colNms.
  }
##
## 4.  Delete obsolete dups
##
  stateReps <- table(Abbr$Name)
  dups <- names(stateReps)[stateReps>1]
  Dups <- which(Abbr$Name %in% dups)
  Dup2 <- regexpr('Obsolete', Abbr[Dups, "Status"])
  Del <- Dups[Dup2>0]
  if(length(Del)>0){
      Abbr <- Abbr[-Del,]
  }
##
## 5. Names
##
  Ab <- Abbr
  nc <- ncol(Abbr)
  nn <- length(Names)
  if(nn == nc){
      names(Abbr) <- Names
  } else {
      warning('Number of columns of the table read = ', nc,
              ' != length(Names) = ', nn, '; ignoring Names')
  }
##
## 6.  Clean
##
  if(clean){
      for(i in 1:nc){
          si <- subNonStandardCharacters(Abbr[,i])
#         delete leading "_"
          si2 <- sub('^_', '', si)
#         delete trailing footnote
          si3 <- strsplit(si2, '\\[')
          si4 <- sapply(si3, function(y){
              if(length(y)>0)y[1] else ''
          } )
          Abbr[, i] <- si4
      }
  }
##
## 7.  Done
##
  Abbr
}
