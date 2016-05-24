readCookPVI. <- function(url.=
"https://en.wikipedia.org/wiki/Cook_Partisan_Voting_Index",
      UShouse=readUShouse(), USsenate=readUSsenate(), ...){
##
## 1.  readCookPVI()
##
  CookPVI <- readCookPVI(url.)
##
## 2.  merge with UShouse
##
  keyPVI <- with(CookPVI$House, paste(State, District, sep='.'))
  houseKey <- with(UShouse, paste(state, district, sep='.'))
  rownames(UShouse) <- houseKey
  House <- cbind(UShouse[keyPVI, ], CookPVI$House[, c('PVInum', 'PVIchar')])
##
## 3.  merge with USsenate
##
  CookSenate <- CookPVI$Senate
  rownames(CookSenate) <- CookSenate$State
  Senate <- cbind(USsenate, CookSenate[USsenate$State, -1])
  rownames(Senate) <- rownames(USsenate)
##
## 4.  Done
##
  list(House=House, Senate=Senate)
}

readCookPVI <- function(url.=
"http://en.wikipedia.org/wiki/Cook_Partisan_Voting_Index"){
##
## 1.  download content
##
#  library(RCurl)
  Start <- paste(date(), ': readCookPVI(',
                 url., ')', sep='')
  cat(Start)
  startTime <- proc.time()
  Url. <- try(RCurl::getURL(url.))
  et <- max(proc.time()-startTime, na.rm=TRUE)
  Read <- paste('|', nchar(Url.), 'bytes read in',
                round(et, 2), 'seconds\n')
  cat(Read)
  if(class(Url.)=='try-error'){
      cat(Url.)
      stop()
  }
##
## 2.  readHTMLTable
##
#  library(XML)
  Wikitbls <- XML::readHTMLTable(Url., stringsAsFactors=FALSE)
##
## 3.  Find House and Senate tables
##
  len <- sapply(Wikitbls, length)
  Wikitbs <- Wikitbls[len>0]
  ns <- sapply(Wikitbs, dim)
  house <- which((434<ns[1,]) & (ns[1,]<450))
  if((nh <- length(house))!=1){
      stop('URL source changed: Found ', nh,
           ' tables with between 434 and 450 rows.  oops.')
  }
  House <- Wikitbs[[house]]
  district <- sub('st$|nd$|rd|th', '', House$District)
  Dct <- sub('At-large', '0', district)
  House$District <- Dct
#
  senate <- which(ns[1,]==50)
  if((nsen <- length(senate)) != 1){
      stop('URL source changed:  Found ', nsen,
           ' tables with 50 rows.  oops.')
  }
  Senate <- Wikitbs[[senate]]
##
## 4.  Parse PVI columns
##
  pvi <- function(x){
      x. <- strsplit(x, ' !')
      PVIn <- sapply(x., '[', 1)
      PVInum <- as.numeric(PVIn)
      PVIchar <- sapply(x., '[', 2)
      list(PVInum=PVInum, PVIchar=PVIchar)
  }
  House. <- cbind(House[c('State', 'District')],
                  pvi(House[['PVI']]),
                  stringsAsFactors=FALSE)
  House.$PartyOfRepresentative <- factor(House[['Party of\nRepresentative']])
#
  Senate. <- cbind(Senate['State'], pvi(Senate[['PVI']]),
                   stringsAsFactors=FALSE)
  Senate.$PartyOfGovernor <- factor(Senate[['Party of\nGovernor']])
  Senate.$PartyInSenate <- factor(Senate[['Party\nin Senate']])
  houseBal <- pvi(Senate[['House\nbalance']])
  Senate.$houseBalanceNum <- houseBal[[1]]/10
  Senate.$houseBalanceChar <- houseBal[[2]]
##
## 8.  Done
##
  list(House=House., Senate=Senate.)
}

