readUShouse <- function(url.="http://www.house.gov/representatives/",
   nonvoting=c('American Samoa', 'District of Columbia',
               'Guam', 'Northern Mariana Islands', 'Puerto Rico',
               'Virgin Islands'),
   fixNonStandard=subNonStandardNames, ...){
#, USstateAbbreviations=readUSstateAbbreviations() ){
##
## 1.  download content
##
#  library(RCurl)
  Start <- paste(date(), ': readUShouse(', url., ')', sep='')
  cat(Start)
  startTime <- proc.time()
  house.gov <- try(RCurl::getURL(url.))
  et <- max(proc.time()-startTime, na.rm=TRUE)
  Read <- paste('|', nchar(house.gov), 'bytes read in',
                round(et, 2), 'seconds\n')
  cat(Read)
  if(class(house.gov)=='try-error'){
      cat(house.gov)
      stop()
  }
##
## 2.  find "state"
##
#  st <- gregexpr('state_', house.gov)[[1]] # finds 75
  st <- gregexpr('\t<h2 id=\"', house.gov)[[1]]
#  substring(house.gov, st, st+28)
  st. <- sapply(st, function(x){
      hgi <- substring(house.gov, x)
      x2 <- (x+regexpr('\">', hgi))
  } )
  stCodes <- substring(house.gov, st+9, st.-2)
# 2014-12-06: state_al (Alabama representatives)... 
#    name_y (names beginning with y)
  st. <- strsplit(stCodes, "_")
  St. <- sapply(st., "[", 2)
##
## 3.  Convert to tables
##
#  library(XML)
  House.gov <- XML::readHTMLTable(house.gov, stringsAsFactors=FALSE)
  names(House.gov) <- stCodes
##
## 4.  Rbind tables with the same headers
##
  headers <- lapply(House.gov, names)
  h. <- sapply(headers, paste, collapse=":")
  H. <- table(h.)
  st2 <- St.[1:H.[1]]
#
  ns <- sapply(House.gov, nrow)
  St2 <- rep(st2, ns[1:length(st2)])
#
  nh <- length(H.)
  out <- vector(mode='list', length=nh)
  for(ih in 1:nh){
    sel <- (h.==names(H.)[ih])
    out[[ih]] <- do.call(rbind, House.gov[sel])
  }
  names(out) <- names(H.)
##
## 5.  If 2 tables, use the first but add the state names
##
  if(nh>2){
      warning(nh, " > 2 different tables found;  I'm confused")
      return(out)
  }
  n. <- sapply(out, nrow)
  if(n.[1] != n.[2]){
    warning('2 tables found with differing numbers of rows:  ', 
            paste(n., collapse=' and '), '; returning the larger')
#      stop('2 tables found with different numbers of rows')
#    return(out)
  }
#
  n2 <- max(n.)
  i2 <- which(n.==n2)[1]
  i1 <- 3-i2
  n1 <- n.[i1]
#  Dist <- out[[2]]$District
  Dist <- out[[i2]]$District
  D. <- strsplit(Dist, ' ')
  state <- sapply(D., function(x){
      nx <- length(x)
      Di <- (x[nx]=='District')
      x. <- paste(x[seq(length=nx-Di-1)], collapse=' ')
      x.
  } )
# but in the wrong order
#  State <- character(n.[1])
  State <- character(n2)
#  for(i. in 1:n.[1]){
  for(i. in 1:n2){
      j. <- which(out[[i2]]$Name==out[[i1]]$Name[i.])
      if((nj <- length(j.)) != 1){
        warning('for row ', i., ' found ', nj, ' matches')
      } else State[i.] <- state[j.]
  }
#
#  USPS <- grep('USPS', names(USstateAbbreviations))
#  USPScodes <- USstateAbbreviations[, USPS]
#  good <- (USPScodes!='')
#  stateAbbr <- USstateAbbreviations[good,]
#  USPScds <- USPScodes[good]
#  rownames(stateAbbr) <- USPScds
  ST2 <- toupper(St2)
#  which(!(ST2 %in% USPScds))
# state <- stateAbbr[ST2, "Name"]
#
# Out <- cbind(data.frame(State = State, state = ST2), out[[1]][1:5])
  Out <- cbind(data.frame(State=State, state=ST2),
               out[[i2]][1:5])
  Out$Party <- factor(Out$Party)
# Out$Committees <- out[[1]][["Committee Assignment"]]
#  Out$Committees <- out[[i2]][["Committee Assignment"]]
  Out$nonvoting <- (Out$State %in% nonvoting)
#
  surnm <- parseName(Out$Name, TRUE)
  O <- cbind(Out, as.data.frame(surnm, stringsAsFactors=FALSE))
##
## 6.  convert District to "district",
##     "At Large", "At-Large", "AtLarge" to "0"
##
  Dist <- which(names(O) %in% "District")
  if(length(Dist)!=1){
      print(names(O))
      stop('UShouse does not contain a unique column named "District"')
  }
  names(O)[Dist] <- "district"
  zero <- (O[, Dist] %in% c("At Large", "At-Large", "AtLarge"))
  O[zero, Dist] <- "0"
##
## 7.  rename Party -> party because it's "D" & "R",
##     and standard "Party" = "Democrat", "Republican", ...
##
  Pty <- which(names(O) %in% "Party")
  if(length(Dist)!=1){
      print(names(O))
      stop('UShouse does not contain a unique column named "Party"')
  }
  names(O)[Pty] <- "party"
##
## 8.  nonStandard?
##
  O$surname <- fixNonStandard(O$surname, ...)
  O$givenName <- fixNonStandard(O$givenName, ...)
##
## 8.  Done
##
  O
}
