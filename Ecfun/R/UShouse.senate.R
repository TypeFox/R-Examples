UShouse.senate <- function(house=readUShouse(), senate=readUSsenate()){
##
## 1.  reformat house
##
  nh <- nrow(house)
#  cat('hrow(house) = ', nh, '\n')
  Party <- as.character(house$party)
  Party[Party=='D'] <- 'Democratic'
  Party[Party=='R'] <- 'Republican'
  np <- length(Party)
  if(np != nh){
      stop('nrow(house) = ', nh, '; length(Party) = ', np,
           ':  NOT equal')
  }
#
  hs <- with(house, data.frame(Office=factor(rep('House', nh)),
                   State=State, state=factor(toupper(state)),
                   district=district, Party=factor(Party),
                   surname=surname, givenName=givenName,
                   nonvoting=nonvoting, stringsAsFactors=FALSE) )
##
## 2.  reformat senate
##
  ns <- nrow(senate)
#  cat('nrow(senate) =', ns, '\n')
  hS <- rep('Senate', ns)
#
  sn <- with(senate, data.frame(Office=factor(hS),
                       State=State, state=state,
                       district=Class, Party=Party,
                       surname=surname,
                       givenName=givenName,
                       nonvoting=FALSE,
                       stringsAsFactors=FALSE) )
##
## 3.  rbind
##
  HS <- rbind(hs, sn)
  HS
}

