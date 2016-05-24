mergeVote <- function(x, vote, Office="House", vote.x,
                      check.x=TRUE){
##
## 1.  parse vote.x
##
  nameOfx <- deparse(substitute(x))
  nameOfVote <- deparse(substitute(vote))
  nx <- nrow(x)
  nv <- nrow(vote)
  nmx <- names(x)
  nmv <- names(vote)
  votey <- grep('vote', nmv, value=TRUE)
  if(length(votey)<1)
      stop('No vote column found in the vote data.frame = ',
           deparse(substitute(vote)))
#
  if(missing(vote.x)){
      vote.x <- grep('vote', names(x), value=TRUE)
      if(length(vote.x)<1)vote.x <- votey
  }
  if(!(vote.x %in% names(x)))
      x[, vote.x] <- rep('notEligible', nx)
##
## 2.  Office
##
  if(!('Office' %in% nmv))
      vote <- cbind(vote, Office=Office)
##
## 3.  keys
##
  lnmx <- tolower(nmx)
  lnmv <- tolower(nmv)
  surnmx <- nmx[grep('surname', lnmx)]
  surnmv <- nmv[grep('surname', lnmv)]
  givenx <- nmx[grep('givenname', lnmx)]
  givenv <- nmv[grep('givenname', lnmv)]
  stx <- nmx[grep('state', lnmx)]
  stv <- nmv[grep('state', lnmv)]
  distx <- nmx[grep('district', lnmx)]
  distv <- nmv[grep('district', lnmv)]
  keyx <- paste(x$Office, x[[surnmx]], sep=":")
  keyv <- paste(vote$Office, vote[[surnmv]], sep=":")
  keyx2 <- paste(keyx, x[[givenx]], sep=":")
  keyv2 <- paste(keyv, vote[[givenv]], sep=':')
  keyx. <- paste(x$Office, x[[stx]], x[[distx]], sep=":")
  keyv. <- paste(vote$Office, vote[[stv]], vote[[distv]], sep=":")
##
## 4.   record votes
##
  vote.notFound <- integer(0)
  voteFound <- rep(0, nv)
  for(iv in 1:nv){
      jv <- which(keyx == keyv[iv])
      if(length(jv)<1){
          jv <- which(keyx. == keyv.[iv])
          if(length(jv)!=1)
              vote.notFound <- c(vote.notFound, iv)
      }
      if(length(jv)>1){
          jv <- which(keyx2 == keyv2[iv])
          if(length(jv)!=1)
              jv <- which(keyx.==keyv.[iv])
#              vote.notFound <- c(vote.notFound, keyv2[iv])
          if(length(jv)!=1){
              vote.notFound <- c(vote.notFound, iv)
          }
      }
      if(length(jv)==1) {
          x[jv, vote.x] <- as.character(vote[iv, votey])
          voteFound[iv] <- jv
      }
  }
##
## 5.  Check
##
  if(check.x){
      Votex <- which(x[, vote.x] != 'notEligible')
      oops <- which(!(Votex %in% voteFound))
      if(length(oops)>0){
          print(x[oops,])
          stop('People found voting in x = ', nameOfx[1],
               '\n  not found in the data.frame vote = ',
               nameOfVote[1],
               '\n  look for and fix the error(s) printed above.')
      }
  }
##
##  Done
##
  x[, vote.x] <- factor(x[, vote.x])
  if((no <- length(vote.notFound))>0){
      cat(no, 'rows of vote not found:\n')
      print(vote[vote.notFound,] )
      stop('Unable to find vote in x')
  }
  x
}

