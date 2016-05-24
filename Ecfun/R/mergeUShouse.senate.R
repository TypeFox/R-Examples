mergeUShouse.senate <- function(x, UScongress=UShouse.senate(),
                                newrows='amount0',
        default=list(member=FALSE, amount=0, vote="notEligible",
                     incumbent=TRUE) ){
##
## 0.  Check district
##
  noDist <- sum(is.na(x$district))
  if(noDist>0)
      stop(noDist, ' NAs in x$district')
##
## 1.  keys
##
#  X <- x
#  x$district[x$district=='0'] <- 'At Large'
#
  keyx <- with(x, paste(Office, state, district, sep=":"))
  keyy <- with(UScongress, paste(Office, state, district, sep=":"))
##
## 2.  notx
##
  notx <- !(keyy %in% keyx)
  huh <- !(keyx %in% keyy)
  if((nhuh <- sum(huh))>0){
      cat(nhuh, 'Districts in x not found in UScongress;  the first is:\n')
      print(x[huh,][1,])
      stop('District coding problem in x')
  }
  notx. <- (notx & !UScongress$nonvoting)
#
  Y <- UScongress[notx., ]
  keyy. <- keyy[notx.]
##
## 3.  Add default columns to Y
##
  nd <- length(default)
  nmx <- names(x)
  for(id in 1:nd){
      found <- regexpr(names(default)[id], tolower(nmx))
      for(f in nmx[found>0]){
          Y[, f] <- default[id]
      }
  }
##
## 4.  newrows
##
  if(!(newrows %in% nmx)){
      x <- cbind(x, FALSE)
      nx <- ncol(x)
      names(x)[nx] <- newrows
      nmx <- names(x)
  }
  Y[, newrows] <- TRUE
##
## 4.  rbind
##
  nmx.Y <- (nmx %in% names(Y))
  nmxY <- which(!nmx.Y)
  if((Noops <- length(nmxY))>0){
      warning(Noops, " column(s) of x not in Y;  first = ",
              nmx[nmxY], '.  Discarding.')
  }
  nmY.x <- nmx[nmx.Y]
  xY <- rbind(x[nmx.Y], Y[nmY.x])
##
## 5.  Replace 'Democrat' with 'Democratic' in xY$Party
##
  Pty <- xY$Party
  Pty[Pty=='Democrat'] <- 'Democratic'
  xY$Party <- factor(Pty)
##
## 6.  incumbent?
##
  if('incumbent' %in% names(xY)){
      oops <- is.na(xY$incumbent)
      keyo <- with(xY[oops, ], paste(Office, state, district, sep=":"))
#     see notx above
#      oddDist <- which(!(keyo %in% keyy))
#      if((nod <- length(oddDist))>0){
#          msg <- paste(nod, 'district(s) in x not in UScongress:')
#          cat(msg, '\n')
#          print(xY[oddDist, ])
#          warning(msg)
#      }
      kexY <- which(keyo %in% keyy)
      match1 <- function(a,b){
          a1 <- strsplit(a, ' ')[[1]][1]
          b. <- strsplit(b, ' ')[[1]]
          a1 %in% b.
      }
      for(ix in kexY){
          USci <- UScongress[keyy == keyo[ix], ]
          xYix <- xY[oops,][ix,]
          matchSur <- match1(xYix$surname, USci$surname)
          matchGiv <- match1(xYix$givenName, USci$givenName)
          xY[oops,][ix, 'incumbent'] <- (matchSur & matchGiv)
      }
  }
##
## 7.  done
##
  xY
}
