parseName <- function(x, surnameFirst=(median(regexpr(',', x))>0),
          suffix=c('Jr.', 'I', 'II', 'III', 'IV', 'Sr.', 'Dr.', 
                   'Jr', 'Sr'),
          fixNonStandard=subNonStandardNames, 
          removeSecondLine=TRUE, 
          namesNotFound="attr.replacement", ...){
##
## 1.  length(x)<1?
##
  nx <- length(x)
  if(nx<1){
      warning("length(x) == 0")
      mat0 <- matrix(character(0), nrow=0, ncol=2,
         dimnames=list(NULL, c("surname", "givenName")))
      return(mat0)
  }
  x00 <- x
# Force eval of default if(missing(surnameFirst))
# else x changes and surnameFirst may change 
  miSNF <- missing(surnameFirst)
  surnameFirst <- surnameFirst
  x <- x01 <- tis::stripBlanks(x00)
##
## 2.  secondLine?
##
  if(removeSecondLine){
    X2. <- strsplit(x01, '\n')
    x2 <- sapply(X2., function(xx){
      xx2 <- xx[-1]
      paste(xx2, collapse='\n')
    } ) 
    no2 <- sapply(X2., length)
    x2[no2<2] <- NA
    x <- sapply(X2., '[', 1)
    nch0 <- (nchar(x01)<1)
    x[nch0] <- x01[nch0] 
    x <- tis::stripBlanks(x)
  } else x2 <- NULL
##
## 3.  Drop "(AL)", etc.
##
  dropEndParen <- function(x){
      endParen <- grep(')$', x)
      if(length(endParen)>0){
          x.endP <- x[endParen]
          nch <- nchar(x.endP)
          nch3 <- pmax(1, nch-3)
          openP <- substring(x.endP, nch3, nch3)
          nch3. <- ifelse(openP=='(', pmax(1, nch3-1), nch)
          x.woP <- substring(x.endP, 1, nch3.)
          x[endParen] <- x.woP
      }
#      if(require(tis)){
#          x <- stripBlanks(x)
#      } else {
#          warning('need stripBlanks{tis} to delete leading and trailing',
#                  ' blanks;  not available')
#      }
      x
  }
  x0 <- x
  x <- dropEndParen(x)
##
## 4.  surnameFirst
##
  if(surnameFirst){
      Sep <- regexpr(',', x)
      oops <- which(Sep<0)
      if(miSNF & (no <- length(oops))>0){
          err <- paste(no, ' of ', length(x), ' elements of x ',
                       'are NOT in surnameFirst format.  ',
                       'The first is ', x[oops[1]],
               '.  Assume they are in (givenName surname) format.')
          warning(err)
#          fix0 <- parseName(x[oops])
#          fix1 <- paste(fix0[, 2], fix0[, 1], sep=', ')
#          x[oops] <- fix1
      }
      sur. <- substring(x, 1, Sep-1)
#     drop (AL), etc.
      sur <- dropEndParen(sur.)
#
      giv <- substring(x, Sep+1)
#      if(require(tis)){
#          giv <- stripBlanks(giv)
#      } else {
#          warning('need stripBlanks{tis} to delete leading and trailing',
#                  ' blanks;  not available')
#      }
      Sur <- fixNonStandard(sur, namesNotFound=namesNotFound, ...)
      Sur. <- strsplit(Sur, ' ')
      look4suf <- sapply(Sur., tail, n=1)
      suf <- (look4suf %in% suffix)
      surname <- rep('', nx)
      lx <- sapply(Sur., length)-suf
      for(ix in 1:nx){
          xi <- Sur.[[ix]][seq(length=lx[ix])]
          surname[ix] <- paste(xi, collapse=' ')
      }
#
      Giv <- fixNonStandard(giv, namesNotFound=namesNotFound, ...)
      Giv[suf] <- paste(Giv[suf], look4suf[suf], sep=', ')
      out <- cbind(surname, Giv)
      colnames(out) <- c('surname', 'givenName')
      nNF <- c(attr(Sur, 'namesNotFound'), 
               attr(Giv, 'namesNotFound'))
      if(!is.null(nNF))attr(out, 'namesNotFound') <- nNF
      return(out)
  }
##
## 5.  strsplit on either blank or comma
##
  x. <- strsplit(x, ' |,')
##
## 6.  Suffix?
##
  look4suf <- sapply(x., tail, n=1)
  suf <- (look4suf %in% suffix)
##
## 7.  parse
##
  givenName <- surname <- rep('', nx)
  lx <- sapply(x., length)-suf
# Check for ending comma 
  endComma <- regexpr(',$', x)
  for(ix in 1:nx){
    if(lx[ix]<1) next
    if(lx[ix]<2) {
#   only one name 
      if(endComma[ix]>0){
        surname[[ix]] <- x.[[ix]]
      } else {
        givenName[[ix]] <- x.[[ix]]
      }
    } else {
      xi <- x.[[ix]][1:lx[ix]]
#      ni <- length(xi)
      ni0 <- nchar(xi)
      xi. <- xi[ni0>0]
      surname[ix] <- tail(xi., 1)
      ni <- length(xi.)
      givenName[ix] <- paste(head(xi., ni-1), collapse=' ')
      if(suf[ix])
          givenName[ix] <- paste(givenName[ix], look4suf[ix],
                                 sep=', ')
    }
  }
##
## 8.  fixNonStandard
##
  Sur <- fixNonStandard(surname, namesNotFound=namesNotFound, ...)
  Giv <- fixNonStandard(givenName, namesNotFound=namesNotFound, ...)
##
## 9.  Done
##
  out <- cbind(surname=Sur, givenName=Giv)
  nchx2 <- nchar(x2)
  nchx2[is.na(x2)] <- 0 
  if(any(nchx2>0))attr(out, 'secondLine') <- x2
  nNF <- c(attr(Sur, 'namesNotFound'), 
           attr(Giv, 'namesNotFound'))
  if(!is.null(nNF))attr(out, 'namesNotFound') <- nNF
  out
}

