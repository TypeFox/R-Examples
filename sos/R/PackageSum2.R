PackageSum2 <- function(x,
      fields=c("Title", "Version", "Author", "Maintainer", "Packaged",
          'helpPages', 'vignette'), lib.loc=NULL, ...){
  UseMethod('PackageSum2')
}

PackageSum2.findFn <- function(x,
      fields=c("Title", "Version", "Author", "Maintainer", "Packaged",
          'helpPages', 'vignette'), lib.loc=NULL, ...){
  PackageSum2(attr(x, 'PackageSummary'), fields, lib.loc, ...)
}

PackageSum2.list <- function(x,
      fields=c("Title", "Version", "Author", "Maintainer", "Packaged",
          'helpPages', 'vignette'), lib.loc=NULL, ...){
  PackageSum2(x$PackageSummary, fields, lib.loc, ...)
}

PackageSum2.data.frame <- function(x,
      fields=c("Title", "Version", "Author", "Maintainer", "Packaged",
          'helpPages', 'vignette'), lib.loc=NULL, ...){
##
## 1.  Create character matrix for fields
##
  nf <- length(fields)
  nx <- nrow(x)
#  xout <- x
#  for(ic in seq(1, length=nf))xout[[fields[ic]]] <- rep('', nx)
  xP <- as.character(x$Package)
  xnew <- matrix('', nx, nf, dimnames=list(xP, fields))
##
## 2.  installed packages?
##
  instPkgs <- .packages(TRUE)
##
## 3.  Get packageDescription for each package
##
  for(ip in seq(1, length=nx)){
    if(xP[ip] %in% instPkgs){
      pkgDesci <- packageDescription(x$Package[ip], lib.loc=lib.loc)
      pkgHelp <- try(help(package=x$Package[ip], lib.loc=lib.loc,
                          help_type='text'))
      if(class(pkgHelp) != 'try-error'){
        for(ic in seq(1, length=nf)){
          if(fields[ic] == "Packaged"){
            if(is.null(pkgDesci$Packaged)){
              pkgd <- (strsplit(pkgDesci$Built, ';')[[1]][3])
            } else
              pkgd <- (strsplit(pkgDesci$Packaged, ';')[[1]][1])
            xnew[ip, ic] <- pkgd
#          xout$Packaged[ip] <- pkgd
          } else
            if(fields[ic] %in% names(pkgDesci)){
              xnew[ip, ic] <- pkgDesci[[fields[ic]]]
#            xout[ip, fields[ic]] <- pkgDesci[[fields[ic]]]
            } else {
              if(fields[ic] == 'helpPages'){
                helpInfo2 <- pkgHelp$info[[2]]
                nhr <- (length(helpInfo2) -
                    sum(substring(helpInfo2, 1, 1) == ' '))
                xnew[ip, ic] <- nhr
              } else {
                if(fields[ic] == 'vignette') {
#  The following says nrow(vinfo) has zero length in vignette
#  using R version 3.0.1 Patched (2013-06-21 r63003)
#  using platform: i386-pc-solaris2.10 (32-bit) in a CRAN test
#  For x$Package[ip] == 'base', dim(vignette(.)$results) = 0  4
#                  vig <- (vignette(package=x$Package[ip])$results)
                  xPi <- x$Package[ip]
                  VIG <- vignette(package=x$Package[ip])
                  vig <- VIG$results
                  if(nrow(vig)>1){
                    clps <- paste(vig[, 'Item'], collapse=', ')
                    xnew[ip, ic] <- paste(nrow(vig), clps, sep=':  ')
                  } else if(nrow(vig)>0)
                      xnew[ip, ic] <- vig[ 1, 'Title']
                }
              }
            }
        }
      }
    }
  }
##
## 4.  Parse Packaged and update Date
##
  xnew. <- as.data.frame(xnew, stringsAsFactors=FALSE)
  if(nx>0){
    pkgd <- xnew[, 'Packaged']
    nch <- nchar(pkgd)
    dateCh <- pkgd[nch>0]
    nd <- length(dateCh)
    dateC <- dateCh
    if(nd>0){
      bl1 <- (substring(dateCh, 1,1)==' ')
      dateCh[bl1] <- substring(dateCh[bl1], 2)
      num. <- regexpr('[0123456789]', dateCh)
      dateC[num.<9] <- substring(dateCh[num.<9], 1, 10)
      mo <- substring(dateCh[num.>8], 5)
      dd <- as.Date(mo, '%b %d %H:%M:%S %Y')
      dateC[num.>8] <- as.character(dd)
    }
#    pkgd[nch>0] <- dateC
#    xnew.$Packaged <- pkgd
    x$Date[nch>0] <- dateC
  }
  xnew.$Packaged <- NULL
##
## 5.  Done
##
  fixSpace <- function(x){
      x. <- gsub('\n', ' ', x)
      nx <- nchar(x.)
      repeat{
          x. <- gsub('  ', ' ', x.)
          nx2 <- nchar(x.)
          if(all(nx2==nx))break
          nx <- nx2
      }
      x.
  }
  x. <- cbind(x, xnew.)
  rownames(x.) <- 1:nx
  x.$Title <- fixSpace(x.$Title)
  x.$Author <- fixSpace(x.$Author)
  x.$Maintainer <- fixSpace(x.$Maintainer)
  x.$vignette <- fixSpace(x.$vignette)
  x.
}
