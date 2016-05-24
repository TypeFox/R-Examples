match.data.frame <- function(x, y, by, by.x=by, by.y=by,
                             grep., split, sep=':'){
##
## 1.  Check lengths of by, grep., split
##
  if((missing(by.x) || missing(by.y)) && missing(by)){
      by <- names(x)
  }
#
  kx <- length(by.x)
  nx <- nrow(x)
  if(kx<1){
      warning('length(by.x)==0;  nothing to match.  Returning NAs')
      return(rep(NA, nx))
  }
  chk.x <- which(!(by.x %in% names(x)))
  if(length(chk.x)>0){
      stop('by.x not in names(x);  first error = ',
           by.x[chk.x[1]])
  }
 #
  if((ky <- length(by.y))!= kx){
      stop('length(by.x) = ', kx, ' != length(by.y) = ', ky)
  }
  chk.y <- which(!(by.y %in% names(y)))
  if(length(chk.y)>0){
      stop('by.y not in names(y);  first error = ',
           by.y[chk.y[1]])
  }
#
  if(missing(grep.)){
      grep. <- rep(NA, kx)
  } else {
      if((kg <- length(grep.)) == 1){
         grep. <- rep(grep., kx)
     } else if(kg != kx)
         stop('length(by.x) = ', kx, ' != length(grep.) = ', kg)
  }
#
  if(missing(split)){
      split <- c(NA, ' ')[1+!is.na(grep.)]
  } else {
      if((ks <- length(split))==1){
          split <- rep(split, kx)
      } else if(ks != kx)
          stop('length(by.x) = ', kx, ' != length(split) = ', ks)
  }
##
## 2.  fullMatch
##
  oops <- which(is.na(grep.) != is.na(split))
  if(length(oops)>0){
      stop('grep. cannot be NA when split is not and vice versa;',
           '  first error in position number ', oops[1],
           ' where grep. = ', grep[oops[1]],
           ' and split = ', split[oops[1]])
  }
#
  fullMatch <- (is.na(grep.) & is.na(split))
  afM <- any(fullMatch)
  if(afM){
      matchx <- c(as.list(x[, by.x[fullMatch], drop=FALSE]), sep=sep)
      keyfx <- with(x, do.call(paste, matchx))
#
      matchy <- c(as.list(y[, by.y[fullMatch], drop=FALSE]), sep=sep)
      keyfy <- with(y, do.call(paste, matchy))
  }
  anf <- any(!fullMatch)
  if(anf){
      parMatch <- ((1:kx)[!fullMatch])
      x. <- mapply(function(z, s){
          z. <- strsplit(z, s)
          sapply(z., '[', 1)
      }, x[, by.x[parMatch], drop=FALSE],
                   split[parMatch])
#      y. <- lapply(by.y[parMatch], function(z){
#          strsplit(z, split)
#      } )
  }
##
## 3.  iterate over each row of x
##
  xyi <- rep(NA, nx)
  ny <- nrow(y)
  iy. <- 1:ny
  for(ix in 1:nx){
      if(afM){
          iy <- which(regexpr(keyfx[ix], keyfy, fixed=TRUE)>0)
          if((kiy <- length(iy))>0){
              if(kiy<2) {
                  xyi[ix] <- iy
              }
          }
      } else iy <- iy.
      if(anf & (kiy>1)){
          for(j in parMatch){
              xij <- x.[ix, by.x[j]]
              xiju <- do.call(grep.[j], list(xij, y[iy, by.y[j]],
                              fixed=TRUE) )
              if((kiju <- length(xiju))<1) next
#             match not found
              if(kiju<2){
                  xyi[ix] <- iy[xiju]
                  next
#                 unique match
              } else iy <- iy[xiju]
#             multiple matches;  continue
          }
      }
  }
##
## 4.  done
##
  xyi
}
