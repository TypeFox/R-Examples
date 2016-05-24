Newdata <- function(data, x, n, na.rm=TRUE){
##
## 1.  check data
##
  DF <- inherits(data, 'data.frame')
  Dat <- as.data.frame(data)
  vars <- names(Dat)
##
## 2.  check x
##  
  if(missing(x)){
    x <- names(Dat)
  } else {
    if(is.numeric(x)){
      x <- vars[x]
    } else {
      oops <- which(!(x %in% vars))
      if(length(oops)>0){
        msg0 <- paste('Some elements of x are not', 
                  'in colnames(data):')
        if(length(oops)>4){
          msg1 <- paste(msg0, paste(oops[1:4], collapse=', '), 
                        '...')
        } else msg1 <- paste(msg0, paste(oops, collapse=', '))
        msg2 <- paste(msg1, 'not in colnames(data) =')
        if(length(vars)>4){
          msg3 <- paste(msg2, paste(vars[1:4], collapse=', '), 
                        '...')
        } else {
          msg3 <- paste(msg2, paste(vars, collapse=', '))
        }
        stop(msg3)
      }
    }
  }
##
## 3.  check n
##
  nx <- length(x)
  if(missing(n)){
    n <- rep(2, nx) 
  } else {
    if(length(n)<1) {
      n <- rep(2, nx) 
    } else {
      if(length(n)<2){
        if(is.na(n)){
          n <- rep(2, nx)
        } else {
          n <- rep(n, nx)
        }
      } else {
        if(length(n) != length(x)){
          stop('length(n) = ', length(n), 
               ' != length(x) = ', length(x))
        }
      }
    }
  }
  if(is.null(names(n))){
    names(n) <- x
  } else {
    oopn <- which(!(names(n) %in% x))
    if(length(oopn)>0){
      stop('names(n) = ', names(n)[oopn[1]], 
           ' not matched in x; x[1] =', 
           x[1])
    }
  }
##
## 4.  loop over vars  
##
  Lat <- as.list(Dat)
  for(v in vars){
    if(v %in% x){
      if(n[v]<1){
        Lat[[v]] <- NULL 
      } else {
        v. <- data[, v]
        if(canbeNumeric(v.)){
          vn <- as.numeric(v.)
          if(n[v]<2){
            vm <- median(vn)
            attributes(vm) <- attributes(v.)
            Lat[[v]] <- vm 
          } else {
            vrng <- range(vn)
            vNew <- seq(vrng[1], vrng[2], length=n[v])
            attributes(vNew) <- attributes(v.)
            Lat[[v]] <- vNew 
          }
        } else {
          vt <- sort(table(v.), TRUE)
          if(length(vt) <= n[v]){
            Vt <- names(vt)
          } else {
            Vt <- names(vt)[1:n[v]]
          }
          iv <- 1:length(Vt)
          Vt. <- v.[iv]
          for(i. in iv){
            jv <- which(v. == Vt[i.])
            Vt.[i.] <- v.[jv[1]]
          }
          Lat[[v]] <- Vt.
        }
      }
    } else {
      v. <- data[, v]
      if(canbeNumeric(v.)){
        vn <- as.numeric(v.)
        vm <- median(vn)
        attributes(vm) <- attributes(v.)
        Lat[[v]] <- vm 
      } else {
        vt <- sort(table(v.), TRUE)
        Vt <- names(vt)[1]
        jv <- which(v. == Vt)
        Lat[[v]] <- v.[jv[1]]
      }
    }
  }
##
## 5.  create newDat 
##
  Lat$stringsAsFactors <- FALSE 
  newDat <- do.call(expand.grid, Lat)
##
## 6.  done 
##  
  newDat
}
