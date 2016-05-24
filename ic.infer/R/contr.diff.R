contr.diff <- function (n, contrasts=TRUE){
    if (length(n) <= 1) {
        if (is.numeric(n) && length(n) == 1 && n > 1) 
            levels <- 1:n
        else stop("invalid choice for n in contr.diff")
    }
    else levels <- n
    lenglev <- length(levels)

    ## definition of contrast matrix
    if (contrasts){
      if (lenglev==2) cont <- matrix(c(0,1),ncol=1)
       else {
          cont <- diag(1, lenglev)
          for (i in 2:(lenglev))
            for (j in 1:(i-1))
              cont[i,j] <- 1
                  cont <- cont[,(2:lenglev)]
      }
       rownames(cont) <- levels
       colnames(cont) <- paste(levels[2:lenglev],"-",levels[1:(lenglev-1)],sep="")
    }
    else{
        cont <- array(0, c(lenglev, lenglev), list(levels, levels))
        cont[col(cont) == row(cont)] <- 1
    }
    cont
}
