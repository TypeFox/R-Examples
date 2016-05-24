.SPNET.COLLAPSE.SEP <- ';'

.extract.multiple.strings <- function(x, sep = .SPNET.COLLAPSE.SEP) {
  out <- c()
  for (k in x) {
    spl <- strsplit(x = k, split = sep)[[1]]
    spl <- spl[nzchar(spl)] ## remove "" strings
    spl <- sub("^ +", "", spl)  ## remove leading spaces
    spl <- sub(" +$", "", spl)  ## remove trailing spaces
    out <- c(out, spl)
  }
  return(out)
  
  #   sub("[[:space:]]+$", "", str) ## white space, POSIX-style
  #   sub("\\s+$", "", str, perl = TRUE) ## Perl-style white spac
}

# a1 <- 'President'
# a2 <- ';President'
# a3 <- 'President;'
# a4 <- "President;Team leader"
# a5 <- "President ; Team leader"
# .extract.multiple.strings(a1)
# .extract.multiple.strings(a2)
# .extract.multiple.strings(a3)
# .extract.multiple.strings(a4)
# .extract.multiple.strings(a5)
# .extract.multiple.strings(c(a4, a5))

.expand.multiple.names <- function(x, sep = .SPNET.COLLAPSE.SEP) {
  nam <- names(x)  
  out <- c()
  out.nam <- character(0)
  for (k in 1:length(x)) {
    nam.temp <- .extract.multiple.strings(nam[k])
    for (i in 1:length(nam.temp)) {
      out <- c(out, x[k])
      out.nam <- c(out.nam, nam.temp[i])
    }
  }
  names(out) <- out.nam
  return(out)
}


# b1 <- c("Président;Chef de groupe" = 2, "Chef de groupe" = 4, "Porteur du projet" = 5) 
# .expand.multiple.names(b1)

setReplaceMethod(
  f ="[",
  signature ="SpatialNetwork",
  definition = function(x, i, j, value){
    if (length(value) > 1) {
      value <- paste(value, collapse = .SPNET.COLLAPSE.SEP)
    }
    if(missing(i)) {
      if(missing(j)) {
        stop("error: you have to specifie i or j")
      } else {
        if (inherits(j, 'character')) {
          j <- match(j, names(x))
        }
        stopifnot(all(!is.na(j)))
        stopifnot(min(j) >= 1 && max(j) <= ncol(x))
        for(k in j){
          x@.Data[[k]][i] <- value
        }
        return(x)
      }
    } else { # i not missing
      if (missing(j)){
        j <- 1:ncol(x)
      }
      if (inherits(j, 'character')) {
        j <- match(j, names(x))
      }
      stopifnot(all(!is.na(j)))
      stopifnot(min(j) >= 1 && max(j) <= ncol(x))
      for(k in j){
        x@.Data[[k]][i] <- value
      }
      return(x)
    }
  }
)

# a <- spnet::spnet.create(data.frame('NODE' = 1:4, 'POSITION' = 8:11))
# a[1,'NODE'] <- 6
# a
# a[,'NODE'] <- c(1,3)
# a[1,] <- c(2,3)


.position.multiple.symbols <- function(x, n = 1, cex = 1, space = 0.5){
  if(n==1) return(matrix(x, ncol = 2))
  
  x.x <- x[1]
  x.y <- x[2]
  
  x.all <- x.x + seq(from = 0, to = space*cex*(n-1), by = space*cex)
  
  if(n%%2 == 0) { ## even
    x.all <- x.all - space*cex*((n/2 - 1) + 0.5)
  }
  if(n%%2 == 1) { ## odd
    x.all <- x.all - space*cex*((n-1)/2)
  }
  
  y.all <- rep(x.y, n)
  
  out <- matrix(c(x.all, y.all), nrow=n, byrow=FALSE)
  return(out)
}

# .position.multiple.symbols(c(0,1), n = 1)
# .position.multiple.symbols(c(0,1), n = 2)
# .position.multiple.symbols(c(0,1), n = 2, cex = 2)
# .position.multiple.symbols(c(0,1), n = 3)
# .position.multiple.symbols(c(0,1), 3, cex = 2)
