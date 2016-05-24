simtable.bab <-  function(args, nosim = NULL, maxiter = NULL){
  if (!is.null(nosim)) args$nosim <- nosim
  if (!is.null(maxiter)) args$maxiter <- maxiter
  if (args$nosim >  args$maxiter) {
    warning("update.bab, nosim > maxiter, setting maxiter = nosim")
    args$maxiter <- args$nosim
  }
  nocol <- length(args$y)
  chain <- matrix(0, args$nosim, nocol + 1)
  i <- 1
  j <- 0
  while ((i <= args$nosim) & (j < args$maxiter)){
     j <- j + 1
     shuffle <- sample(1 : args$n1)
     conde1.permute <-  args$conde1[shuffle]
     ##this is equal to P %*% condv1 %*% t(P)
     condv1.permute <-  args$condv1[shuffle, shuffle]
     ##get and unshuffle y1.new
     ##note conde1.permute and condv1.permute
     ##are now the sequential means and variances
     y1.new.permute <- .Call("multinormfull",
                             conde1.permute,
                             condv1.permute,
                             args$tdf,
                             PACKAGE="exactLoglinTest")
     y1.new <- y1.new.permute[order(shuffle)]
     y2.new <- args$x2invt %*% (args$s - t(args$x1) %*% y1.new)
     ##though technically y.new has to be an interger we
     ##coerce here since the calculation
     ##of y.new is done as double
     y.new <- round(c(y1.new, y2.new))

     # HJ
     N.na <- sum(is.na(y.new))
     if (i==1 & N.na > 0)
	warning("i = ", i, ": 'y.new' has ", N.na, " NAs")

     if (N.na==0 & all(y.new >= 0)){
       d <- y.new
       ##importance weights on the log scale
       w <- args$dens(y.new) - rounded.tprob(y1.new.permute,
                                             conde1.permute,
                                             diag(condv1.permute),
                                             args$tdf)
       ##the following subtracts off a constant from
       ##all of the importance weights
       ##the constant is the weight of the first
       ##simulated table
       chain[i,] <- c(d, w)
       i <- i + 1
     }
   }
  if (i == 1)
    warning("No valid tables found")
  else {
    rval <- cbind(chain[1 : (i - 1), order(args$ord)], chain[1 : (i-1), nocol +1])
    #return(chain[1 : (i - 1),order(args$ord)])
    return(rval)
  }
}







