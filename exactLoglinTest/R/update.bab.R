##this is the workhorse program for the
##booth and butler method
##it requires the input to be of the form
##constructed by mcexact
update.bab <- function(object, ...){
  bab(object, ...)
}

bab <-  function(args, nosim = NULL, maxiter = NULL, savechain = FALSE){
  if (!is.null(nosim)) args$nosim <- nosim
  if (!is.null(maxiter)) args$maxiter <- maxiter
  if (args$nosim >  args$maxiter) {
    warning("update.bab, nosim > maxiter, setting maxiter = nosim")
    args$maxiter <- args$nosim
  }
  if (is.null(args$startiter)) args$startiter <- 1
  if (is.null(args$sumdw)) args$sumdw <- 0
  if (is.null(args$sumdwsq)) args$sumdwsq <- 0
  if (is.null(args$sumw)) args$sumw <- 0
  if (is.null(args$sumwsq)) args$sumwsq <- 0
  if (savechain){
    nocol <- length(args$stat(rowlabels = TRUE))
    args$chain <- matrix(0, args$nosim, nocol + 1)
    colnames(args$chain) <- c(args$stat(rowlabels = TRUE),  "log imp weight")
  }
  else args$chain <- NULL
  perpos <- 0

  i <- args$startiter
  j <- 0
  while ((i - args$startiter < args$nosim) & (j < args$maxiter)){
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
       perpos <- perpos + 1
       d <- args$stat(y = y.new, mu = args$mu.hat, rowlabels = FALSE)
       ##importance weights on the log scale
       w <- args$dens(y.new) - rounded.tprob(y1.new.permute, conde1.permute, diag(condv1.permute), args$tdf)
       ##the following subtracts off a constant from
       ##all of the importance weights
       ##the constant is the weight of the first
       ##simulated table
       if (savechain) args$chain[i - args$startiter + 1,] <- c(d, w)
       
       if (i == 1) args$impconst <- w
       w <- exp(w - args$impconst)
       ##the following are the partial sums required for the
       ##importance sampling estimate
       args$sumdw <- args$sumdw + (d >= args$dobs) * w
       args$sumdwsq <- args$sumdwsq + (d >= args$dobs) * w ^ 2
       args$sumw <- args$sumw  + w
       args$sumwsq <- args$sumwsq + w ^ 2
       i <- i + 1
     }
   }
  if (i == args$startiter)
    warning("No valid tables found")
  else if (savechain){
    args$chain <- args$chain[1 : (i - args$startiter),]
  }
 # if ((i - args$startiter) < (args$nosim - 1))
 #  warning("Maximum iterations reached yet desired number of simulated tables not attained.")
  args$startiter <- i
  theta <- args$sumdw / args$sumw
  setheta <- sqrt((1 - 2 * theta) * args$sumdwsq + (theta ^ 2) * args$sumwsq) / (i - 1);
  args$phat <- theta
  args$mcse <- setheta
  args$perpos <- perpos / min(args$nosim, args$maxiter)
  return(args)
}
