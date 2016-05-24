###this is the workhorse program for the booth and butler method
###it requires the input to be of the form constructed by mcexact
update.cab <- function(object,...){
  cab(object, ...) 
}

cab <- function(args, nosim = NULL, batchsize = NULL, savechain = FALSE, p = NULL, flush = FALSE){
  ##error checking and initializing
  if (!is.null(p)){
    if (!is.double(p)) stop("p must be double")
    else if ((p < 0) | (p > 1)) stop("p must be in [0,1]")
  }
  if (!is.null(batchsize)) args$batchsize <- batchsize
  if (!is.null(nosim)) args$nosim <- nosim
  if (is.null(args$y.start)) y.start <- args$y
  if (is.null(args$startiter) | flush) args$startiter <- 1
  if (is.null(args$phat) | flush) phat <- 0
  else phat <- args$phat
  if (is.null(args$mhap) | flush) args$mhap <- 0
  if (is.null(args$bmsq) | flush) bmsq <- 0
  else bmsq <- args$bmsq
  if (is.null(args$nobatches) | flush) nobatches <- 0
  else nobatches <- args$nobatches
  if (is.null(args$current.batchmean)) current.batchmean <- 0
  else current.batchmean <- args$current.batchmean
  if (savechain){
    nocol <- length(args$stat(rowlabels = TRUE))
    args$chain <- matrix(0, args$nosim, nocol)
    colnames(args$chain) <- args$stat(rowlabels = TRUE)
  }
  else args$chain <- NULL
  
  perpos <- 0
  y.old <- y.start
  y1.old <- y.start[1 : args$n1]  
  for (i in args$startiter : (args$startiter + args$nosim - 1)){    
    shuffle <- sample(1 : args$n1)
    conde1.permute <-  args$conde1[shuffle]
    condv1.permute <-  args$condv1[shuffle, shuffle]
    ##k is the number of elements to be left fixed
    k <- rbinom(1, args$n1 - 1, args$p)
    ##separate y1 into those that stay the same and
    ##those that get updated
    y1.old.permute <- y1.old[shuffle]
    if (k > 0)
      staysfixed <- y1.old.permute[1 : k]
    else
       staysfixed <- NULL
    getsupdated <- y1.old.permute[(k + 1) : args$n1]
    ##multinorm calculates the required mean for
    ##going backwards
    temp <- .Call("multinorm",
                  conde1.permute,
                  condv1.permute,
                  as.double(staysfixed),
                  as.double(y1.old.permute),
                  args$tdf,
                  as.integer(k),
                  PACKAGE="exactLoglinTest")
    ##the ones that got updated
    y1.new.permute <- temp[[1]]
    conde1.old.permute <- temp[[2]]
    changed <- y1.new.permute[(k + 1) : args$n1]

    tf2 <- tf3 <- FALSE

    # HJ
    N.na <- sum(is.na(changed)) 
    if (i==1 & N.na > 0)
	warning("i = ", i, ": 'changed' has ", N.na, " NAs")
    tf1 <- N.na==0 & all(changed >= 0)

    if (tf1){
      y1.new <- y1.new.permute[order(shuffle)]
      y2.new <- round(args$x2invt %*% (args$s - t(args$x1) %*% y1.new))
      tf2 <- all(y2.new >= 0) 
      if (tf2){
        perpos <- perpos + 1
        y.new <- c(y1.new, y2.new)
        target.new <- args$dens(y.new)
        target.old <- args$dens(y.old)
        mean.new <- conde1.permute[(k + 1) : args$n1]
        mean.old <- conde1.old.permute[(k + 1) : args$n1]
        var.old.new <- diag(condv1.permute)[(k + 1) : args$n1]
        cand.new <- rounded.tprob(changed    , mean.new, var.old.new, args$tdf)
        cand.old <- rounded.tprob(getsupdated, mean.old, var.old.new, args$tdf) 
        w <- target.new +  cand.old - target.old - cand.new
        tf3 <- log(runif(1)) <=  w
      }
    }
    if (all(tf1, tf2, tf3)){
      y.old <- y.new
      y1.old <- y1.new
      args$mhap <- args$mhap + 1
    }
    else {
      y.new <- y.old
      ##not necessary but just to remind you
      ##y.old remains y.old for the next iteration
      ##y1.old remains y1.old for the next iteration
    }
    d <- args$stat(y = y.new, mu = args$mu.hat, rowlabels = FALSE)
    if (savechain) args$chain[i - args$startiter + 1,] <- d
    
    phat <- (phat * (i - 1) + (d >= args$dobs)) / i
    ##upddate the batchmean estimate
    batch.iter <- i %% args$batchsize + 1
    if (batch.iter == 1) current.batchmean <- 0
    current.batchmean <- (current.batchmean * (batch.iter - 1) + (d >= args$dobs)) / batch.iter
    if (batch.iter == args$batchsize) {
      bmsq <- bmsq + current.batchmean ^ 2
      nobatches <- nobatches + 1
    }
  }
  args$startiter <- i + 1
  ##keep the current batchmean in case simulation is restarted
  args$current.batchmean <- current.batchmean
  args$bmsq <- bmsq
  args$nobatches <- nobatches
  args$phat <- phat
  args$mcse <- sqrt((bmsq /  nobatches - phat ^ 2) / nobatches)
  args$y1.start <- y1.new
  args$perpos <- perpos / args$nosim
  return(args)
}
