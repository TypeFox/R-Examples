simtable.cab <- function(args, nosim = NULL, p = NULL, y.start = NULL){
  ##error checking and initializing
  if (!is.null(p)){
    if (!is.double(p)) stop("p must be double")
    else if ((p < 0) | (p > 1)) stop("p must be in [0,1]")
  }
  if (!is.null(nosim)) args$nosim <- nosim
  if (is.null(y.start)) y.start <- args$y
  else if (t(args$x) %*% y.start != args$s)
    stop("invalid starting value")
           
  nocol <- length(args$y)
  chain <- matrix(0, args$nosim, nocol)

  y.old <- y.start
  y1.old <- y.start[1 : args$n1]  
  for (i in 1 : args$nosim){    
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
    chain[i,] <- y.new
  }
  return(chain[,order(args$ord)])
}
