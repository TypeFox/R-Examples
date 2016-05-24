# Author: D.R. Bickel
# Modications: P. Poncet
hsm <-
function(x,                  # sample (the data)
         bw = NULL,          # bandwidth (fraction of the observations to consider)
         k,                  # length of the intervals
         tie.action = "mean",
         tie.limit = 0.05,
         ...)
{
##################################################
# Robertson and Cryer's / FSM / HSM mode estimator
# FSM = fraction-of-sample mode
# HSM = half-sample mode
##################################################

  if (!missing(k) & is.null(bw)) {
    bw <- (k+1)/length(x)
  } else if (missing(k) & is.null(bw)) {
    bw <- 1/2
  }

  if (is.numeric(bw)) {
    if (bw <= 0 | bw > 1) stop("argument 'bw' must belong to (0, 1]")
  }
  
  y <- sort(x)
   
  while (length(y) >= 4) {
    ny <- length(y)
    if (is.function(bw)) {
      k <- ceiling(bw(ny, ...)*ny) - 1
    } else {
      k <- ceiling(bw*ny) - 1
    }

    inf <- y[1:(ny-k)]
    sup <- y[(k+1):ny]
    diffs <- sup - inf
    i <- which(diffs==min(diffs))
    
    ## Ties?
    if(length(i) > 1) i <- .deal.ties(ny, i, tie.action, tie.limit) 
   
    if (diffs[i]==0) {
      y <- y[i]
    } else {
      y <- y[i:(i+k)]
    }
    #y <- ifelse(diffs[i]==0, y[i], y[i:(i+k)])
  }
  if (length(y) == 3) {
    z <- 2*y[2] - y[1] - y[3]
    M <- switch(as.character(sign(z)), "-1" =  mean(y[1:2]), "1" = mean(y[2:3]), "0" = y[2])
  } else {
    M <- mean(y)
  }

  ## Output
  return(M)
}

#HSM <- hsm


#! Recursive estimator
#mlv.hsm.rec <-
#function(x,  # sample
#         bw, # fraction of the observations to consider
#         tie.action = "mean",
#         tie.limit = 0.05)
#{
###############################
# Bickel's / HSM mode estimator
# HSM = half-sample mode
###############################

#  if (missing(bw)) bw <- 1/2
#  if (bw <= 0 || bw > 1) stop("Argument 'bw' must belong to (0, 1].")

#  aux <-
#  function(y, ny, k)
#  {
#    if (ny == 3) {
#      z <- 2*y[2] - y[1] - y[3]
#      return(switch(as.character(sign(z)), "-1" =  mean(y[1:2]), "1" = mean(y[2:3]), "0" = y[2]))
#    } else {
#      if (ny < 3) {
#        return(mean(y))
#      } else {
#        diffs <- y[(k+1):ny] - y[1:(ny-k)]
#        i <- which(diffs==min(diffs))
#        if(length(i) > 1) i <- .deal.ties(ny, i, tie.action, tie.limit)        
#        y <- ifelse(diffs[i]==0, y[i], y[i:(i+k)])
#        ny <- length(y)
#        Recall(y, ny, ceiling(ny*bw)-1)
      
#      }
#    }
#  }
  
  ## Output
#  nx <- length(x)
#  return(aux(sort(x), nx, ceiling(bw*nx)-1))
   
#}
