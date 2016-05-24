
# This function computes the projections of new data points onto the curve, VERY SLOW!

lpc.project.spline <- function(lpcsl, newdata, num.knots=100, optimize=TRUE) {
    squared.distance <- function(t,coords,sps) {
      dist <- 0
      for (i in 1:length(sps))
        dist <- dist + (coords[i]-sps[[i]](t))^2
      dist
    }
    nbranch <- length(lpcsl)-1      

    fit.spline <- lpc.fit.spline(lpcsl, num.knots=num.knots)
    knots.coords  <- fit.spline[[1]]
    knots.pi      <- fit.spline[[2]] 
    newdata <- as.matrix(newdata)
   
    closest.branch <- numeric(nrow(newdata))
    closest.or.pi <- numeric(nrow(newdata))
    closest.idx <- numeric(nrow(newdata))
    closest.dist <- numeric(nrow(newdata))
 
     for (i in 1:nrow(newdata)) {
      sqdist <- Inf
      for (r in 0:nbranch) {
       
        cur.sqdist <- apply((knots.coords[[r+1]]-newdata[i,])^2,2,sum)
        if (min(cur.sqdist)<sqdist) {
          closest.branch[i] <- r
          closest.idx[i] <- which.min(cur.sqdist)
          closest.or.pi[i] <- knots.pi[[r+1]][closest.idx[i]]
          sqdist <- min(cur.sqdist)
         }
      }
      if (optimize) {
        from <- knots.pi[[closest.branch[i]+1]][max(1,closest.idx[i]-1)]
        to <- knots.pi[[closest.branch[i]+1]][min(num.knots,closest.idx[i]+1)]
        if (from==to){
            closest.dist[i] <- sqrt(squared.distance(closest.or.pi[i],coords=newdata[i,],sps=lpcsl[[closest.branch[i]+1]]$splinefun))
        }  else {           
           opt <- optimize(squared.distance,lower=from,upper=to,coords=newdata[i,],sps=lpcsl[[closest.branch[i]+1]]$splinefun)
           closest.or.pi[i] <- opt$minimum
           closest.dist[i] <- sqrt(opt$objective)
        }
     } else {
        closest.dist[i] <- sqrt(squared.distance(closest.or.pi[i],coords=newdata[i,],sps=lpcsl[[closest.branch[i]+1]]$splinefun))
     }
    }
    closest.pi <- lpc.curve.length(lpcsl,closest.or.pi,closest.branch)
    closest.coords <- lpc.spline.eval(lpcsl, closest.or.pi, closest.branch)

    result <- list(closest.branch=closest.branch,closest.pi=closest.pi,closest.or.pi=closest.or.pi,closest.coords=closest.coords,closest.dist=closest.dist)
    result
  }

