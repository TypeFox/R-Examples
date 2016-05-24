############################# Utils

.mergegaps <- function(gaps, support){
 if(is.null(gaps)) return(NULL)
 if(is.null(support)) return(gaps)
 if(length(support)==0) return(gaps)
 
 mm <- rbind(cbind(gaps[ ,1],1),
             cbind(gaps[ ,2],2),
             cbind(support,3))
 mm <- mm[order(mm[,1]),]
 jj <- 0
 ein <- FALSE
 gaps.new <- mm
 for(j in 1:nrow(mm))
  {if (mm[j,2] == 1 & !ein){
       jj <- jj + 1
       ein <- TRUE
       gaps.new[jj,1] <- mm[j,1]
       }
   if (mm[j,2] == 3 &  ein){
       gaps.new[jj,2] <- mm[j,1]
       jj <- jj + 1
       gaps.new[jj,1] <- mm[j,1]
       }
   if (mm[j,2] == 2)       {
       ein <- FALSE
       gaps.new[jj,2] <- mm[j,1]
       }
   }
 ga <- gaps.new[1:jj,,drop = FALSE]
 if(jj>0){
    ln <- seq(jj)
    ga <- matrix(ga[ln[ga[,1] != ga[,2]], ], ncol = 2)
    }
 else ga <- NULL
 return(ga)

}

.mergegaps2 <- function(gaps1, gaps2){
 if(is.null(gaps1)) return(NULL)
 if(is.null(gaps2)) return(NULL)
 mm <- rbind(cbind(gaps1[ ,1],1), 
             cbind(gaps1[ ,2],2),
             cbind(gaps2[ ,1],3),
             cbind(gaps2[ ,2],4))
 mm <- mm[order(mm[ ,1]), ]
 jj <- 0
 state1 <- 0
 state2 <- 0
 gaps.new <- mm
 for(j in 1:nrow(mm))
  {if (mm[j,2] == 1 && state2 == 0){
        state1 <- 1
       }
   if (mm[j,2] == 1 && state2 == 1){
        jj <- jj + 1
        gaps.new[jj,1] <- mm[j,1]
        state1 <- 1
       }
   if (mm[j,2] == 2 && state2 == 1){
        gaps.new[jj,2] <- mm[j,1]
        state1 <- 0
       }
   if (mm[j,2] == 2 && state2 == 0){
        state1 <- 0
       }
   if (mm[j,2] == 3 && state1 == 0){
        state2 <- 1
       }
   if (mm[j,2] == 3 && state1 == 1){
        jj <- jj + 1
        gaps.new[jj,1] <- mm[j,1]
        state2 <- 1
       }
   if (mm[j,2] == 4 && state1 == 1){
        gaps.new[jj,2] <- mm[j,1]
        state2 <- 0
       }
   if (mm[j,2] == 4 && state1 == 0){
        state2 <- 0
       }
   }
 erg <- if (jj > 0) gaps.new[1:jj, ,drop=FALSE] else NULL
 return(.consolidategaps(erg))

}

.consolidategaps <- function(gaps){
 if(is.null(gaps)) return(NULL)
 if(nrow(gaps)==0) return(NULL)
 if(nrow(gaps)==1) return(gaps)
 jj <- 1
 for(j in 1:(nrow(gaps)-1)){
   if (.isEqual(gaps[jj,2],gaps[j+1,1]))
      gaps[jj,2] <- gaps[j+1,2]
   else jj <- jj+1   
 }     
 return(gaps)
}

.pmixfun <- function(mixDistr, mixCoeff, leftright = "right"){
  l <- length(mixCoeff)
  return(function(q, lower.tail = TRUE, log.p = FALSE ){
  p0 <- as.vector(
           matrix(unlist(
              lapply(mixDistr, 
                  function(x){
                       p.lr <- if(pmatch(leftright, table=c("left","right"), 
                                         nomatch = 2)==2) p.l(x) 
                               else x@p
                       do.call(p.lr, list(q = q, lower.tail = lower.tail)) 
                  } 
              )),
              ncol = l, nrow = length(q)) %*% mixCoeff
           )
  if(log.p) p0 <- log(p0)
  return(p0)               
   })
}
.dmixfun <- function(mixDistr, mixCoeff, withStand = FALSE, supp = NULL){
  l <- length(mixCoeff)
  if(withStand) {su <- sum(as.vector(matrix(unlist(lapply(mixDistr, function(y)
                                               do.call(y@d, list(x = supp)))),
                                 ncol = l, nrow = length(supp)) %*% mixCoeff))
                }else su <- 1
 
  return(function(x, log = FALSE ){
         d0 <- 
         as.vector(matrix(unlist(lapply(mixDistr, function(y)
            do.call(y@d, list(x = x)))),
         ncol = l, nrow = length(x)) %*% mixCoeff)/su
         if(log) d0 <- log(d0)
         return(d0)
         })
}

.rmixfun <- function(mixDistr, mixCoeff){
  l <- length(mixCoeff)
  return(function(n){
  Un <- outer(sample(1:l, size = n, replace = TRUE, prob = mixCoeff),
                          1:l, function(x,y) x == y)
  mal <- lapply(mixDistr, function(x) x@r(n))
  ma <- matrix(unlist(mal), ncol = l, nrow = n)
  return(rowSums(Un*ma))})

}


.qmixfun <- function(mixDistr, mixCoeff, Cont = TRUE, pnew, gaps = NULL, 
                     leftright = "left"){
  l <- length(mixCoeff)
  if(l==0) return(NULL)
  loup <- .loupmixfun(mixDistr)

  n <- getdistrOption("DefaultNrGridPoints")
  up <- if(is.finite(loup$qu)) loup$qu else  1000 #min(loup$qu,1000)#getdistrOption("LARGE")
  lo <- if(is.finite(loup$ql)) loup$ql else -1000 #max(loup$ql,-1000)#getdistrOption("LARGE")

  h <- (up-lo)/n
  suppsA <- NULL
  for (i in 1:l){
    if(!is(try(su0 <- support(mixDistr[[i]]), silent=TRUE),"try-error"))
       suppsA <- c(suppsA,su0)
    }

  xseq <- c(seq(from = lo, to = up, by = h),
            suppsA,
            suppsA-getdistrOption("DistrResolution"))
  xseq <- sort(unique(xseq))          
  if(length(xseq)<2)
     xseq <- c(min(lo,up)-0.1,max(lo,up)+0.1)

  px.l <- pnew(xseq, lower.tail = TRUE)
  px.u <- pnew(xseq, lower.tail = FALSE)
  qnew <- .makeQNew(xseq, px.l, px.u, TRUE, lo, up, Cont = Cont)
  if(!is.null(gaps)) 
      qnew <- .modifyqgaps(pfun = pnew, qfun = qnew, gaps = gaps, 
                           leftright = leftright)
  return(qnew)
}


.loupmixfun <- function(mixDistr){
    if(length(mixDistr)==0) return(list(qL = NA, ql = NA, qU = NA, qu = NA))
    if(length(mixDistr)==1){
      q1 <- q(mixDistr[[1]])
      return(list(qL = q1(p = 0, lower.tail = TRUE),
                  ql = q1(p = getdistrOption("TruncQuantile"), lower.tail = TRUE),
                  qU = q1(p = 0, lower.tail = FALSE),
                  qu = q1(p = getdistrOption("TruncQuantile"), lower.tail =FALSE)
                  ))
    }
    qL0 <- as.vector(unlist(lapply(mixDistr, function(x)
                         do.call(x@q,list(p = 0, lower.tail = TRUE)))))
    qL1 <- as.vector(unlist(lapply(mixDistr, function(x)
                         do.call(x@q,list(p = getdistrOption("TruncQuantile"),
                                 lower.tail = TRUE)))))
    qL  <- min(qL0); ql <- min(qL1)

    qU0 <- as.vector(unlist(lapply(mixDistr, function(x)
                         do.call(x@q,list(p = 0, lower.tail = FALSE)))))
    qU1 <- as.vector(unlist(lapply(mixDistr, function(x)
                         do.call(x@q,list(p = getdistrOption("TruncQuantile"),
                                 lower.tail = FALSE)))))
    qU  <- max(qU0); qu <- max(qU1)
    return(list(qL = qL, ql = ql, qU = qU, qu = qu))
    }

.del0dmixfun <- function(mixDistr){
  dac <- mixDistr@mixDistr[[1]]@d
  if(!is.null(dac)){
      dnew <- function(x, log = FALSE, ...){
               d0 <- dac(x, log = log, ...)
               d0[.isEqual(x,0)] <-  0
               return(d0)
            }
      mixDistr@mixDistr[[1]]@d <- dnew
  }
  return(mixDistr)
}

.ULC.cast <- function(x){
         if( is(x,"AbscontDistribution"))
             x <- as(as(x,"AbscontDistribution"), "UnivarLebDecDistribution")
         if(is(x,"DiscreteDistribution"))
             x <- as(as(x,"DiscreteDistribution"), "UnivarLebDecDistribution")
         if(!is(x,"UnivarLebDecDistribution"))
            x <- as(x,"UnivarLebDecDistribution")
         return(x)
}
