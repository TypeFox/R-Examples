#######################################################################
##
## Function: minimum.entropy()
## Author  : Jonathan Wand <wand@stanford.edu>
##
#######################################################################
minimum.entropy.calc <- function(obj , debug=0 ) {

  cmax <- obj$max
  if (obj$n.interval == 0) { return(NULL) }

  freq.non       <- as.numeric(obj$scalar$N)
  idx <- obj$interval$from != obj$interval$to 
  out.span.count <- as.numeric(obj$interval$N[ idx ])
  tties <- subset( obj$interval, subset=idx , select=c("from","to"))
  rnames <- row.names(tties)
  mties <- as.matrix(tties)
  
  ## how many span cat are there...
  ncombo <- eval(parse(text=paste("(",tties$to, "-", tties$from,"+1)" ,sep="", collapse="+")))

  ridx <- cidx <- nidx <- rep(NA,ncombo)

#  print(ncombo)
#  print(length(out.span.count))
        
  k <- 0
  for (i in 1:length(out.span.count)) {
    #print(tties$from[i]:tties$to[i])
    for (j in tties$from[i]:tties$to[i]) {
      k <- k+1
      cidx[k] <- j
      nidx[k] <- out.span.count[i]
      ridx[k] <- i
    }
  }

  cassign <- rep(NA,length(out.span.count))

  for (i in 1:cmax) {

    if (debug > 1) 
      print( cbind(ridx,cidx,nidx))

    tmp <- tapply( nidx, cidx, sum )
    ct <- as.numeric(names(tmp))

    if (debug > 1) {
      cat("apply\n")
      print(tmp)
    }
    
    ttot <- freq.non
    for (j in ct ) {
      ttot[j] <- freq.non[j] + tmp[ ct == j ]
    }
    ## what C are we going to load on 
    cc <- which.max(ttot)
    ## which idx of expansion of ties does this then pull out
    ee <- cc == cidx
    ## which leads to this index of the row of ties...
    rr <- ridx[ee]
    ## we have a new c assignment for tie now
    cassign[ rr ] <- cc

    if (debug > 1)  {
      cat("test",cc,"\n")
      print(ttot)
      print(ee)
      print(rr)
      cat("cassign\n")
      print(cassign)
    }
    
    if (all(!is.na(cassign))) {
      names(cassign) <- rnames
      return(cassign)
    }
    
    ## now delete used rows from cidx,ridx,nidx and do it again
    cidx <- cidx[ !(ridx %in% rr) ]
    nidx <- nidx[ !(ridx %in% rr) ]
    ridx <- ridx[ !(ridx %in% rr) ] ## do this last...
    ## and from freq.non
    freq.non[cc] <- 0
  }

  
}
