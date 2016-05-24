"plot.dirichlet" <-
function(x,      # object of "dirichlet" class 
                              t=4,      # growth over 2 years (8 quarters)
                              brand=1:x$nbrand,
                              ## Plot Tuning 
                              incr=1,   # plot by increment of 0.1 of a quarter.
                              result=NULL, # Use known result to save time on plot tuning.
                              ...
                              ) {

  tseq <- seq(0,t,incr)
  nt <- length(tseq)
  nb <- length(brand)
  nc <- rainbow(nb)                     # separate color for each brand

  if (is.null(result)) {
    ## store results for Pen and Buy Rate
    r.pen <- matrix(0,nrow=nb,ncol=nt) # 1st column is time=0 (pen=0)
    r.buy <- matrix(1,nrow=nb,ncol=nt)

    ## speed up with function 
    bmat <- function(j){
      ct <- tseq[j]                        # Current Time Period
      x$period.set(ct)                   # change the time period to "ct"
      
      ## Observed Market Share (Per Capita Purchase Rate)
      ## Theoreticl Brand Penetration is:                                        
      r.pen[,j] <<- sapply(brand, x$brand.pen)
      ## Purchase Rate of the Brand
      r.buy[,j] <<- sapply(brand, x$brand.buyrate) 
    }
    
    sapply(2:nt, bmat)                       # j=1 is "time period =0"
  }
  else {
    r.pen <- result$pen
    r.buy <- result$buy
  }

  ## Penetration Plot
  par(mfrow=c(1,2))
  
  plot(c(0,rep(t,nb)), c(0,r.pen[,nt]), type="n",
       main=paste("Theoretical Penetration Growth of Retailer Over",t,"Quarters"),
       xlab="Quarters",ylab="Penetration")
  abline(h=seq(0.1,1,0.05),col="gray90")
  for (i in 1:nb) {
    lines(tseq, r.pen[i,], lty=i,lwd=2,col=nc[i])
  }
  legend("topleft",x$brand.name,lty=1:nb,lwd=3,col=nc)

  ## Buy Rate Plot
  plot(matrix(tseq,nrow=nb,ncol=nt,byrow=T), r.buy, type="n",
       ylim=c(1,max(r.buy)+0.1),
       main=paste("Theoretical Shopping Rate Growth Over", t, "Quarters"),
       xlab="Quarters",ylab="Shopping Frequency")
  abline(h=seq(2,30,0.5),col="gray90")        # make horizong lines 
  for (i in 1:nb) {
    lines(tseq, r.buy[i,], lty=i,lwd=2,col=nc[i])
  }
  legend("topleft",x$brand.name,lty=1:nb,lwd=3,col=nc)

  par(mfrow=c(1,1))
  
  rownames(r.pen) <- x$brand.name
  rownames(r.buy) <- x$brand.name
  colnames(r.pen) <- tseq
  colnames(r.buy) <- tseq
  
  list(pen=r.pen, buy=r.buy)
}







