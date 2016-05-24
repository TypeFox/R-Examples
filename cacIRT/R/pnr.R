pnr <- function(resp, bw.g = NULL, alpha = .5){  
  #pnr = poly.non.ram

  N<-dim(resp)[1]
  I<-dim(resp)[2]
  opt <- as.numeric(levels(factor(resp)))
  M <- length(opt) 
  
  tots <- rowSums(resp, na.rm=T)
  high <- max(tots, na.rm=T)	
  low  <- min(tots, na.rm=T)
  x.points <- low:high
  n.x <- length(x.points)
  
  density <- hist(tots, breaks = seq(low-.5,(high+.5), by=1), plot=FALSE)$density
  
  kden <- density(tots, bw = "nrd", n = n.x, from = low, to = high)
  
  if(is.null(bw.g)){ 
    bw.g <- bw.nrd(tots)
  }
  
  #adaptive bandwidth
  g.mean <- exp(mean(log(kden$y))) 
  lambda <- (kden$y/g.mean)^-alpha		
  
  ks<-array(NA, dim = c(n.x, M, I ) )
  for(j in 1:I) {
    
    ts.j <- tots

    for(m in 1:M){
      an <-opt[m]
      temp <- ifelse(resp[,j]==an, temp<-1, temp<-0)			
      for(i in 1:n.x){
        
        ks[ i, m, j] <- sum( 
          exp(-((ts.j - x.points[i])/(bw.g*lambda[i]))^2/2)
          *temp) /
          sum( 
            exp(-((ts.j - x.points[i])/(bw.g*lambda[i]))^2/2)
            )	
      }}}
  
  return(list(x.points, density, ks))
  }

