lpc <- function(X,h, t0=mean(h),  x0=1,   way = "two", scaled=TRUE,  weights=1, pen=2,
              depth=1, control=lpc.control()){
  iter    <-  control$iter
  boundary <- control$boundary
  convergence.at <- control$convergence.at
  thresh  <- control$pruning.thresh
  rho0    <- control$rho0
  cross   <- control$cross
  mult    <- control$mult

  Xi      <- as.matrix(X)
  N       <- dim(Xi)[1]
  d       <- dim(Xi)[2]
  if (N %% length(weights) !=0){
       stop("The length of the vector of weights is not a multiple of the sample size.")
  } else {
      weights <-rep(weights, N %/%   length(weights))
  }     
 
  countb   <- 0                # counter for branches 17/07/2007
  Lambda   <- matrix(0,0,5)    # stores the parameters (lambda, countb) as well as  depth, order, and  initialization number of the corresponding branch.
  s1       <- apply(Xi, 2, function(dat){ diff(range(dat))})  # range
  if (missing(h)){ if (!scaled){h   <- s1/10} else {h<- 0.1}} # bandwidth by default: 10% of range
  if (length(h)==1){h <- rep(h,d)}
  ms.h    <- if (!is.null(control$ms.h)){control$ms.h} else {h}
  ms.sub  <- control$ms.sub 
  if (scaled){        # scales the data to lie by its range
        Xi <- sweep(Xi, 2, s1, "/")
   } 
  if (d==1){
      stop("Data set needs to consist of at least two variables!")
  } 

  # Selecting and scaling starting points:
  if (length(x0)==1 && x0==1){
    
      n  <- sample(N,1);
      x0 <- matrix(ms.rep(Xi, Xi[n,],ms.h, plotms=0)$final, nrow=1)
      
   }  else if (length(x0)==1 && x0==0 ){
          if (N<= ms.sub) {
              sub<- 1:ms.sub
          }  else {    
              Nsub <- min(max(ms.sub, floor(ms.sub*N/100)), 10*ms.sub)
              #print(Nsub)
              sub <- sample(1:N, Nsub)
          }
          x0<- ms(Xi, ms.h, subset=sub, plotms=0)$cluster.center
   }  else {
        if (is.null(x0)){
             if (is.null(mult)){stop("One needs to allow for at least one starting point.")}
             x0<- matrix(0, nrow=0, ncol=d)
        } else {    
             x0 <- matrix(x0, ncol=d, byrow=TRUE)
        }   
      if (scaled){ x0 <- sweep(x0, 2, s1, "/") } # scales the starting point just as the scaled data
  }

 
  if(!is.null(mult)){  
      if (dim(x0)[1] < mult) {n <- runif(mult-dim(x0)[1],1,N+1)%/%1; x0 <- rbind(x0,Xi[n,])} #
      if (dim(x0)[1] > mult) {x0<-x0[1:mult]} 
  }
  x0       <- as.matrix(x0)    # putting in matrix format; just in case.... 
  mult     <- dim(x0)[1]       # final number of starting points used
  X1       <- matrix(0,mult,d) # collects starting points for branches of depth 1.
  X2<-X3   <- matrix(0,0,d)    # collect starting points for branches of depth 2 and 3.
  S2<-S3   <- matrix(0,0,d)    # collects candidates for starting points of depth 2 and 3.
  saveall  <- matrix(0,0,d)    # stores the LPC

  
  gapsize <- 1.5 # Note: In the original paper, this was set to 2, in order to avoid higher-depth branches to run back into the main branch. 
                 # This  effect is now partly achieved through "jweights", see below.
  
  for (j in 1:mult){
    
     xo           <- x0[j,]                       # defines appropriate starting point for the jth curve
     X1[j,]       <- xo                           # adds starting point to list  
     curve0       <-  followx(Xi, xo, h, t0, iter, way, weights, pen, phi =1, 0,rho0,  boundary, convergence.at,  cross ) # computes LPC of depth 1
     saveall      <- rbind(saveall,curve0[[1]])	  # stores LPC
     l            <- dim(curve0[[5]])[1]	  # number of candidates for junctions
     Lambda       <- rbind(Lambda, cbind(curve0[[6]], countb,1,1,j))
     countb       <- countb + 1
     dimnames(Lambda)[[2]]<- c("lambda", "branch", "depth", "order", "init")

     
     if (depth >1 && l>=1 ){							  # constructs branches of depth 2.

       for (s in 1:l){

          x          <- curve0[[5]][s,]				        # candidates for junctions on initial curve.
          center.x   <- cov.wt(Xi, wt= kernd(Xi,x,h)*weights)           # computes local covariance and mean at x
          eigen.cov  <- eigen(center.x[[1]])  				# eigenvalues and eigenvec's of local cov matrix
          eigen.vecd <- eigen.cov[[2]][,2]				# second local eigenvector
          c.x        <- center.x[[2]]                                   # local center of mass around x
          new.x      <- c.x+ gapsize*t0 * eigen.vecd 		# double stepsize to escape from initial curve 
          #print( kdex(Xi,new.x, h))
          if (kdex(Xi,new.x, h) >= thresh){ 				# Pruning
             X2     <- rbind(X2, t(new.x))				# adds new starting point to list          
             x      <- new.x
             jweights<-  1-kernd(Xi, c.x, h)/ kernd(c.x,c.x,h)
             curve1 <- followx(Xi, x, h, t0, iter, way ="one", weights*jweights, pen,  phi=2, lasteigenvector= eigen.vecd, rho0, boundary,convergence.at,  cross)
             saveall<- rbind(saveall, curve1[[1]])    # stores LPC
             S2     <- rbind(S2,curve1[[5]])          # stores candidates for further junctions
             Lambda <- rbind(Lambda, cbind(curve1[[6]], countb, 2, 2,j))
             countb <- countb + 1
          }
          new.x <- c.x- gapsize*t0 * eigen.vecd 		# go in opposite direction
           # print( kdex(Xi,new.x, h))      
          if (kdex(Xi,new.x, h) >= thresh){				# Pruning
             X2     <- rbind(X2, t(new.x))
             x      <- new.x
             jweights<-  1-kernd(Xi, c.x, h)/ kernd(c.x,c.x,h)
             curve1 <- followx(Xi, x, h, t0, iter, way="back", weights*jweights, pen, phi=2, lasteigenvector=-eigen.vecd, rho0,boundary, convergence.at,  cross)
             saveall<- rbind(saveall, curve1[[1]])
             S2     <- rbind(S2, curve1[[5]])
             Lambda <- rbind(Lambda, cbind(curve1[[6]], countb, 2,2,j))
             countb <- countb + 1     
          } # end if kdex.....
      } # end for (s)
     }# end if (depth)

     k <- dim(S2)[1]							      # no. of candidates for starting points of depth 3.
     if (depth >2 && k>=1 ){					# constructs branches of depth 3
         for (s in 1:k){

            x         <- S2[s,]			  # candidates for junctions on curve of depth 2.
            center.x  <- cov.wt(Xi, wt= kernd(Xi,x,h)*weights)
            eigen.cov <- eigen(center.x[[1]])
            eigen.vecd<- eigen.cov[[2]][,2]
            c.x        <- center.x[[2]]   
            new.x     <- c.x+ gapsize*t0 * eigen.vecd
             # print( kdex(Xi,new.x, h))
            if (kdex(Xi,new.x, h)  >= thresh){
                  X3      <- rbind(X3, t(new.x))
                  x       <- new.x
                  jweights<-  1-kernd(Xi, c.x, h)/ kernd(c.x,c.x,h)
                  curve1  <- followx(Xi, x , h, t0, iter, way ="one", weights*jweights, pen, phi=3, lasteigenvector= eigen.vecd, rho0, boundary, convergence.at, cross)
                  saveall <- rbind(saveall, curve1[[1]])
                  S3      <- rbind(S3,curve1[[5]])
                  Lambda <- rbind(Lambda, cbind(curve1[[6]], countb, 3,2,j))
                  countb <- countb + 1 
                }
            new.x <- c.x- gapsize*t0* eigen.vecd
             # print( kdex(Xi,new.x, h) )       
            if (kdex(Xi,new.x, h)   >= thresh){
                  X3     <- rbind(X3, t(new.x))
                  x      <- new.x
                  jweights<-  1-kernd(Xi, c.x, h)/ kernd(c.x,c.x,h)
                  curve1 <- followx(Xi, x, h, t0, iter, way="back", weights*jweights, pen, phi=3, lasteigenvector=-eigen.vecd, rho0, boundary, convergence.at, cross)
                  saveall<- rbind(saveall, curve1[[1]])
                  S3     <- rbind(S3,curve1[[5]])
                  Lambda <- rbind(Lambda, cbind(curve1[[6]], countb, 3,2,j))
                  countb <- countb + 1
            } # end if kdex....
         } # end for (s)
     } # end if (depth)
  } # end for (j)

#print(X1)
  
   # if (d==2){
      starting.points <- rbind(X1,X2,X3) 
      dimnames(starting.points)<-list(c(rep(1, dim(X1)[1]), rep(2, dim(X2)[1]), rep(3, dim(X3)[1]) ), NULL)
   # }
   # else {    # removed 05/05/11
   #   starting.points <- as.matrix(X1)
   #   dimnames(starting.points)<-list(c(rep(1, dim(X1)[1])), NULL)
   #}


  
  fit<-  c(LPC=list(saveall),
           Parametrization= list(Lambda),
           h=list(h),
           t0=t0,
           starting.points=list(starting.points),
           data=list(Xi),
           scaled=list(scaled),
           weights=list(weights),
           control=list(control),
           Misc=list(list( rho=curve0[[4]], scaled.by = if (scaled) s1 else rep(1,d), adaptive.band=curve0[[7]] , angles=curve0[[3]]))
           ) 
   class(fit)<-"lpc"
   fit
  
} # end lpc function
