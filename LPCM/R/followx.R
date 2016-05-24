
followx<-function(Xi, x0, h, t0, iter, way,  weights, pen=2,  phi =1, lasteigenvector=0, rho0=0.4, boundary=0.005, convergence.at= 0.000001, cross=TRUE){

   N  <- dim(Xi)[1]
   d  <- dim(Xi)[2]
   
   if (way=="one"){ b<-1}else if (way %in% c("two","back")) {b<-2}
   else {stop("way has to be one of: one, two, back.")}

   highrhopoints <- matrix(0,0,d)					# here points with high local 2nd eigenvalue will be stored
   save.xd      <- eigen.vecd<-matrix(0, b*iter,d)			# here the LPC will be stored
   cos.alt.neu  <- cos.neu.neu  <- rep(1, b*iter)
   lambda       <- rep(NA, b*iter)                                     # Parametrization
   rho          <- rep(0, b*iter)					# construction of empty vectors
   c0           <- rep(1,b*iter)
   x <- x0

   
   if (way == "one" || way =="two"){

     for (i in iter:1){
      
       center.x <- cov.wt(Xi, wt= kernd(Xi,x,c0[i]*h)*weights)           # local covariance and mean at x
       mu.x     <- center.x[[2]] 				         # Local center of mass at x (mean shift).
       save.xd[i,]  <- mu.x						  # stores the first point of that branch.

       if (i==iter){
          lambda[i] <- 0
       } else {
          lambda[i] <-  lambda[i+1] + sqrt(sum( (mu.x-save.xd[i+1,])^2))
       }# measures distance for parameterization


       eigen.cov              <- eigen(center.x[[1]]) 				  # local covariance at x
       eigen.vecd[i,]         <- eigen.cov[[2]][,1]				  # local first eigenvector at x

       rho[i]<-eigen.cov[[1]][2]/eigen.cov[[1]][1]                               # computes rho
       if ((i< iter) && (rho[i] > rho0) && (rho[i+1]<= rho0)){#print(i); print(x)
                         highrhopoints<-rbind(highrhopoints, x)
       }                                                                         # compares with rho0

       if (i == iter && lasteigenvector[1] !=0){ cos.alt.neu[i]<- sum(lasteigenvector* eigen.vecd[i,]) }
       if (i < iter){ cos.alt.neu[i]<- sum(eigen.vecd[i+1,]* eigen.vecd[i,])}
									   # angle between current and previous eigenvector
       if (cos.alt.neu[i]<0){ eigen.vecd[i,]<- - eigen.vecd[i,]
                              cos.neu.neu[i] <- - cos.alt.neu[i]
       } else {
                    cos.neu.neu[i] <- cos.alt.neu[i]
       }
        							         # signum flipping
       
       if (pen>0 && i<iter){
             a<-(abs(cos.alt.neu[i]))^pen
             eigen.vecd[i,]<-a* eigen.vecd[i,] + (1-a)*eigen.vecd[i+1,]
        }									 # angle penalization

       if (pen>0 && i == iter && lasteigenvector[1] !=0 ){
             a<-(abs(cos.alt.neu[i]))^pen
             eigen.vecd[i,]<-a*  eigen.vecd[i,] + (1-a) * lasteigenvector
        }									 # angle penalization


     new.x    <- center.x[[2]] + eigen.vecd[i,]*t0			 # Update
     x        <- new.x
      if (i!=iter && i!=1){
         if (!cross){
             prox <- which(sqrt(d)*distancevector(as.matrix(save.xd[iter:i,]), mu.x)<= mean(h))       
             if (length(prox)!=(diff(range(prox))+1)){break}
         } 
        
         if ( abs(lambda[i]-lambda[i+1])/abs(lambda[i]+lambda[i+1]) <convergence.at){break}
         if ( abs(lambda[i]-lambda[i+1])/abs(lambda[i]+lambda[i+1]) <boundary   ){c0[i-1]<-0.995*c0[i]}
         #else {h<-ifelse(1.01*h<h0,1.01*h,h0)}
         else {c0[i-1]<-min(1.01*c0[i],1)}
    }
      #print(c0[i]*h)
     }
 }

 if (way=="back" || way =="two" ){					 # opposite direction

#    h <- h0
    x <- x0
  
    for (i in (iter+1):(2*iter)){
     
      center.x<-cov.wt(Xi, wt= kernd(Xi,x,c0[i]*h)*weights)
      mu.x <- center.x[[2]]
      save.xd[i,]            <- mu.x
      if (i==iter+1){
          lambda[i] <- -sqrt(sum( (mu.x-save.xd[iter,])^2))
      }  else {
          lambda[i] <- lambda[i-1] - sqrt(sum( (mu.x-save.xd[i-1,])^2))
      }                                  # measures distance for parametrization
      eigen.cov              <- eigen(center.x[[1]]) 				 # local eigenvalues and eigenvectors
      eigen.vecd[i,]         <- eigen.cov[[2]][,1]

      rho[i]<-eigen.cov[[1]][2]/eigen.cov[[1]][1]
      if (i > iter +1 && rho[i] > rho0 && rho[i-1]<= rho0){#print(i); print(x)
              highrhopoints<-rbind(highrhopoints,x)
      }
      if (i == (1 + iter) && lasteigenvector [1] != 0){ cos.alt.neu[i]<- -sum(lasteigenvector* eigen.vecd[i,])}
      if (i >= (2 + iter)){ cos.alt.neu[i]<- sum(eigen.vecd[i-1,]* eigen.vecd[i,])}
      if (cos.alt.neu[i]<0){
                       eigen.vecd[i,]<- - eigen.vecd[i,]
                       cos.neu.neu[i] <- - cos.alt.neu[i]
      } else {
                    cos.neu.neu[i] <- cos.alt.neu[i]
      }

        

      if (pen>0 && i>=(2+iter)){
              a<-(abs(cos.alt.neu[i]))^pen
              eigen.vecd[i,]<- a* eigen.vecd[i,]+(1-a)* eigen.vecd[i-1,]
      }

      if (pen>0 && i == 2 +iter && lasteigenvector[1] !=0 ){
            a<-(abs(cos.alt.neu[i]))^pen
            eigen.vecd[i,]<-a * eigen.vecd[i,] + (1-a)*lasteigenvector
      }

      new.x    <-  center.x[[2]] - eigen.vecd[i,]*t0
      #lambda[i]<-  -sqrt(sum( (x0-new.x)^2))
      x        <-  new.x
      if ( i!= iter &&  i!=b*iter &&  !is.na(lambda[i-1])   && abs(min(lambda[i],lambda[i-1]))!=0){
      if (!cross){
            prox <- which(sqrt(d)*distancevector(as.matrix(save.xd[1:i,]), mu.x)<= mean(h))
            if (length(prox)!=(diff(range(prox))+1)){break}
         }   
        
         if (abs(lambda[i]-lambda[i-1])/abs(lambda[i]+lambda[i-1]) <convergence.at){break}
         if (abs(lambda[i]-lambda[i-1])/abs(lambda[i]+lambda[i-1]) <boundary  ){c0[i+1]<-0.995*c0[i]}
         #else {h<-ifelse(1.01*h<h0,1.01*h,h0)}
         else {c0[i+1]<-min(1.01*c0[i],1)}
      }
      #print(c0[i]*h) 
  } # for (i)


 # if (way=="back"){save.xd<-save.xd[(iter+1):(2*iter),]  ;
 #                  eigen.vecd<- eigen.vecd[(iter+1):(2*iter),]
 #                  }  # removed 17/07/2007
 }# if (way)

   #print(lambda)
   filter <- !is.na(lambda); if (way =="two" ){filter[iter+1]<- FALSE}###

   list(save.xd[filter,],eigen.vecd[filter,],cos.neu.neu[filter], rho[filter], highrhopoints, round(-lambda[filter]+max(lambda[filter]), digits=4), c0)
  
}# function
