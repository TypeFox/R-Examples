  
  
 ##################################
 #### inner-products matrix Gw ####
 ##################################

 ## Description:  Compute the nxn inner-products matrix Gw.
 ##               If G is positive semidefinite the squared distances 
 ##               Delta become Euclidean
 ##
 
 Gcalc<-function(n,weights,Delta){  #trure la n
 
   In <- diag(nrow=n) # identity matrix
   onesn <- matrix(rep(1,n),nrow=n) #ones vector
   weights <- as.matrix(weights)
   
   dw <-Delta%*%weights  # weighted average for each row of Delta
   vw <-t(weights)%*%dw  # scalar constant

   # G matrix of inner-products calculated as the difference between Delta and 
   # the averages for rows and columns of Delta plus the sum of the constant vw.
   G <-(-0.5)*(Delta- onesn%*%t(dw)-dw%*%t(onesn)+onesn%*%t(onesn)*vw[1]) 
   return(G)
}                   