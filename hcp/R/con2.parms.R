con2.parms <-
function(x,y,n,j0,k0,e10,e20){

th <- matrix(0,100,7)

# Iteration 0

 th[1,1] <- e10
 th[1,2] <- e20
 bc <- beta2.calc(x,y,n,j0,k0,e10,e20)
 th[1,3:7] <- bc$B

# Iterate to convergence (100 Iter max)

 for (iter in 2:100){
 m <- iter-1
 ec <- eta2.calc(x,y,n,j0,k0,th[m,3:7])
 th[iter,1] <- ec$eta1
 th[iter,2] <- ec$eta2
 bc <- beta2.calc(x,y,n,j0,k0,ec$eta1,ec$eta2)
 th[iter,3:7] <- bc$B
 theta <- th[1:iter,]
 #delta <- abs(th[iter,]-th[m,])
 delta <- abs(th[iter,]-th[m,])/th[m,]
 if( (delta[1]<.001) & (delta[2]<.001) & (delta[3]<.001)
   & (delta[4]<.001) & (delta[5]<.001) & (delta[6]<.001)
   & (delta[7]<.001) )
 break
 }
 list(theta=theta)
}
