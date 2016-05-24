`Jacobian` <-
function(zeta){
if(length(zeta)==1) 
    return(1)
J<-JacobianK(zeta,1)
if (length(zeta)==2) 
    return(J)
for (j in 2:(length(zeta)-1)){
         J<-J%*%JacobianK(zeta,j)
         }
J
}

