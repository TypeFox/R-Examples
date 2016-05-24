##########################################
#This function computes the log of a     #
#inverse Wishart distribution            #
##########################################

logdiwish <-
function(W, v, S)
   {
   #This function is from the package MCMCpack by  
   #Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park
   if(is.pos.def(W))
     {
     k <- nrow(S)
     gammapart <- 1
     for (i in 1:k) {
       gammapart <- gammapart * gamma((v + 1 - i)/2)
       }
     denom <- gammapart * 2^(v * k/2) * pi^(k * (k - 1)/4)
     detS <- det(S)
     detW <- det(W)
     hold <- S %*% solve(W)
     tracehold <- sum(hold[row(hold) == col(hold)])
     num <- detS^(v/2) * detW^(-(v + k + 1)/2) * exp(-1/2 * tracehold)
     ans <- log(num)-log(denom)
     }
   else
     {ans <- -Inf}
   return(ans)
   }

