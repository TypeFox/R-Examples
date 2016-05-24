CMPStdErrors <-
function(x,beta,nu,max=100){

   weight <- weights(x,beta,nu,max)
   Ibeta <- InfoMatrix.beta(x,weight)
   Ibetanu <- InfoMatrix.betanu(x,beta,nu,max)
   Inu <- InfoMatrix.nu(x,beta,nu,max)
   info <- InfoMatrix(Ibeta,Ibetanu,Inu)
   CMP.SEs <- sqrt(diag(solve(info)))

return(CMP.SEs)

}

