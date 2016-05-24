"bisect" <-
function(t,za,zb,target,drft=0,upper=FALSE){
   tol <- 0.000001
   dl <- 0.25
   gotlo <- 0
   gothi <- 0
   prev <- 0
   pr <- 0
   while ({abs(pr-target) > tol}&{abs(drft-prev) > tol/10}){
      glan.out <- glan(t,za,zb,drft)
      if (upper){
         pr <- sum(glan.out$qpos)
       }
      if (!upper){
         pr <- glan.out$pr
       }
      if (pr > target+tol){
         hi <- drft
         drft <- drft-dl
         gothi <- 1
       }
      if (pr < target-tol){
         lo <- drft
         drft <- drft+dl
         gotlo <- 1
       }
      if ({gothi==1}&{gotlo==1}){
         prev <- drft
         drft <- (lo+hi)/2
       }
    }
   if ({abs(drft-prev) <= tol/10}&{abs(pr-target) > tol}){
      warning("Convergence problem")
    }
   return(drft)
 }

