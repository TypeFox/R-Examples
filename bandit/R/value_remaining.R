value_remaining <- function(x,
  n,
  alpha = 1,
  beta = 1,
  ndraws = 10000
 ){
   post = sim_post(x,n,alpha,beta,ndraws)
   postWin = prob_winner(post)
   iMax = which.max(postWin)
   thetaMax = apply(post,1,max)
   #value_remaining:
   vR = (thetaMax-post[,iMax])/post[,iMax]
   return(vR)
 }
