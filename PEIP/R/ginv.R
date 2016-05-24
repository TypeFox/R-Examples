ginv <-
function(G, x, tol=1e-12 )
  {
    ###  need something to make up for the lame-o
    ##   matlab code that does this h = G\x to get the inverse
    
    return(solve(t(G)%*%G , t(G) %*% x , tol)) 

    
  }
