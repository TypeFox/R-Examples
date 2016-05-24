p0tos1 <-
function( p0, block.obj, Xa = NA, Y = NA, r2.cut = 0.01, s1high, s1low, max.iter = 1e2, tol = 1e-1 ){
 if(isS4(block.obj)){
   Xa <- as.double.snp.data(block.obj@gtdata)
   Y <- block.obj@phdata$dx[sample( length(block.obj@phdata$dx), size = length(block.obj@phdata$dx), replace = FALSE )]
 }
 
flag <- 0; s1old <- 1e6;
s1 <- s1high; iter <- 0

while( !flag ){

  lasso.obj <- ld_lasso( block.obj = NA, Xa = Xa, Y = Y, s1 = s1, s2 = NULL, r2.cut = r2.cut, form = 1, ytype = 1, block.cood = NA, solve = TRUE )
  beta <- lasso.obj$qp$solution[1:ncol(Xa)]
  p.lasso <- sum(abs(beta)>1e-6)

  if( sum( abs( lasso.obj$y ) > 1e-6 ) < p0 ){
    return( NA )
  }
  
  if( p.lasso == p0 ){
    flag <- 1
  }else if( p.lasso < p0 ){
    s1low <- s1
    s1 <- mean(c(s1low,s1high))
  }else if( p.lasso > p0){
    s1high <- s1    
    s1 <- mean(c(s1low,s1high))
  }

  iter <- iter+1
  if( iter == max.iter ){
    return( NA )
  }
  if( abs(s1-s1old)<tol ){
    return(s1)
  }
  s1old <- s1
}

return( s1 )

}

