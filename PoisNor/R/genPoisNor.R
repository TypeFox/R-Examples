genPoisNor <-
function(n, no.norm,  no.pois, cmat.star, lamvec, sd.vec, mean.vec) {
  
  if ( length(lamvec) != no.pois ){
    stop("Dimension of the lambda vector does not match the number of Poisson variables!\n")
  } 
  
  if ( length(mean.vec) != no.norm ){
    stop("Dimension of the mean vector does not match the number of normal variables!\n")
  } 
  
  
  if (no.norm==0){
    
    if(!is.null(mean.vec) | !is.null(sd.vec) ){
      
      stop("mean.vec and sd.vec must be NULL when no.norm is 0!\n")
      
    }
      
  }else{
      
    if (min(sd.vec)<0){
      stop("Standard deviation values cannot be negative!\n")
    } 
    
    if ( length(sd.vec) != no.norm ){
      stop("Dimension of the sd vector does not match the number of normal variables!\n")
    } 
    
  }

    
  if (nrow(cmat.star)!=(no.pois+no.norm) | nrow(cmat.star)!=(no.pois+no.norm) ){
      stop("Dimension of cmat.star and number of variables do not match!\n")
  }
  
  if(!is.positive.definite(cmat.star)){
    stop( "Intermediate correlation matrix is not positive definite!\n")
  }
    
  if (no.pois==0){
    YY = rmvnorm(n, rep(0,ncol(cmat.star) ), cmat.star)
    YY = t( t(YY)*sd.vec+mean.vec )
  }
  
  if (no.norm==0){
      
    XX=rmvnorm(n, rep(0,(no.pois)), cmat.star)
    YY=NULL;
    
    for( i in 1:length(lamvec)  ) {
      UU = pnorm(XX[,i])
      XXpois = qpois(UU,lamvec[i])
      YY=cbind(YY, XXpois) 
    }
  }
  
  if (no.norm>0 & no.pois>0) {

  	XX=rmvnorm(n, rep(0,(no.pois+no.norm)), cmat.star)
  	YY=NULL;
  
  	for( i in 1:length(lamvec)  ) {
  		UU = pnorm(XX[,i])
  		XXpois = qpois(UU,lamvec[i])
  		YY=cbind(YY, XXpois) 
  	}
  
  	YY = cbind(YY, XX[, (length(lamvec)+1) : ncol(cmat.star) ] )
  	rm(XX)

  	YY[,(length(lamvec)+1) : ncol(cmat.star)] = t( t(YY[,(length(lamvec)+1) : ncol(cmat.star)])*sd.vec+mean.vec )
  }
  colnames(YY)<-NULL
	return(YY)
}
