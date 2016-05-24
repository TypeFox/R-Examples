# $Id: binomix.R 166 2014-07-13 21:41:15Z khliland $



binomixEstimate <- function( pan.matrix, K.range=3:5, core.detect.prob=1.0, verbose=TRUE ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  y <- table( factor( colSums( pan.matrix ), levels=1:dim( pan.matrix )[1] ) )
  bic.tab <- matrix( NA, nrow=length(K.range), ncol=3 )
  colnames( bic.tab ) <- c( "Core.size", "Pan.size", "BIC" )
  rownames( bic.tab ) <- paste( K.range, "components" )
  mix.list <- vector( "list", length( K.range ) )
  for( i in 1:length( K.range ) ){
    if( verbose ) cat( "binomixEstimate: Fitting", K.range[i], "component model...\n" )
    lst <- binomixMachine( y, K.range[i], core.detect.prob )
    bic.tab[i,] <- lst[[1]]
    mix.list[[i]] <- lst[[2]]
  }
  if( bic.tab[length(K.range),3] == min( bic.tab[,3] ) ) warning( "Minimum BIC at maximum K, increase upper limit of K.range" )
  binomix <- list( BIC.table=bic.tab, Mix.list=mix.list )
  class( binomix ) <- c( "Binomix", "list" )
  return( binomix )
}



binomixMachine <- function( y, K, core.detect.prob=1.0 ){
  n <- sum( y )
  G <- length( y )
  ctr <- list( maxit=300, reltol=1e-6 )
  np <- K - 1
    
  pmix0 <- rep( 1, np )/K            # flat mixture proportions
  pdet0 <- (1:np)/(np+1)              # "all" possible detection probabilities
  p.initial <- c( pmix0, pdet0 )      # initial values for parameters
  # the inequality constraints...
  A <- rbind( c( rep( 1, np ), rep( 0, np ) ), c( rep( -1, np ), rep( 0, np ) ), diag( np+np ), -1*diag( np+np ) )
  b <- c( 0, -1, rep( 0, np+np ), rep( -1, np+np ) )
  
  # The estimation, maximizing the negative truncated log-likelihood function
  est <- constrOptim( theta=p.initial, f=negTruncLogLike, grad=NULL, method="Nelder-Mead", control=ctr, ui=A, ci=b, y=y, core.p=core.detect.prob )
  
  estimates <- numeric( 3 )
  names( estimates ) <- c( "Core.size", "Pan.size", "BIC" )
  estimates[3] <- 2*est$value + log(n)*(np+K)                         # the BIC-criterion
  p.mix <- c( 1 - sum( est$par[1:np] ), est$par[1:np] )               # the mixing proportions
  p.det <- c( core.detect.prob, est$par[(np+1):length( est$par )] )   # the detection probabilities
  ixx <- order( p.det )
  p.det <- p.det[ixx]
  p.mix <- p.mix[ixx]
    
  theta_0 <- choose( G, 0 ) * sum( p.mix * (1-p.det)^G )
  y_0 <- n * theta_0/(1-theta_0)
  estimates[2] <- n + round( y_0 )
  ixx <- which( p.det >= core.detect.prob )
  estimates[1] <- round( estimates[2] * sum( p.mix[ixx] ) )
    
  mixmod <- matrix( c( p.det, p.mix ), nrow=2, byrow=T )
  rownames( mixmod ) <- c( "Detection.prob", "Mixing.prop" )
  colnames( mixmod ) <- paste( "Comp_", 1:K, sep="" )

  return( list( estimates, mixmod ) )
}



negTruncLogLike <- function( p, y, core.p ){
  np <- length( p )/2
  p.det <- c( core.p, p[(np+1):length(p)] )
  p.mix <- c( 1-sum( p[1:np] ), p[1:np] )
  G <- length( y )
  K <- length( p.mix )
  n <- sum( y )
    
  theta_0 <- choose( G, 0 ) * sum( p.mix * (1-p.det)^G )
  L <- -n * log( 1 - theta_0 )
  for( g in 1:G ){
    theta_g <- choose( G, g ) * sum( p.mix * p.det^g * (1-p.det)^(G-g) )
    L <- L + y[g] * log( theta_g )
  }
  return( -L )
}


plot.Binomix <- function( x, type="pan", cex=2, ncomp=NA, ... ){
  Binomix <- x
  if( is.na(ncomp) ){
    ncomp <- which( Binomix$BIC.table[,3] == min( Binomix$BIC.table[,3] ) )[1]
  } else {
    ncomp <- which( as.numeric( gsub( " components", "", rownames( Binomix$BIC.table ) ) ) == ncomp )
    if( length( ncomp ) != 1 ) stop( "Specified value of ncomp has not been fitted for this model" )
  }
  layout( matrix( c(1,1,1,1,1,1,2), nrow=1 ) )
  cpar <- par()$mar
  par( mar=c(2,2,2,2) )
  typ <- grep( type, c( "pan", "single" ) )
  fit <- Binomix$Mix.list[[ncomp]]
  dprob <- as.character( round( fit[1,]*1000 )/1000 )
  if( length( typ ) == 0 ){
    stop( "Unknown type specified" )
  } else if( typ == 1 ){
    pie( fit[2,], clockwise=T, col=panColor( fit[1,] ), labels=dprob, radius=1.0, cex=cex )
  } else {
    eg <- fit[2,]*fit[1,]
    pie( eg/sum( eg ), clockwise=T, col=panColor( fit[1,] ), labels=dprob, radius=1.0, cex=cex )
  }
  
  par( mar=c(1,1,1,3) )
  p <- (0:100)
  plot( rep( 0, 101 ), p, cex=0, col="white", xlim=c(0,1), ylim=c(0,100), 
        xaxt="n", yaxt="n", xlab="", ylab="" )
  box( lwd=3, col="white" )
  cols <- panColor( p/100 )
  xl <- rep( 0, 100 )
  xr <- rep( 1, 100 )
  yb <- (0:100)
  yt <- 1:101
  rect( xl, yb, xr, yt, col=cols, border=cols )
  axis( side=4, at=seq(0,100,10), labels=as.character( seq(0,100,10)/100) )
  par( mar=cpar )
}

panColor <- function( p ){
  level <- pretty( c(0,1), 100 )
  nlevel <- length( level )
  crp <- colorRampPalette( c("pink","orange","green","cyan","blue") )(nlevel)
  return( crp[1+round( 100*p )] )
}

summary.Binomix <- function( object, ... ){
  ncomp <- which( object$BIC.table[,3] == min( object$BIC.table[,3] ) )[1]
  cat( "Minimum BIC model at", rownames( object$BIC.table )[ncomp], "\nFor this model:\n" )
  cat( "Estimated core size:", object$BIC.table[ncomp,1], "clusters\n" )
  cat( "Estimated pangenome size:", object$BIC.table[ncomp,2], "clusters\n" )
}

str.Binomix <- function( object, ... ){
  cat( "Binomial mixture model object with", dim( object$BIC.table )[1], "fitted models\n" )
}

