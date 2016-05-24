get.cp <- function( s2low, s2high, s2.vec.length, block.obj, Xa = NA, Y = NA, s1, r2.cut, block.cood, B = 20 ){

  if( isS4(block.obj)){
    Xa <- as.double.snp.data(block.obj@gtdata)
    Y <- block.obj@phdata$dx
  }
  
s2.vec <- 10^( seq(log10(s2low), log10(s2high), length.out = s2.vec.length) )

p <- ncol(Xa)
jj <- 0; Ntot <- length(s2.vec)*B
cp.vec <-  c(); beta0.mat <- c()

for( s2 in s2.vec ){
  
  if(isS4(block.obj)){
    ld.lasso.obj <- ld_lasso( block.obj = block.obj, s1 = s1, s2 = s2, r2.cut = r2.cut, form = 3, ytype = 1, block.cood = block.cood )
  }else{
    ld.lasso.obj <- ld_lasso( block.obj = NA, Xa = Xa, Y = Y, s1 = s1, s2 = s2, r2.cut = r2.cut, form = 3, ytype = 1, block.cood = block.cood )
  }

  y0 <- ld.lasso.obj$y
  beta0 <- ld.lasso.obj$qp$solution[1:ncol(Xa)]
  beta0.mat <- rbind( beta0.mat, beta0 )

  ystar.mat <- c()
  betastar.mat <- c()

  for( b in 1:B ){

    jj <- jj + 1

    if( Ntot > 10 ){
      if( jj%%floor(Ntot/10) == 0 ){
        cat( c( 100*round(jj/Ntot,2), "% " ), sep = "" )
      }
    }
    
    boot.indx <- sample(x = nrow(Xa), size = nrow(Xa), replace = TRUE)
    Xstar <- Xa[boot.indx,]
    Ystar <- Y[boot.indx]

    ld.lasso.boot.obj <- ld_lasso( block.obj = NA, Xa = Xstar, Y = Ystar, s1 = s1, s2 = s2, r2.cut = r2.cut, form = 3, ytype = 1, block.cood = block.cood )
    ystar <- ld.lasso.boot.obj$y
    betastar <- ld.lasso.boot.obj$qp$solution[1:ncol(Xstar)]
    ystar.mat <- rbind( ystar.mat, ystar )
    betastar.mat <- rbind( betastar.mat, betastar )

  }

  a <- ( y0 - beta0 )%*%( y0 - beta0 )

  df <- 0
  for( j in 1:ncol(ystar.mat) ){
    df <- df + cov(ystar.mat[,j], betastar.mat[,j])
  }
  cp <- a - p + 2*df
  cp.vec <- c( cp.vec, cp )
}
return(list( s2.vec = s2.vec, cp.vec = cp.vec, beta0.mat = beta0.mat, s1 = s1 ) )

}
