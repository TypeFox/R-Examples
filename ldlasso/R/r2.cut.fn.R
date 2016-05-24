r2.cut.fn <-
function( block.obj, block.cood, Xa = NA, Y = NA, maxcol = 5000, r2.cut.min = 0, r2.cut.max = 0.5, r2.vec.length = 25 ){

if( isS4(block.obj) ){
    Xa <- as.double.snp.data(block.obj@gtdata)
    Y <- block.obj@phdata$dx
  }
  
  r2.cut.vec <- seq(r2.cut.min, r2.cut.max, length.out=r2.vec.length)
  ncol.vec <- c(); nrow.vec <- c();
  for( r2.cut in r2.cut.vec ){
    lasso.obj <- ld_lasso( block.obj = NA, Xa = Xa, Y = Y, s1 = 0, s2 = 0, r2.cut = r2.cut, form = 3, ytype = 1, solve = FALSE, block.cood = block.cood )
    ncol.vec <- c( ncol.vec, ncol(lasso.obj$A) )
    nrow.vec <- c( nrow.vec, nrow(lasso.obj$A) )
  }
  r2.cut <- min( r2.cut.vec[ncol.vec < maxcol] )
}

