ld_lasso_method <-
function( block.obj, block.cood = NA, Xa = NA, Y = NA, bpmap = NA, maxcol = 5e3, p.frac = 0.10, B = 5, s2low = 5e-3, s2high = 5e1, s2.vec.length = 4, null = FALSE ){
  if( isS4(block.obj)){
  
    Xa <- as.double.snp.data(block.obj@gtdata)
    Y <- block.obj@phdata$dx
    bpmap <- (block.obj@gtdata@map-block.obj@gtdata@map[1])/1e3
  }
  

  if( null ){
    Y <- Y[sample( length(Y), size = length(Y), replace = FALSE )]
  }
  
  block.bounds.vec <- block.bounds( map = bpmap , block.cood = block.cood )
  cat( "finding r2 cutoff...\n" )
  r2.cut <- r2.cut.fn( block.obj = NA, Xa = Xa, Y = Y, block.cood = block.cood, maxcol = maxcol, r2.cut.min = 0, r2.cut.max = 0.5, r2.vec.length = 25 )
  cat( "r2 cutoff = ", r2.cut, "\n", sep = "" )
  cat( "finding s1...\n" )

  if(isS4(block.obj)){
    s1 <- get.s1( i = floor(p.frac*ncol(Xa)), block.obj, r2.cut = r2.cut, s1high = 10, s1low = 1e-6, max.iter = 1e3, tol = 1e-3 )
  }else{
    s1 <- get.s1( i = floor(p.frac*ncol(Xa)), block.obj, Xa = Xa, Y = Y, r2.cut = r2.cut, s1high = 10, s1low = 1e-6, max.iter = 1e3, tol = 1e-3 )
  }


  cat( "s1 = ", s1, "\n", sep = "" )
  cat( "finding s2 " )

  if( isS4(block.obj)){
    cp.obj <-  get.cp( block.obj = block.obj, s2low = s2low, s2high = s2high, s2.vec.length = s2.vec.length, s1 = s1, r2.cut = r2.cut, block.cood = block.cood, B = B )
  }else{
    cp.obj <-  get.cp( block.obj = NA, Xa = Xa, Y = Y, s2low = s2low, s2high = s2high, s2.vec.length = s2.vec.length, s1 = s1, r2.cut = r2.cut, block.cood = block.cood, B = B )
  }

  s2star <- min(cp.obj$s2.vec[cp.obj$cp.vec==min(cp.obj$cp.vec, na.rm = TRUE)])
  cat( "\ns2 = ", s2star, "\n", sep = "" )
  cat( "solving ld lasso...\n" )
  if( isS4(block.obj)){
    ld.lasso.obj1 <- ld_lasso( block.obj = block.obj, s1 = s1, s2 = s2star, r2.cut = r2.cut, form = 3, ytype = 1, block.cood = block.cood )
    ld.lasso.obj2 <- ld_lasso( block.obj = block.obj, s1 = s1, s2 = min(cp.obj$s2.vec, na.rm = TRUE), r2.cut = r2.cut, form = 3, ytype = 1, block.cood = block.cood )
    ld.lasso.obj3 <- ld_lasso( block.obj = block.obj, s1 = s1, s2 = max(cp.obj$s2.vec, na.rm = TRUE), r2.cut = r2.cut, form = 3, ytype = 1, block.cood = block.cood )
  }else{
    ld.lasso.obj1 <- ld_lasso( block.obj = NA, Xa = Xa, Y = Y, s1 = s1, s2 = s2star, r2.cut = r2.cut, form = 3, ytype = 1, block.cood = block.cood )
    ld.lasso.obj2 <- ld_lasso( block.obj = NA, Xa = Xa, Y = Y, s1 = s1, s2 = min(cp.obj$s2.vec, na.rm = TRUE), r2.cut = r2.cut, form = 3, ytype = 1, block.cood = block.cood )
    ld.lasso.obj3 <- ld_lasso( block.obj = NA, Xa = Xa, Y = Y, s1 = s1, s2 = max(cp.obj$s2.vec, na.rm = TRUE), r2.cut = r2.cut, form = 3, ytype = 1, block.cood = block.cood )
  }
    
  beta1 <- ld.lasso.obj1$qp$solution[1:ncol(Xa)]
  beta2 <- ld.lasso.obj2$qp$solution[1:ncol(Xa)]
  beta3 <- ld.lasso.obj3$qp$solution[1:ncol(Xa)]
  chi2.vec <- chi2( Xa = Xa, Y = Y)
  log10p <- log10( pchisq( chi2.vec, lower.tail = FALSE, df = 1 ) )
  cat( "finished\n" )
  return( list( beta1 = beta1, beta2 = beta2, beta3 = beta3, s2star = s2star, cp.obj = cp.obj, log10p = log10p, bpmap = bpmap, block.bounds.vec = block.bounds.vec, s1 = s1, B = B, s2.vec.length = s2.vec.length ) )
}

