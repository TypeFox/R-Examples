ld_lasso <-
function( block.obj, block.cood = NA, Xa = NA, Y = NA, s1, s2, r2.cut = 0.5, delta = 1e-10, form = 3, ytype = 1, solve = TRUE ){

  if( isS4(block.obj) ){
    Xa <- as.double.snp.data(block.obj@gtdata)
    Y <- block.obj@phdata$dx
  }

  if ( dim(Xa)[2] <= 2  ){
        return()
  }else{
	p <- dim(Xa)[2]
        if( is.na(block.cood[1]) ){
          block.cood <- c(1,rep(0,ncol(Xa)-1),1)
        }
        index.mat <- which(cor(Xa)^2 > r2.cut & lower.tri(matrix(1, p, p)) & block.map.matrix(block.cood), arr.ind = TRUE )
	r <- dim(index.mat)[1]
        if( r == 0 ){
          return("r is zero")
        }else{
          D <- matrix( 0, nrow = r, ncol = p )
          r2.vec <- c()
          for( i in 1:r ){
		D[i,index.mat[i,1]] <- -1
		D[i,index.mat[i,2]] <- 1
                r2.vec <- c( r2.vec, cor(Xa[,index.mat[i,1]],Xa[,index.mat[i,2]])^2 )
          }
          if( ytype == 1 ){
          X2 <- Xa[Y == 1,] # 1 is case
	  X1 <- Xa[Y == 0,] # 0 is control
          f2 <- colSums(X2)/2/dim(X2)[1]
	  f1 <- colSums(X1)/2/dim(X1)[1]
	  n0 <- sum( Y == 0 )
	  n1 <- sum( Y == 1 )
          OR <- ifelse( ( f2 == 0 & f1 == 0) | (f2 == 1 & f1 == 1), 1, f2/(1-f2)/(f1/(1-f1)) )
	  OR <- ifelse( OR == Inf, 1e6, OR )
	  OR <- ifelse( OR == 0, 1/1e6, OR )
          y <- log(OR)
          var_y <- ifelse( f1 == 0 | f2 == 0, 1e6, ( n0*f1*(1-f1) + n1*f2*(1-f2) ) / ( 2*n0*n1*f1*f2*(1-f1)*(1-f2) ) )
          y <- y/sqrt(var_y)
        }else if ( ytype == 2 ){
          y <- Y
          
        }
          r2 <- ifelse( r2.vec == 0, 1e-6, r2.vec )
          if( form == 1 ){
            cpcc.vec <- 1e6*rep(1,length(r2))
          }else if( form == 2 ){
            cpcc.vec <- -s2*log(r2) + delta
            s1 <- 1e6
          }else if( form == 3 ){
            cpcc.vec <- -s2*log(r2) + delta
          }else{
            stop("bad constraint form")
          }
          #cpcc.vec <- rep(10, length(r2))
          One <- rep( 1, p )
          I <- diag(1, nrow = p)
          Zero <- rep( 0, p)
          ZeroMat <- matrix( 0, nrow = p, ncol = p )
          ZeroMat_minus1 <- matrix( 0, nrow = p - 1, ncol = p )
          ZeroMat_r <- matrix( 0, nrow = r, ncol = p )
          A1 <- cbind( I, -I, I )
          A2 <- cbind( ZeroMat, I, ZeroMat )
          A3 <- cbind( ZeroMat, ZeroMat, I )
          A4 <- c( Zero, -One, -One )
          A5 <- cbind( ZeroMat_r, D, D )
          A6 <- cbind( ZeroMat_r, -D, -D )
          I_3p <- diag( 1, nrow = 3*p )
          yc = c( y, ifelse( y >= 0, y, 0 ), ifelse( y <= 0, -y, 0 ) )
          A <- t(rbind( A1, A2, A3, A4, A5, A6 ))
          b0 <- c( rep( Zero, 3 ), -s1, rep( -cpcc.vec , 2 )  )
          if( solve ){
            qp <- solve.QP(Dmat = I_3p, dvec = yc, Amat = A, bvec = b0, meq = p, factorized = FALSE)
            return(list( qp = qp, y = y, A = A, r2 = r2, b0 = b0, OR = OR ))
          }else{
            return(list( qp = NULL, y = y, A = A, r2 = r2, b0 = b0, OR = OR ))
          }
          
        }
   }
}

