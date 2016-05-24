################################################################################
## manylm.multiwfit                                                           ## 
################################################################################

manylm.multiwfit <- function(x, y, w, offset = NULL, tol=1.0e-050,
  singular.ok = TRUE, ...)
{


 y <- as.matrix(unabund(y))
 w <- as.matrix(w)
 
if (is.null(n <- nrow(x)))
        stop("'x' must be a matrix")
    if (n == 0)
        stop("0 (non-NA) cases in 'x'")
    p <- ncol(x)
    
    if (p == 0) {
        return(list(coefficients = numeric(0), residuals = y,
            fitted.values = 0 * y, rank = 0, df.residual = NROW(y)))
    }
    
    ny <- ncol(y)

    if (!is.null(offset)) y <- y - offset
    if (nrow(y) != n)
        stop("dimensions of 'x' and 'y' do not match")

    if( ncol(y) != ncol(w) )  stop("dimensions of 'y' and 'w' do not match")

    if (length(list(...))){
	  tmp <- paste(names(list(...)), sep = ", ")
        warning("extra arguments ", tmp ,
            " are just disregarded.")
	}

    storage.mode(x) <- "double"
    storage.mode(y) <- "double"

    ok      <- w == 0
    n.ok    <- n - colSums( w == 0 )
    if( any(n.ok == 0) )
      stop("No informative observations in the weights for variable ",
    which(n.ok==0), ". All 'weights' are equal to Zero.")
    
################### BEGIN Fit multivariate linear model  ##################
####### does the same as model.fit, but I can return txX and hat.X ########

results         <- list()
results$assign  <- attr(x, "assign")
xnames          <- colnames(x)

if(is.null( xnames ))
    xnames <- paste("x", 1:p, sep = "")

qrx <- list()
qrx[[1]]  <- qr( x * sqrt(w[,1]) )

pivot <- qrx[[1]]$pivot

rank        <- qrx[[1]]$rank
p1          <- 1:rank

if (rank < p & !singular.ok ) 
      stop("singular fit encountered") 
  
x.ind <- x[,pivot[p1], drop=FALSE]

coefficients  <- matrix(nrow=p,ncol=ny)

txy      <-  t( x.ind ) %*% ( y*w )

chol2invqrx  <-  list()

hat.X  <- txX <- list()

for( i in 1:(ncol(w)) ) {
   qrx[[i]]         <- qr( x * sqrt(w[,i]) )
   chol2invqrx[[i]] <- chol2inv( qrx[[i]]$qr[p1,p1, drop=FALSE] )
   coefficients[ pivot[p1], ] <-  chol2invqrx[[i]] %*% txy
   txX[[i]] <- t(x.ind*w[,i]) %*% x.ind
}

dimnames(coefficients) <- list(xnames,dimnames(y)[[2]])


fitted.values <- x %*% coefficients

residuals     <- y - fitted.values
# Ensure the numerical stability for weights=0 cases.
residuals[ok] <- y[ok]
colnames(residuals) <- colnames(y)

df.residual   <- n.ok - rank

################## END Fit multivariate linear model ##################

results$coefficients        <- coefficients
results$residuals           <- residuals
results$fitted.values       <- fitted.values
results$weights             <- w
results$rank                <- rank
results$qr                  <- qrx # QR decomposition of x * sqrt(w)
results$df.residual         <- df.residual

return(results)
}
 
