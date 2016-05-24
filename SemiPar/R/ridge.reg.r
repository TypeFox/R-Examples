###### R-function: ridge.reg ##########

# For computing a generalized ridge regression
# fit for a given X matrix, y vector, weight
# vector and ridge vector. The C=[X Z] notation 
# for mixed model representation is used in this 
# function.

# Last changed: 21 SEP 2004

ridge.reg <- function(C.mat,y,weights=rep(1,nrow(as.matrix(y))),ridge.vec)
{
   # Form the augmented C.mat-matrix and y-matrix.

   sqrt.wts <- sqrt(as.vector(weights))
 
   Ca <- rbind(sqrt.wts*C.mat,diag(sqrt(ridge.vec)))

   ya <- rbind(as.matrix(sqrt.wts*y),matrix(0,ncol(C.mat),ncol(as.matrix(y))))

   # Solve this using QR methods.

   qr.fit.out <- lm.fit(Ca,ya)

   return(list(coef=qr.fit.out$coef,qr=qr.fit.out$qr,ridge.vec=ridge.vec))
}

########## End of ridge.reg ##########
