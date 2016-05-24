###### S-function: glmeAux ##########

# For obtaining df and covariance matrix
# of fixed and random parameter estimates and predictions in 
# genealized linear mixed models.

# Last changed: 29 JUL 2002

glmeAux <- function(X,Z,G,block.inds=NA,ridge.reg.fit,family)
{   
   # Set the sample size, n.

   n <- nrow(X)
   
   # Obtain the degrees of freedom of the overall fit, starting
   # with the  Q1-matrix.

   Q1 <- qr.Q(ridge.reg.fit$qr)[1:n,]
   
   Q1TQ1 <- t(Q1)%*%Q1

   df.fit <- sum(diag(Q1TQ1))

   # Obtain  error degrees of freedom.

   df.var <- sum(diag(Q1TQ1%*%Q1TQ1))

   df.res <- n - 2*df.fit + df.var

   # Extract further information about the model

   # Obtain the C-matrix and R.inv-matrix

   R <- qr.R(ridge.reg.fit$qr)

   Cmat.w <- Q1%*%R

   I.mat <- diag(rep(1,nrow(R)))

   R.inv <- backsolve(R,I.mat)

   # Obtain the unscaled covariance matrix, starting
   # with the ridge vector.

   ridge.vec <- ridge.reg.fit$ridge.vec

   cov.mat <- R.inv %*%(I.mat - t(R.inv)%*%(ridge.vec*R.inv))%*%t(R.inv)

   # Calculate df for jth component of model

   df <- numeric()
   for (j in 1:length(block.inds))
   {
      curr.inds <-  block.inds[[j]]
      df[j] <- sum(diag(R.inv[curr.inds,]%*%t(Q1)%*%Cmat.w[,curr.inds]))
   }

   # Package auxiliary quantities object and return

   glmeAux.object <- list(cov.mat=cov.mat,df=df,
                          block.inds=block.inds,random.var=G,
                          df.fit=df.fit,df.res=df.res)

   return(glmeAux.object)

}

########## End of glmeAux ##########


