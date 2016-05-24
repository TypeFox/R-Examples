########### S-function: lmeAux #################

# For obtaining df and conditional (on variance components) variance 
# of fixed and random parameter estimates and predictions in linear
# mixed models.

# Last changed: 20 SEP 2004

lmeAux <- function(X,Z,G,resid.var,block.inds=NA,indiv=1,rho=0)
{
   # Required Arguments:
   #
   # X is the design matrix of the fix effects.
   # Z is the design matrix of the random effects.
   # G is the covariance matrix of the random effects.
   # resid.var is a scalar with the variance of the regression error.
   # (regression errors need to be iid to use this)
   #
   # block.inds is an optional argument containing a list of column 
   # indices. Each element of the list corresponds to a diffent component
   # in an addititve model.  NOTE THAT THIS LIST CORRESPONDS
   # TO THE COLUMNS OF C = [X,Z], NOT JUST THE COLUMNS OF Z. The default
   # for block.inds is list(1:number of columns of [X,Z]).
   #
   # Optional Arguments 
   # 
   # These should only be specified when an AR(1) structure has been
   # fit to var(e).
   #
   # indiv is the vector that numbers the individuals. It is
   #         assumed to be sorted into homogenous groupings. i.e.: 
   #         1 1 1 2 2 3 3 3 3 3 is OK
   #         1 2 3 1 2 3 1 3 3 3 is not
   # indiv should only be specified when an ar(1) structure on var(e) was fit
   #  
   # rho = AR(1) parameter
   #         rho is zero if it is left blank
   #
   # Mixed Model: y = X beta + Z * u + e 
   # We define G so that var(u) = G and resid.var = var(e_i)
   # 
   # lmeAux.object contains:
   #    cov.mat       - Var(\beta,\uhat-u)
   #    df.fit        - df of the whole fit
   #    df            - df of each additive component
   
   # Compute indices

   n <- nrow(as.matrix(X))

   if (!is.null(Z))
   {
      q <- ncol(as.matrix(Z)) 
      p <- ncol(as.matrix(X))

      # Calculate default block.inds if necessary

      if (is.list(block.inds)==0)
         block.inds <- list(1:(p+q))

      # Make appropriate matrices: C.mat is the "design" matrix

      if (qr(G)$rank == ncol(G))
         G.tilde <- rbind(matrix(0,p,p+q),
                        cbind(matrix(0,q,p),resid.var*solve(G)))
 
      if (qr(G)$rank < ncol(G))
         G.tilde <- rbind(matrix(0,p,p+q),
                        cbind(matrix(0,q,p),diag(rep(1000000,q))))

      C.mat  <- cbind(X,Z)
   }

   if (is.null(Z)) C.mat <- X

   CTC         <- t(C.mat)%*%C.mat
   CTRC        <- compute.CTRinvC(X,Z,indiv,rho)

   if(is.null(Z))
     G.tilde <- matrix(0,nrow(CTRC),ncol(CTRC))

   Ridge <- CTRC + G.tilde

   # Roundoff errors sometimes  
   # make Ridge and CTRC seem asymmetric

   Ridge       <- (Ridge + t(Ridge))/2
 
   # Compute the Cholesky factorization and the inverse

   R.ridge     <- chol(Ridge)   

   R.ridge.inv <- backsolve(R.ridge,diag(rep(1,nrow(R.ridge))))

   ridge.inv   <- R.ridge.inv%*%t(R.ridge.inv)
 
   # Compute the covariance matrix and the matrix for the df
   
   df.mat      <-   ridge.inv%*%CTC

   cov.mat     <-   resid.var*ridge.inv         
                
   # Calculate df for jth term in each component of model
   # Must first partition/index accordingly

   df <- numeric()
   for (j in 1:length(block.inds))
   {
      curr.inds <-  block.inds[[j]]
      df[j] <- sum(diag(as.matrix(df.mat[curr.inds,curr.inds])))
   }

   # Calculate df_fit

   df.fit <- sum(df)

   # Calculate df_res

   df.res <- n - 2*df.fit + sum(diag(df.mat%*%df.mat))

   # Form and return auxiliary quantities object

   lmeAux.object <- list(cov.mat=cov.mat,df=df,
                          block.inds=block.inds,
                          resid.var=resid.var,random.var=G,
                          df.fit=df.fit,df.res=df.res)

   return(lmeAux.object)
}

########## End of lmeAux #################
