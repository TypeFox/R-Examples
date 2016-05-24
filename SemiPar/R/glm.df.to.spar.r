
########## S function: glm.df.to.spar ##########

# Obtains the smoothing parameter corresponding
# to a specified degrees of freedom 
# for penalized likelihood estimation 
# Assumes single predictor

# Last changed: 17 JUL 2002

  glm.df.to.spar <- function(df,y,X,Z,family)
   {
      # Get initial smoothing parameter from
      # scatterplot smooth with df=4
        spar.init <- df.to.spar(df=4,X,Z)

      # Get corresponding fitted values 
        C.mat <- cbind(X,Z)        
        CTC <- t(C.mat)%*%C.mat
        D.mat <- spar.init*diag(c(rep(0,ncol(X)),rep(1,ncol(Z))))
   
        R <- chol(CTC+D.mat)

        R.inv <- backsolve(R,diag(rep(1,nrow(R))))   
        DRinv <- D.mat%*%R.inv 
        M <- t(R.inv)%*%DRinv

        svd.out <- svd(M)
        s.vec <- svd.out$d
        U <- svd.out$u

        A <- C.mat%*%solve(R,U)
        B <- t(U)%*%solve(t(R),t(C.mat))

        By <- B%*%y
        r.mat <- 1/(1+spar.init*s.vec)

        mu.fitted <- t(as.vector(By)*t(A))%*%r.mat   
        mu.fitted <-  gen.corr.range(mu.fitted,family)

        # Transform to get weight matrix

        if (family=="poisson")
             Wt <- mu.fitted

        if (family=="binomial")
             Wt <- mu.fitted*(1-mu.fitted)

        # Tranform X and Z by multiplying by sqrt of weight matrix

        Wt.half <- diag(as.vector(sqrt(Wt)))
  
        Wt.half.X <- Wt.half %*%X
        Wt.half.Z <- Wt.half %*%Z

        spar.final <- df.to.spar(df,Wt.half.X,Wt.half.Z)
        #spar.final <- sqrt(spar.final)

        return(spar.final)
    }

########## End of glm.df.to.spar ##########

