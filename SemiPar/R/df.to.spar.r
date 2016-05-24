########## S function: df.to.spar ##########

# Obtains the smoothing parameter corresponding
# to a specified degrees of freedom for a penalized
# spline scatterplot smooth.

# Last changed: 27 AUG 2001

df.to.spar <- function(df,X,Z)
{
   # Form C, C^TC and D matrices

   C.mat <- cbind(X,Z)

   CTC <- t(C.mat)%*%C.mat
   D <- diag(c(rep(0,ncol(X)),rep(1,ncol(Z))))

   # Check for df in valid range.

   if (df < ncol(X)) stop("df cannot be less than the number of columns in X")

   if (ncol(C.mat)<df)  stop("df cannot exceed the number of columns in (X,Z)")


   if (df==ncol(X)) return(Inf)

   if (ncol(C.mat)==df) return(0)


   # Use Cholesky/SVD idea to obtain fast df function

   R <- chol(CTC)

   R.inv <- backsolve(R,diag(rep(1,nrow(R))))   

   DRinv <- D%*%R.inv 

   M <- t(R.inv)%*%DRinv

   Lamvec <- svd(M)$d

   df.func <- function(arg,df,Lamvec)
      return(sum(1/(1+arg*Lamvec))-df)

   # Obtain approximation to root using Wand (1999)

   tr.G <- sum(diag(M))
      
   approx.root <- (ncol(C.mat)-df)/tr.G

   root.lower <- 0   

   func.lower <- df.func(root.lower,df,Lamvec)

   root.upper <- approx.root

   upper.found <- FALSE

   while (upper.found==FALSE)
   {
      func.upper <-  df.func(root.upper,df,Lamvec)

      if (func.lower*func.upper>0)
              root.upper <- 2*root.upper
      else
          upper.found <- TRUE

   }

   # Solve using uniroot()

   root <- uniroot(df.func,c(root.lower,root.upper),df=df,
                   Lamvec=Lamvec,tol = 0.00001*(.Machine$double.eps^.25))$root
  
   return(root)
}

########## End of df.to.spar ##########

