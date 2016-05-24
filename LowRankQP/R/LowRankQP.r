########## R function: LowRankQP  ##########

# For low rank quadratic programming.

# Last changed: 13 SEP 2005

LowRankQP <-function(Vmat,dvec,Amat,bvec,uvec,method="PFCF",verbose=FALSE,niter=200)
{
   # Some Type Checking
   typeError <- FALSE
   if ( nrow(Vmat)!=length(dvec) )
   {
        print("ERROR: nrow(Vmat)!=length(dvec)")
        typeError <- TRUE
   }
   if ( nrow(Vmat)!=ncol(Amat) )
   {
        print("ERROR: nrow(Vmat)!=ncol(Amat)")
        typeError <- TRUE
   }
   if ( nrow(Vmat)!=length(uvec) )
   {
        print("ERROR: nrow(Vmat)!=length(uvec)")
        typeError <- TRUE
   }
   if ( nrow(Amat)!=length(bvec) )
   {
        print("ERROR: nrow(Amat)!=length(bvec)")
        typeError <- TRUE
   }
   if (typeError) stop("ERROR: check input dimensions.")

   n <- nrow(Vmat)
   m <- ncol(Vmat)
   p <- nrow(Amat)
   
   alpha <- as.array(matrix( 0.0, n, 1 ))
   beta  <- as.array(matrix( 0.0, p, 1 ))
   xi    <- as.array(matrix( 0.0, n, 1 ))
   zeta  <- as.array(matrix( 0.0, n, 1 ))

   # Create numerical version of method for C call.

   if (method=="LU")   methodNum <- 1
   if (method=="CHOL") methodNum <- 2
   if (method=="SMW")  methodNum <- 3
   if (method=="PFCF") methodNum <- 4

   res <- .C("LowRankQP", n, m, p, as.integer(methodNum), as.integer(verbose),
         as.integer(200), Vmat, dvec, t(Amat), bvec, uvec, alpha, beta, xi, 
         zeta, PACKAGE="LowRankQP")

   alpha <- res[[12]]
   beta  <- res[[13]]
   xi    <- res[[14]]
   zeta  <- res[[15]]

   list(alpha=alpha, beta=beta, xi=xi, zeta=zeta)
}

######## End of LowRankQP ##########

