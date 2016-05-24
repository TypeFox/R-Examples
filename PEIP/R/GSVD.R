GSVD <-
function(A,B)
{
    ##  the author thanks Berend Hasselman  and Lapack authors
    ## for help in preparing some of these wrappers
    ## some of the code here has been adapted from package
    ## geigen, by Hasselman

    if(!is.matrix(A)) stop("Argument A should be a matrix")
    if(!is.matrix(B)) stop("Argument B should be a matrix")
  
    
    # A=U*E1*Q'
    # B=V*E2*Q'
  

    dimA = dim(A)
    dimB = dim(B)
    if(dimA[1]==0) stop("Matrix A has zero rows/columns")
    if(dimB[1]==0) stop("Matrix B has zero rows/columns")
    
    if(!all(is.finite(A))) stop("Matrix A may not contain infinite/NaN/NA")
    if(!all(is.finite(B))) stop("Matrix B may not contain infinite/NaN/NA")

    ncolA = dimA[2]
    nrowA = dimA[1]
    ncolB = dimB[2]
    nrowB = dimB[1]

        
    lwork = max(c(3*ncolA,nrowA,nrowB))+ncolA

    ju = 2
    jv = 2
    jq = 2

  
    z <- .Fortran("zdggsvd",
                  as.integer(ju),
                  as.integer(jv),
                  as.integer(jq),
                  as.integer(nrowA),
                  as.integer(ncolA),
                  as.integer(nrowB),
                  integer(1),
                  integer(1),
                  as.double(A),
                  as.integer(nrowA),
                  as.double(B),
                  as.integer(nrowB),
                  double(ncolA),
                   double(ncolA),
                   double(nrowA*nrowA),
                  as.integer(nrowA),
                  double(nrowB*nrowB),
                  as.integer(nrowB),
                   double(ncolA*ncolA),
                  as.integer(ncolA),
                  double(lwork),
                  integer(ncolA),
                  integer(1) ,dup=FALSE,
                  PACKAGE="PEIP")

    K=z[7][[1]]
    L=z[8][[1]]
    U=z[15][[1]]
    V=z[17][[1]]


    Q=z[19][[1]]
    ALPHA=z[13][[1]]

        BETA=z[14][[1]]

    R=matrix(z[9][[1]],ncol(A),nrow=nrow(A),byrow=FALSE)

     
    U=matrix(U,ncol=nrow(A),nrow=nrow(A),byrow=FALSE)
     
    V=matrix(V,ncol=nrow(B),nrow=nrow(B),byrow=FALSE)
     
      Q=matrix(Q,ncol=ncol(A),nrow=ncol(A),byrow=FALSE)
     

    D1=mat.or.vec(nrow(A),K+L)
    D2=mat.or.vec(nrow(B),K+L)

     

    oR=mat.or.vec((K+L),ncol(A))


 ###    ORDALPHA = order(ALPHA)
 ###    ORDBETA = order(BETA)

###ALPHA = ALPHA[ORDALPHA]
###BETA  = BETA[ORDBETA]

     
###U = U[ ,ORDALPHA]
     
 ###    V = V[ ,ORDBETA[ORDBETA<=nrow(B)]  ]
     
    if(K > 0)
    {

        if(K==1)
    { D1[1:K,1:K] =rep(1,K)
    }
    else
    {
      diag(D1[1:K,1:K])=rep(1,K)
    }

         diag(D1[(K+1):(K+L),(K+1):(K+L)])=ALPHA[(K+1):(K+L)]
        diag(D2[1:L,(K+1):(K+L)])=BETA[(K+1):(K+L)]
               

    }

    if(K ==0)
    {

        diag(D1[(K+1):(K+L),(K+1):(K+L)])=ALPHA[(K+1):(K+L)]
        diag(D2[1:L,(K+1):(K+L)])=BETA[(K+1):(K+L)]


    }    

    Ci=ALPHA[(K+1):(K+L)]
    S=BETA[(K+1):(K+L)]



     
    oR[(1):(K+L),(ncol(A)-K-L+1):(ncol(A))]=R[(1):(K+L),(ncol(A)-K-L+1):(ncol(A))]    



    X = t(oR  %*% t(Q))

     #############  matlab style return
     return(list(U=U,V=V, X=X, C=D1,S=D2) )

            
  ##   return(list(U=U,V=V,Q=Q,D1=D1,D2=D2,oR=oR,C=Ci,S=S,K=K,L=L,Z=z))

     
}
