svmmajtuning <- function(X,y, lambda= 1,
        weights.obs = 1, weights.var= 1, standardize = 'interval',
        spline.knots = 0, spline.degree = 1L,
        kernel = vanilladot, kernel.sigma=1 , kernel.degree=1L ,
        kernel.scale=1 , kernel.offset=0,
        hinge = 'absolute', hinge.k = 5,  eps=1e-8,increase.step=20L,
        convergence=1e-4,check.positive=TRUE,
                      ngroup=5,groups=NULL, na.action=na.omit,...)
  {
#PERFORMING CROSS VALIDATION USING A LIST OF PARAMETER SETTINGS
#INITIALISES ALL COMBINATION OF PARAMETER SETTINGS AND
#CALLS THE FUNCTION CROSSVAL.SVMMAJ
    if(!is.factor(y)) y <-factor(y)
    classes      <- sort(levels(y))
    y            <- sign( (y==classes[2]) - .5 )
    n            <-  length(y)

    if(!is.data.frame(X)){
    #OBTAIN THE DIMENSION OF X
        if(!is.matrix(X)) X  <- matrix(X,nrow=n)
        X    <- data.frame(X=X)
    }

    #===================================================
    #CHECKING FOR NONNEGATIVE VALUES
    #---------------------------------------------------
    #CHECK FOR NONPOSITIVE VALUES (IF NEEDED)
    if(check.positive){
        if(any(weights.obs<0)) stop('weights should be nonnegative')
        if(any(lambda<0))
                               stop('lambda should be nonnegative')
    }
   #====================================================
   #INITIALISING WEIGHTS OF OBJECTS
   #----------------------------------------------------
    if(!is.null(names(weights.obs)) && length(weights.obs)==2)
        weights.obs <- weights.obs[sort(names(weights.obs))]

    w <- rep(1,length(y))
    if(length(weights.obs)==2){
       w[y==-1] <- weights.obs[[1]]
       w[y==1]  <- weights.obs[[2]]
    } else if(length(weights.obs)==length(y))
       w[] <- weights.obs
    w <- matrix(w)

    #===================================================
    #INITIALIZE GROUP SIZE
    #---------------------------------------------------
    ngroup <- trunc(ngroup)
    if(is.null(groups)){
      if (ngroup < 2) stop("ngroup should be greater than or equal to 2")
      if (ngroup > n) stop("ngroup should be less than or equal to the number of observations")

      if (ngroup == n) {
      #Jackknife cross validation
          groups <- 1:n
          leave.out <- 1
      } else if (ngroup < n) {
      #Fixed number of groups
           leave.out <- trunc(n/ngroup)
           groups    <- sample(rep(1:ngroup, floor(n/ngroup) )[1:(n-leave.out)])
      }
    } else {
        #Prespecified groups
        ngroup    <- max(groups)
        leave.out <- trunc(n/ngroup)
    }
    #===================================================
    #PERFORM CROSS-VALIDATION
    #---------------------------------------------------
    #DEFINE HINGE FUNCTION
    newHinge          <- getHinge(hinge, hinge.k , eps = eps)

    #INITIALISE DATA TRANSFORMATION
    X			<- lapply(X,transformdata,standardize,spline.knots,spline.degree)
    prop.data	<- lapply(X,attributes)

    expansion   <- sapply(prop.data,`[[`,'dim')[2,]
    if(!is.null( unlist(sapply(prop.data,`[[`,'splineDegree')))){
    	spline.knots	<- max(0,sapply(X,function(x) length(attr(x,'splineInterval')))-1)
    	spline.degree	<- max(unlist(sapply(prop.data,`[[`,'splineDegree')))
    } else {
    	spline.knots=0
    	spline.degree=1
    }

    if(length(weights.var)==length(X))
    	X <- mapply(`*`,X,weights.var)
    X <- data.frame(X)
    X <- data.matrix(X)

    kernelvars <- list(sigma  = kernel.sigma,
                       degree = kernel.degree,
                       scale  = kernel.scale,
                       offset = kernel.offset)
    kernelvars <- kernelvars[names(formals(kernel))]
    G            <- expand.grid(kernelvars)
    nkernel      <- NROW(G)
    if(nkernel==0){ G=data.frame('(no param)'=NA); nkernel=1}
    lambda       <- sort(unique(lambda),decreasing=TRUE)
    nlambda      <- length(lambda)
    M            <- cbind(sapply(G,rep,each=nlambda),lambda=rep(lambda,nkernel))
    colnames(M) <- c( sapply(colnames(M)[1:length(kernelvars)],function(s) paste('kernel',sep='.',s)),
                        'lambda')
    cat(NROW(M),"gridpoints:\n")
    lMat <- vector('numeric',nkernel*nlambda)
    for( j in 1:nkernel  ){
      formals(kernel)<-G[j,,drop=FALSE]
      #DETERMINE EFFICIENT UPDATE
      method            <- getUpdate(X,kernel())
      theta             <- rep(0,method$x+1)
      #DEFINE MATRICES TO BE USED
      Z                 <- cbind('(constant)' = 1,method$Z)

      for( l in 1:nlambda ){
        losspoint <- vector('numeric',ngroup)
        print(M[ (j-1)*nkernel + l,,drop=FALSE])
        for( i in 1:ngroup ){
                outputs <- .svmmaj(Z[groups!=i, ], y[groups!=i],lambda=lambda[l],w=w[groups!=i],
                              increase.step=increase.step, hinge=newHinge,
                            theta=theta, check.positive=FALSE,convergence=convergence,...)
                qhat    <-  Z[groups==i,] %*% outputs$theta
                theta   <- outputs$theta
                losspoint[i] <-  sum(w[groups==i]* newHinge(qhat,y[groups==i])$loss)
                #sum(w[groups==i]*(qhat*y[groups==i]<0 ))

        }
         lMat[(j-1)*nkernel + l] <-  sum(losspoint)
      }
    }
    #=================================================
    #SAVE THE OPTIMAL VALUES
    #-------------------------------------------------
    index            <- which.min(lMat)
    #=================================================
    #RETURN THE RESULTS
    #-------------------------------------------------
    list( missclass.opt = lMat[index],
          param.opt     = M[index,],
          losses       = cbind(M,loss=lMat))
}
