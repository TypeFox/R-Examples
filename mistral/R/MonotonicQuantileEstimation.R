#-----------------------------------------------------------------------------------------------------------------#
#
#                        Quantile estimation under monotonicity constraints :
#
#-----------------------------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------------------------#
# Input :											    	    
#-----------------------------------------------------------------------------------------------------------------#
# f : failure function
# inputDimension : dimension of the input space
# inputDistribution : a list containing the name and parameters of the input distribution
# dir.monot : a vector in {-1,1}^inputDimension such that dirmonot[i] = -1 if f is decreasing and 1 if is increasing according to the ith input
# N.calls : total number of evaluations by f
# p : probability
# method : a string containing the name of the method to call
# X.input : a set of points. Useful only for ' method = "Bounds" '
# Y.input : value of f on X.input
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
# Output :	Return a list containing
#-----------------------------------------------------------------------------------------------------------------#
#  qm : vector of lower bounds of the quantile
#  qm : vector of upper bounds of the quantile
#  q.hat : an estimate of the quantile
#  Um : lower bounds of the probability obtained from the desing of experiments
#  UM : upper bounds of the probability obtained from the desing of experiments
#  XX :set of simulations
#  YY :set of values of f on XX
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#' @export
MonotonicQuantileEstimation <- function(f, inputDimension, inputDistribution, dir.monot, N.calls, p, method, X.input = NULL, Y.input = NULL){

#-------------------------------------------
# Magnitude of p
#-------------------------------------------
  k <- 1:10
  ordre.p <- which.max(floor(p*10^k) > 0)
#-----------------------------------------------------------------------------------------------------------------#
#                         Create the list containing the distributions and the parameters of each inputs
#-----------------------------------------------------------------------------------------------------------------#
  transformtionToUniformSpace <- function(inputDistribution){

    InputDist <- list()
    InputDist <- inputDistribution

    for(i in 1:inputDimension){
      nparam <- length(inputDistribution[[i]][[2]])
        for(j in 1:nparam){
          InputDist[[i]]$q <- paste("q",InputDist[[i]][[1]], sep = "");
          InputDist[[i]]$p <- paste("p",InputDist[[i]][[1]], sep = "");
          InputDist[[i]]$d <- paste("d",InputDist[[i]][[1]], sep = "");
          InputDist[[i]]$r <- paste("r",InputDist[[i]][[1]], sep = "");
        }
    }
    InputDist
  }

  InputDist <- transformtionToUniformSpace(inputDistribution)

#-----------------------------------------------------------------------------------------------------------------#
#	                           Transformation to the uniform space
#-----------------------------------------------------------------------------------------------------------------#
  G <- function(X){
    XU <- numeric()
    for(i in 1:inputDimension){
      if(dir.monot[i] == -1){X[i] <- 1 - X[i]}
        XU[i] <- do.call(InputDist[[i]]$q,c(list(X[i,drop = FALSE]), InputDist[[i]][[2]]))
      }
    return(f(XU))
  }
 
#-----------------------------------------------------------------------------------------------------------------#
#
#                         		 Test if 
#				                  "lower" ou "greater" than y.
#					             If set == 'low' (1) : return TRUE if x <= y				   
#					             If set == 'up'  (2) : return TRUE if x => y		
#		   
#-----------------------------------------------------------------------------------------------------------------#


  is.dominant <- function(x, y, inputDimension, set){

    if((set != 1)&(set != 2)){
      stop("ERROR : set must to be equal to 1 or 2.")
    } 

    dominant <- NULL;

    if( length(x) == inputDimension ){
	if(set == 2){
       if ( sum(x >= y) == inputDimension ){
          return(TRUE)
        }else{
          return(FALSE)
        } 
      }else{
        if( sum(x <= y) == inputDimension ){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
    }
	
    y.1 <- matrix(rep(y, dim(x)[1]), ncol = inputDimension, byrow = TRUE)

    if(set == 2){
      dominant <- rowSums(x >= y.1) == inputDimension
    }else{
      dominant <- rowSums(x <= y.1) == inputDimension
    }

    return(dominant)
  }


#-----------------------------------------------------------------------------------------------------------------#
#                       Return the volume of S using the Monte Carlo sample X.MC
#-----------------------------------------------------------------------------------------------------------------#
  VolumeComputation <- function(X.MC, S, set){ 

   if(set == 1){S <- 1 - S}
    if( (is.null(dim(S))) | (length(S) == inputDimension) ){
      if(set == 1){
        return(1 - prod(S))
      }
      if(set == 2){
        return(prod(S))
      }
    }

    DS <- dim(S)[1] 
    if(inputDimension == 2){
      S    <- S[order(S[,1]),]
      res  <- diag(outer(S[,1], c(0,S[1:(DS - 1), 1]), "-"))
      res1 <- res%*%S[,2]
      res1 <- ifelse(set == 1, 1 - res1, res1)
      return(res1) 
    }


    MC.VOL <- 0
    RES.VOL <- 0
    MM <- dim(X.MC)[1]

    for(i in 1:(DS-1)){
      u  <- apply(S, MARGIN = 1, prod)
      u.max <- which.max(u)
      uu    <- S[u.max, ]
      S     <- S[-u.max, ]
      ss.1 <- is.dominant(X.MC, uu, inputDimension, 1)
      ss.2 <- is.dominant(X.MC, uu, inputDimension, 2)
      RES.VOL <- RES.VOL + sum(ss.1)
      X.MC <- X.MC[which((ss.1 == 0)&(ss.2 == 0) ) , ]
    }

    uu <- S
    if(set == 1){
      tt <- is.dominant(X.MC, uu, inputDimension, 2)
      RES.VOL <- RES.VOL + sum(tt)
      RES.VOL <- 1 - RES.VOL/MM
    }else{
      tt <- is.dominant(X.MC, uu, inputDimension, 1)
      RES.VOL <- RES.VOL + sum(tt)
      RES.VOL <- RES.VOL/MM
    }
    return(RES.VOL) 

  }

  VolumeComputation.quant <- function(X.MC, S, set, p){ 


    RES.VOL <- 0
    MM <- dim(X.MC)[1]

    if(set == 1){
      temp <- apply(X.MC, MARGIN = 1, function(w){prod(w)} )
      temp.vol <- (temp >= p)
    }
    if(set == 2){
      temp <- apply(X.MC, MARGIN = 1, function(w){prod(1-w)} )
      temp.vol <- ( temp >= (1-p) )
    }

    RES.VOL <- sum(temp.vol)
    X.MC <- X.MC[which(!temp.vol), ]


    if( length(S) == inputDimension){
      if(set == 1){
        s.temp <- is.dominant(X.MC, S, inputDimension, 2 )
        RES.VOL <- RES.VOL + sum(s.temp)
        return(1 - RES.VOL/MM)
      }
      if(set == 2){
        s.temp <- is.dominant(X.MC, S, inputDimension, 1 )
        return(RES.VOL/MM)
      }
    }

    DS <- dim(S)[1] 
    if(set == 1){
      S    <- 1 - S
      X.MC <- 1 - X.MC
    }

    for(i in 1:(DS-1)){
      u  <- apply(S, MARGIN = 1, prod)
      u.max <- which.max(u)
      uu    <- S[u.max, ]
      S     <- S[-u.max, ]
      ss.1 <- is.dominant(X.MC, uu, inputDimension, 1)
      ss.2 <- is.dominant(X.MC, uu, inputDimension, 2)
      RES.VOL <- RES.VOL + sum(ss.1)


      X.MC <- X.MC[which((ss.1 == 0)&(ss.2 == 0) ), ]
    }

    uu <- S
    if(set == 1){
      tt <- is.dominant(X.MC, uu, inputDimension, 2)
      RES.VOL <- RES.VOL + sum(tt)
      RES.VOL <- 1 - RES.VOL/MM
    }else{
      tt <- is.dominant(X.MC, uu, inputDimension, 1)
      RES.VOL <- RES.VOL + sum(tt)
      RES.VOL <- RES.VOL/MM
    }
    return(RES.VOL) 

  }

#-----------------------------------------------------------------------------------------------------------------
#	Return the frontier of a set
#-----------------------------------------------------------------------------------------------------------------

  Frontier <- function(S, set){
    if(is.null(dim(S)) |( dim(S)[1] == 1)){
      return(S)
    }
    R <- NULL
    if(set == 1){
      S <- 1 - S
    }
    while( !is.null(dim(S)) ){
      if(dim(S)[1] == 0){
        if(set == 1){
          return(1 - R)
        }else{
          return(R)
        }
      }
      aa  <- apply(S, MARGIN = 1, prod)
      temp <- S[which.max(aa), ]
      R <- rbind(R, temp)
      S <- S[-which.max(aa), ]
      ss <- is.dominant(S, temp, inputDimension, 1)
      if(!is.null(dim(S))){
        S <- S[which(ss == FALSE),]
      }else{
        S <- matrix(S, ncol= inputDimension)
        S <- S[which(ss == FALSE), ]
        R <- rbind(R, S)
        if(set == 1){
          return(1 - R)
        }else{
          return(R)
        }
      }
    }
   
    if(set == 1){
      return(1 - R)
    }else{
      return(R)
    }
  }

  Frontier.quant <- function(S, Y.S, set){
    if(is.null(dim(S)) |( dim(S)[1] == 1)){
      return(S)
    }
    R <- NULL
    R.Y <- NULL

    if(set == 1){
      S <- 1 - S
    }
    while( !is.null(dim(S)) ){
      if(dim(S)[1] == 0){
        if(set == 1){
          return( list(1 - R, R.Y))
        }else{
          return(list(R, R.Y))
        }
      }
      aa  <- apply(S, MARGIN = 1, prod)
      temp <- S[which.max(aa), ]

      R <- rbind(R, temp)
      R.Y <- c(R.Y, Y.S[which.max(aa)])

      S <- S[-which.max(aa), ]
      Y.S <- Y.S[-which.max(aa) ]
      ss <- is.dominant(S, temp, inputDimension, 1)
      if(!is.null(dim(S))){
        S <- S[which(ss == FALSE),]
        Y.S <- Y.S[which(ss == FALSE)]
      }else{
        S <- matrix(S, ncol= inputDimension)
        S <- S[which(ss == FALSE), ]
        Y.S <- Y.S[which(ss == FALSE)]

        R <- rbind(R, S)
        R.Y <- c(R.Y, Y.S)
        if(set == 1){
          return(list(1 - R, R.Y))
        }else{
          return(list(R, R.Y))
        }
      }
    }
   
    if(set == 1){
      return(list(1 - R, R.Y))
    }else{
      return(list(R,R.Y))
    }
  }

  Frontier.quantile <- function(S, set){
    if(is.null(dim(S)) |( dim(S)[1] == 1)){
      return(S)
    }
    R <- NULL
    V <- NULL
    if(set == 1){
      S <- 1 - S
    }
    while( !is.null(dim(S)) ){
      if(dim(S)[1] == 0){
        if(set == 1){
          return( list( 1 - R , NULL))
        }else{
          return(list( R, NULL ))
        }
      }
      aa  <- apply(S, MARGIN = 1, prod)
      temp <- S[which.max(aa), ]
      R <- rbind(R, temp)
      S <- S[-which.max(aa), ]
      ss <- is.dominant(S, temp, inputDimension, -1)
      if(!is.null(dim(S))){
        V <- rbind(V, S[which(ss == TRUE), ])
        S <- S[which(ss == FALSE),]
      }else{
        S <- matrix(S, ncol= inputDimension)
        V <- rbind(V, S[which(ss == TRUE), ])
        S <- S[which(ss == FALSE), ]
        R <- rbind(R, S)
        if(set == 1){
          return( list( 1 - R, 1- V) )
        }else{
          return( list(R, V) )
        }
      }
    }
   
    if(set == 1){
      return( list( 1 - R, 1- V) )
    }else{
      return( list(R, V) )
    }
  }
#-----------------------------------------------------------------------------------------------------------------
#	   Convert an integer into a binary
#-----------------------------------------------------------------------------------------------------------------
    as.binary <- function (x) { 
      base <- 2;
      r <- numeric(inputDimension)
      for (i in inputDimension:1){ 
        r[i] <- x%%base 
	  x <- x%/%base
      } 
      return(r) 
   }

#-----------------------------------------------------------------------------------------------------------------
#		                             One point is simulated around x
#-----------------------------------------------------------------------------------------------------------------
    SIM <- function(x, W){

      B <- 0
      # The space in split into 2^inputDimension - 1 subsets then the volume of these subsets is computed
      B <- apply( matrix(W, ncol = 1), 
                  MARGIN = 1, 
                  function(v){
                    Z <- as.binary(v)
                    v <- 0
                    u <- 0
                    for(j in 1:inputDimension){
                      u[j] <- ifelse(Z[j] == 0, 1 - x[j], x[j])
                    }     
                    return(prod(u))  
                  }
                 )

      B <- cumsum(B)
      U <- runif(1, 0, max(B))

      # One of the subsets is chosen randomly according to its volume
      pos <- ifelse(U < B[1], 1, which.max(B[B <= U]) + 1)

      Z <- as.binary(W[pos])
      A <- 0
      
      # One point is simulated in the chosen subset
      for(i in 1:inputDimension){
        A[i] = ifelse(Z[i] == 0, runif(1, x[i], 1), runif(1, 0, x[i]))  
      }
      return(A)      
    }

#-----------------------------------------------------------------------------------------------------------------
#                   Simulation of CP points in the non-dominated set
#-----------------------------------------------------------------------------------------------------------------


    Sim.non.dominated.space <- function(CP, Z.safe, Z.fail, W){
      CP1 <- 0;
      Y   <- NULL
      Y.temp  <- apply(1 - Z.safe, MARGIN = 1, prod)
      Y.temp1 <- Z.safe[which.max(Y.temp),]
      while(CP1 < CP){
        Y.temp2 <- SIM(Y.temp1, W)
        tts1 <- is.dominant(Z.safe, Y.temp2, inputDimension, 1)
  	ttf1 <- is.dominant(Z.fail, Y.temp2, inputDimension, 2)
    	if( (sum(tts1) == 0 ) & ( sum(ttf1) == 0) ){
    	  Y <- rbind(Y,Y.temp2)
    	  CP1 <- CP1 + 1
    	}
      }
      return(Y)
    }

    Sim.non.dominated.space.quant <- function(CP, Z.safe, Z.fail, W){
      CP1 <- 0;
      Y   <- NULL
      Y.temp  <- apply(1 - Z.safe, MARGIN = 1, prod)
      Y.temp1 <- Z.safe[which.max(Y.temp),]
      while(CP1 < CP){
        Y.temp2 <- SIM(Y.temp1, W)
	if((prod(Y.temp2) <= p) & ( prod(1-Y.temp2) <= (1-p)) ){ #A condition to diminish the number of test
          tts1 <- is.dominant(Z.safe, Y.temp2, inputDimension, 1)
  	  ttf1 <- is.dominant(Z.fail, Y.temp2, inputDimension, 2)
    	  if( (sum(tts1) == 0 ) & ( sum(ttf1) == 0) ){
    	    Y <- rbind(Y,Y.temp2)
    	    CP1 <- CP1 + 1
    	  }
	} 
      }
      return(Y)
    }


#####################################################################
# Empirical quantile estimation
#####################################################################

  Method.quantil.empirical <- function(N.calls, H){
    X <- NULL
    Y <- res <- numeric(N.calls)
    for(i in 1:N.calls){
      X.temp <- runif(inputDimension)
      X <- rbind(X, X.temp)
      Y.temp <- H(X.temp)
      Y[i]   <- Y.temp  
      res[i] <- quantile(Y,p)
    }
    return( list(res, X, Y))
  }

#####################################################################
#####################################################################

  Delimite.frontier <- function(S, Y, set, pp, UU){
    if(set == 1){
      S <- 1 - S
    }

    vol.temp <- apply(S, MARGIN = 1, prod)
    x.temp <- S[ which.min(vol.temp), ]
    S <- S[-which.min(vol.temp), ]

    lower.p <- prod(x.temp)
    R.S <- x.temp
    Y.temp <- Y
    Y.lower.bound <- Y.temp[which.min(vol.temp)]
    Y.temp <- Y.temp[-which.min(vol.temp)]

    Y.rep <- Y.lower.bound

    if( length(S) == inputDimension){
      S <- matrix(S, ncol = inputDimension)
    }

    while(lower.p < pp){
      Y.lower.bound_old <- Y.lower.bound


      vol.temp <- apply(S, MARGIN = 1, function(u){v.temp <- rbind(x.temp,u); 
                                                        v.temp <- Frontier(v.temp, 2);
                                                        if(length(v.temp) == inputDimension){return(prod(v.temp))}
                                                        return( VolumeComputation(UU, v.temp, 2) ) })

      x.temp <- rbind(x.temp, S[ which.min(vol.temp), ])
      S <- S[-which.min(vol.temp), ]

      if( length(S) == inputDimension ){
        S <- matrix(S, ncol = inputDimension)
      }

      Y.rep <- c(Y.rep, Y.temp[ which.min(vol.temp) ])

      V.temp <- Frontier.quant(x.temp, Y.rep, 2)
      x.temp <- V.temp[[1]]
      Y.rep <- V.temp[[2]]
      Y.lower.bound <- Y.temp[ which.min(vol.temp) ]
      Y.temp <- Y.temp[-which.min(vol.temp)]
      lower.p <- min(vol.temp)

      if( dim(S)[1] == 0 ){
        if(set == 1){
          return(min(Y.rep))
        }else{
          return(max(Y.rep))
        }
      }

    }

    if(set == 1){
      return(min(Y.rep))
    }else{
      return(max(Y.rep))
    }

  }

#####################################################################
#####################################################################

  Delimite.frontier.Quant <- function(S, Y, set, pp, UU){
    if(set == 1){
      S <- 1 - S
    }
    pos <- 1:dim(Y)[1]

    vol.temp <- apply(S, MARGIN = 1, prod)
    x.temp <- S[ which.min(vol.temp), ]
    S <- S[-which.min(vol.temp), ]    

    lower.p <- prod(x.temp)
    R.S <- x.temp
    Y.temp <- Y
    Y.lower.bound <- Y.temp[which.min(vol.temp)]
    Y.temp <- Y.temp[-which.min(vol.temp)]
    pos <- pos[-which.min(vol.temp)]

    Y.rep <- Y.lower.bound

    if( length(S) == inputDimension){
      S <- matrix(S, ncol = inputDimension)
    }

    while(lower.p < pp){
      Y.lower.bound_old <- Y.lower.bound
      vol.temp <- apply(S, MARGIN = 1, function(u){v.temp <- rbind(x.temp,u); 
                                                        v.temp <- Frontier(v.temp, 2);
                                                        if(length(v.temp) == inputDimension){return(prod(v.temp))}
                                                        return( VolumeComputation(UU, v.temp, 2) ) })

      x.temp <- rbind(x.temp, S[ which.min(vol.temp), ])
      S <- S[-which.min(vol.temp), ]
      if( length(S) == inputDimension ){
        S <- matrix(S, ncol = inputDimension)
      }

      Y.rep <- c(Y.rep, Y.temp[ which.min(vol.temp) ])

      V.temp <- Frontier.quant(x.temp, Y.rep, 2)
      x.temp <- V.temp[[1]]
      Y.rep <- V.temp[[2]]
      Y.lower.bound <- Y.temp[ which.min(vol.temp) ]
      Y.temp <- Y.temp[-which.min(vol.temp)]
      pos <- pos[-which.min(vol.temp)]
      lower.p <- min(vol.temp)

      if( dim(S)[1] == 0 ){
        if(set == 1){
          return(min(Y.rep))
        }else{
          return(max(Y.rep))
        }
      }

    }

    if(set == 1){
      return( pos[which.min(Y.rep) ])
    }else{
      return( pos[which.max(Y.rep) ])
    }

  }

#####################################################################
#####################################################################

  Method.quantil.MC <- function(N.calls, H){


    if(inputDimension > 2){
      UU <- runif(inputDimension*10^(ordre.p + 2))
      UU <- matrix(UU, ncol=inputDimension, byrow = TRUE)
      D.UU <- dim(UU)[1]
    }

    if(inputDimension == 2){UU <- 0; D.UU <- 0}

    X <- NULL
    Y <- 0
    Res <- 0

    for(i in 1:N.calls){
      X.temp <- runif(inputDimension)
      X <- rbind(X, X.temp)
      Y.temp <- H(X.temp)
      Y[i]   <- Y.temp
      Res[i] <- sort(Y)[1 + floor(i*p)]
    }

    #Search of an upper bound for the quantile
    XS <- Frontier(X, -1)
    vS <- VolumeComputation(UU, XS, 2) 
    if(vS <  p){
      upper.q <- Inf
    }else{
      upper.q <- Delimite.frontier(X, Y, 2,  p, UU)
    }


    #Search of an lower bound for the quantile
    XF <- Frontier(X, 1)
    vF <- 1-VolumeComputation(UU, XF, 1) 
    if(vF < (1 - p)){
      lower.q <- -Inf
    }else{
      lower.q <- Delimite.frontier(X, Y, 1, 1 - p, UU)
    }

   return(list(lower.q, upper.q, Res, X,Y ))
  }

#####################################################################
#####################################################################

  Method.bounds <- function(X.input, Y.input){
    X <- X.input
    Y <- Y.input
    lower.q <- -Inf
    upper.q <- Inf
    init.up <- apply(X, MARGIN =  1, function(w){prod(w)})
    X.up.temp <- X[ which(init.up >= p), ]
    if(length(X.up.temp) > 0){              #If there is somes points in the safety set
      if( dim(X.up.temp)[1] == dim(X)[1] ){ #If there is somes points in the failure set
        return(c(-Inf , min(Y)) )
      }else{ #The points in the safety set are retrieve
        upper.q <- min(Y[which(init.up >= p) ])
        X <- X[which(init.up < p), ]
        Y <- Y[which(init.up < p) ]
        if(sum(Y >= upper.q) > 0){
          pos.temp <- which(Y >= upper.q)
          if(length(X) == inputDimension){
            X <- matrix(X, ncol = inputDimension)
          }
          X <- X[-pos.temp, ]
          Y <- Y[-pos.temp]
        }
        if(length(X) == 0){
          return(c(lower.q,upper.q))
        }
      }
    }
    if(length(X) == inputDimension){
      X <- matrix(X, ncol = inputDimension)
    }

    init.low <- apply(X, MARGIN =  1, function(w){prod(1-w)})
    X.low.temp <- X[which(init.low >= (1 - p) ), ]

    if(length(X.low.temp) == inputDimension){
      X.low.temp <- matrix(X.low.temp, ncol = inputDimension)
    }

    if(dim(X.low.temp)[1] == dim(X)[1]){ #If all points are in the failure set
      return(c(max(Y) , upper.q))
    }
    if((length(X.low.temp) > 0)&(dim(X.low.temp)[1] < dim(X)[1])){#If it remains some points in the failure set
      lower.q <- max(Y[which(init.low >= (1-p)) ])
      X <- X[which(init.low < (1-p)), ]
      Y <- Y[which(init.low < (1-p)) ]
      if(sum(Y <= lower.q) > 0){
        pos.temp <- which(Y <= lower.q)
        if(length(X) == inputDimension){
          X <- matrix(X, ncol = inputDimension)
        }
        X <- X[-pos.temp, ]
        Y <- Y[-pos.temp]
      }
      if(length(X) == 0){
        return(c(lower.q,upper.q))
      }
    }    
    
    if(length(X) == inputDimension){
      return(c(lower.q,upper.q))
    }

    UU <- runif(inputDimension*10^(ordre.p + 2))
    UU <- matrix(UU, ncol=inputDimension, byrow = TRUE)
    D.UU <- dim(UU)[1]

    XS <- Frontier(X, -1)
    vS <- VolumeComputation(UU, XS, 2) 
    if(vS >=  p){
      upper.q <- Delimite.frontier(X, Y, -1,  p, UU)
    }

    XF <- Frontier(X, 1)
    vF <- 1-VolumeComputation(UU, XF, 1) 
    if(vF >= (1 - p)){
      lower.q <- Delimite.frontier(X, Y, 1, 1 - p, UU)
    }
 
    return(c(lower.q, upper.q))
  }

#####################################################################
#####################################################################

  Method.quantil.IS <- function(N.calls, H){

########################################################
#            Empirical cumulation distribution function
########################################################
    Empirical.cdf <- function(yy, YY, Bounds){
      res <- 0
      for(i in 1:length(yy)){
          res[i] <- mean(Bounds[ ,1] + (Bounds[ ,2]-Bounds[ ,1])*(YY <= yy[i]))
      }
      return(res)
    }

    UU <- runif(inputDimension*10^(ordre.p + 2))
    UU <- matrix(UU, ncol=inputDimension, byrow = TRUE)
    D.UU <- dim(UU)[1]

    VAR <- 0
    VAR.2 <- 0
    IC.inf <- 0
    IC.sup <- 0

    IC.inf.2 <- 0
    IC.sup.2 <- 0

    diff.y <- 0

    alpha <- 0.025
    eps   <- 1e-7

    xM <- rep(p^(1/inputDimension), inputDimension)
    UM <- 1 - prod(1 - xM )
    qM <- H(xM)

    xm <- rep( 1 - (1-p)^(1/inputDimension), inputDimension)
    Um <- prod(xm)
    qm <- H( xm)

    Z.safe <- matrix( xM, ncol=inputDimension)
    Z.fail <- matrix(xm , ncol=inputDimension)

    W  <- 1:(2^(inputDimension) - 1)
    XX <- NULL

    X <-  Sim.non.dominated.space.quant(1, Z.safe, Z.fail, W)

    H.X <- H(X)

    Y <- YY <- H.X
  
    Z <- Y.Z <- NULL
    XX <- rbind(XX, X)

    y <- seq(from = qm, to = qM, length = 1000)
    p.hat.temp <- Um + (UM - Um)*(H.X <= y)

    q.hat <- 0
    q.hat.3 <- (qm + qM)/2

    q.hat[2] <- y[which.max(p.hat.temp > p)]
    q.hat.3[2] <- q.hat.3[1] + (q.hat.3[1] - (H.X <= q.hat.3[1]))/2

    qm[2] <- qM[2] <- 0

    A.test <- (prod(X) > p)
    B.test <- (  prod(1-X) > 1- p)

    if(A.test){
      Z.safe <- rbind(Z.safe, X)
      Z.safe <- Frontier(Z.safe, 1)
      qM[2] <- ifelse(H.X <= qM[1], H.X, qM[1] )
      qm[2] <- qm[1]
      UM[2] <- ifelse(length(Z.safe) == inputDimension, 1 - prod(1 - Z.safe), VolumeComputation.quant(UU, Z.safe, 1, p) )
      Um[2] <- Um[1]
    }

    if( B.test ){
      Z.fail <- rbind(Z.fail, X)
      Z.fail <- Frontier(Z.fail, 2)
      qm[2] <- ifelse(H.X >= qm[1], H.X, qm[1] )
      qM[2] <- qM[1]
      Um[2] <- ifelse(length(Z.fail) == inputDimension, prod(Z.fail), VolumeComputation.quant(UU,Z.fail, 2) )
      UM[2] <- UM[1]
    }

    if(!B.test & !A.test){
      Z <- rbind(Z,X)
      Y.Z <- c(Y.Z, H.X)

      qm[2] <- qm[1]
      qM[2] <- qM[1]
      Um[2] <- Um[1]
      UM[2] <- UM[1]
    }

    j <- 3
    while(j <= N.calls ){
      cat(paste(j, " -- "))

      X <-  Sim.non.dominated.space.quant(1, Z.safe, Z.fail, W)
      H.X <- H(X)
      Y <- c(Y, H.X)
      YY <- c(YY, H.X)
      XX <- rbind(XX, X)

      if(prod(X) >= p){
        qM[j]  <- ifelse( H.X < qM[j-1], H.X, qM[j-1] )
        qm[j]  <- qm[j-1]
        Z.safe <- rbind(Z.safe, X)
        Z.safe <- Frontier(Z.safe, 1)
        UM[j]  <- VolumeComputation.quant(UU, Z.safe, 1)
        Um[j]  <- Um[j-1]
      }

      if(prod(1 - X) >= 1 - p){
        qm[j] <- ifelse( H.X > qm[j-1], H.X, qm[j-1] )
        qM[j] <- qM[j-1]
        Z.fail <- rbind(Z.fail, X)
        Z.fail <- Frontier(Z.fail, 2)
        Um[j] <- VolumeComputation.quant(UU,Z.fail, 2, p)
        UM[j] <- UM[j-1]
      }

      if( (prod(X) < p) & (prod(1 - X) < 1 - p) ){

        if(H.X >= qM[j-1] ){
          Z.safe <- rbind(Z.safe, X)
          Z.safe <- Frontier(Z.safe, 1)
          if(length(Z.safe) == inputDimension){
            Z.safe <- matrix(Z.safe, ncol = inputDimension)
          }
          UM[j] <- VolumeComputation.quant(UU, Z.safe, 1, p)
          Um[j] <- Um[j-1]
          qM[j] <- qM[j-1]
          qm[j] <- qm[j-1]
        }#End if H.X >= qM
        if(H.X <= qm[j-1]){
          Z.fail <- rbind(Z.fail, X)
          Z.fail <- Frontier(Z.fail, 2)
          if(length(Z.fail) == inputDimension){
            Z.fail <- matrix(Z.fail, ncol = inputDimension)
          }
          Um[j] <- VolumeComputation.quant(UU, Z.fail, 2, p)
          UM[j] <- UM[j-1]
          qm[j] <- qm[j-1]
          qM[j] <- qM[j-1]
        }#End if H.X <= qm

        #If there is non information about X
        if(   (H.X > qm[j-1])&(H.X < qM[j-1])){
          UM[j] <- UM[j-1]
          Z <- rbind(Z, X)    
          Y.Z <- c(Y.Z, H.X)
          ZZ.temp <- Frontier.quant(Z, Y.Z, 1)
          Z.temp  <- ZZ.temp[[1]]
          YZ.temp <- ZZ.temp[[2]]
          vol.safe.temp <- 1 - VolumeComputation.quant(UU, Z.temp, 1, p)
          if( vol.safe.temp > 1 - p ){
            #Some information can be deduced from X

            if(length(Z) == inputDimension){
              Z <- matrix(Z, ncol = inputDimension)
              lower.q <- Y.Z
            }else{
              lower.q <- Delimite.frontier(Z, Y.Z, set = 1, pp = 1 - p, UU)
            }
            qm[j] <- ifelse( lower.q > qm[j-1], lower.q, qm[j-1] )
            pos.Y <- which(Y.Z == lower.q)
            uu <- Z[pos.Y, ]
            ff <- is.dominant(Z, uu, inputDimension, 1)
            Z  <- Z[which(ff == FALSE), ]
            Y.Z  <- Y.Z[which(ff == FALSE) ]
            Z.fail <- rbind(Z.fail, uu)
            Um[j] <- VolumeComputation.quant(UU, Z.fail, 2, p)
          }else{
            qm[j] <- qm[j-1]
            Um[j] <- Um[j-1]
          }


          if(length(Z) == 0){
            qM[j] <- qM[j-1]
            UM[j] <- UM[j-1]
          }else{
            Z.temp <- Frontier(Z, 2)
            if(length(Z.temp) == inputDimension){
              vol.fail.temp <-  prod(Z.temp)
            }else{
              vol.fail.temp <- VolumeComputation.quant(UU, Z.temp, 2, p)
            } 
            if(vol.fail.temp >  p ){
              if(length(Z) ==inputDimension ){
                Z <- matrix(Z, ncol = inputDimension)
                upper.q <- Y.Z
              }else{
                upper.q <- Delimite.frontier(Z, Y.Z, set =  - 1, pp = p, UU)
              }
              qM[j] <- ifelse( upper.q < qM[j-1], upper.q, qM[j-1] )
              pos.Y <- which(Y.Z == upper.q)
              uu <- Z[pos.Y, ]
              ss <- is.dominant(Z, uu, inputDimension, 1)
              Z <- Z[which(ss == FALSE), ]
              Y.Z <- Y.Z[which(ss == FALSE) ]
              Z.safe <- rbind(Z.safe, uu)
              UM[j] <- VolumeComputation.quant(UU, Z.safe, 1, p)  
            }else{
              qM[j] <- qM[j-1]
              UM[j] <- UM[j-1]
            }#End if vol.fail.temp > p
          }
        }#End if (H.X > qm[j-1])&(H.X < qM[j-1])
      }#End (prod(X) < p) & (prod(1 - X) < 1 - p)

      YY.temp <- sort(YY)
      pTemp <- Empirical.cdf(YY.temp, YY, cbind(Um, UM)[-j, ])
      q.hat.temp <- YY.temp[which.max(pTemp > p)]
      if(q.hat.temp > qM[j]){ 
        q.hat[j] <- qM[j]
      }
      if(q.hat.temp < qm[j]){ 
        q.hat[j] <- qm[j]
      }

      if( (q.hat.temp < qM[j]) & (q.hat.temp > qm[j])){
        q.hat[j] <- q.hat.temp
      }
      if(is.na(q.hat[j])){
        q.hat[j] <- qM[j]
      }
      j <- j + 1
    }
    return(list(cbind(qm, qM, q.hat, Um,UM),XX, YY))
  }

#----------------------------------------------------------------------------------------------------------------- 
#-----------------------------------------------------------------------------------------------------------------  
   if(method == "MonteCarlo"){
     RESULT <- Method.quantil.MC(N.calls, G)
   }

   if(method == "MonteCarloWB"){
     RESULT <- Method.quantil.empirical(N.calls, G)
   }

   if(method == "MonteCarloIS"){
     RESULT <- Method.quantil.IS(N.calls, G)
   }

   if(method == "Bounds"){
     RESULT <- Method.bounds(X.input, Y.input)
   }

  return(RESULT)

}

