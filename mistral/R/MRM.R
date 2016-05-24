## -----------------------------------------------------------------------------
## Fonction MRM
## -----------------------------------------------------------------------------
#' @export
MRM <- function(f, inputDimension, inputDistribution, dir.monot, N.calls, Method, silent = FALSE){


  transformtionToInputSpace <- function(inputDistribution){

    InputDist <- list()
    InputDist <- inputDistribution

    for(i in 1:inputDimension){
      nparam <- length(inputDistribution[[i]][[2]])
        for(j in 1:nparam){
          InputDist[[i]]$q <- paste("q", InputDist[[i]][[1]], sep = "");
          InputDist[[i]]$p <- paste("p", InputDist[[i]][[1]], sep = "");
          InputDist[[i]]$d <- paste("d", InputDist[[i]][[1]], sep = "");
          InputDist[[i]]$r <- paste("r", InputDist[[i]][[1]], sep = "");
        }
    }
    InputDist
  }

  InputDist <- transformtionToInputSpace(inputDistribution)

#-----------------------------------------------------------------------------------------------------------------#
#	                           Transformation in the uniform space  
#-----------------------------------------------------------------------------------------------------------------#
  G <- function(X){
    XU <- numeric()
    for(i in 1:inputDimension){
      if(dir.monot[i] == -1){X[i] <- 1 - X[i]}
        XU[i] <- do.call(InputDist[[i]]$q,c(list(X[i, drop = FALSE]), InputDist[[i]][[2]]))
      }
    return(f(XU))
  }

#-----------------------------------------------------------------------------------------------------------------#
#                Initialisation of the method : a dichotomie on the Uniforme space
#-----------------------------------------------------------------------------------------------------------------#

  Intersect <- function(inputDimension, FUNC){

    a     <- 2
    k     <- 2
    res   <- list()
    u.new <- 0
    temp  <- 0
    u.dep <- list()
    out   <- list()
    comp  <- 2

    u.dep[[1]]  <- rep(1/2, inputDimension)
    temp        <- FUNC(u.dep[[1]])
    u.dep[[2]]  <- sign(temp)
    cp 	        <- 1 						# compteur d'appel ? G
    u.other     <- u.dep
    u.new       <- u.dep[[1]]
    LIST        <- list()
    LIST[[1]]   <- u.dep
    list.set    <- LIST[[1]][2]


    if(temp > 0){
      u.new <- u.dep[[1]] - 1/(a^k)
    }else{
      u.new <- u.dep[[1]] + 1/(a^k)
    }

    eps <- ( u.new - u.dep[[1]] )%*%( u.new - u.dep[[1]] )

    if( ( u.other[[2]] != sign(temp)) & (eps > 1e-7) ){
      u.other <- list( u.dep[[1]], sign(temp) )
    }

    k          <- k + 1
    sign.0     <- sign(temp)
    sign.other <- - sign.0

    while(sign(temp) != sign.other){         
 
      u.dep[[1]] <- u.new 
      temp       <- FUNC(u.dep[[1]])
      u.dep[[2]] <- sign(temp)
      cp 	 <- cp + 1

      if(temp > 0){
        u.new <- u.dep[[1]] - 1/(a^k)
      }else{
        u.new <- u.dep[[1]] + 1/(a^k)
      }

      eps <- ( u.new - u.dep[[1]] )%*%( u.new - u.dep[[1]] )
      k              <- k + 1
      LIST[[comp]]   <- u.dep
      list.set[comp] <- LIST[[comp]][[2]]
      comp           <- comp + 1
    }

    return(LIST)
    list.set <- as.numeric(list.set)

    if( abs(sum(list.set)) == length(LIST)){

      res[[1]] <- LIST[[length(LIST)]][[1]]
      res[[2]] <- LIST[[length(LIST)]][[2]]
      out <- list(res, cp)
      return(out)

    }else{

      u.dep[[1]] <- LIST[[max(which(list.set == -1))]][[1]]
      u.dep[[2]] <- LIST[[max(which(list.set == -1))]][[2]]

      u.other[[1]] <- LIST[[max(which(list.set == 1))]][[1]]
      u.other[[2]] <- LIST[[max(which(list.set == 1))]][[2]]

      res[[1]] <- rbind(u.dep[[1]], u.other[[1]])
      res[[2]] <- c(u.dep[[2]], u.other[[2]])

      out <- list(res, cp)

      return(out)
    }

  }

#-----------------------------------------------------------------------------------------------------------------#
#
#                         		 Test if the points of the set "x" are
#				                  "smaller" or "greater"  than "y".
#					             If set == 1 : Return TRUE if x[i, ] <= y				   
#					             If set == 2 : Return TRUE if x[i, ] => y		
#		   
#-----------------------------------------------------------------------------------------------------------------#

  is.dominant <- function(x, y, inputDimension, set){

    dominant <- NULL;

    if( is.null(dim(x)) ){
	if(set == -1){
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
	
    y.1 <- NULL
    Y.2 <- NULL
    y.1 <- rep(y,dim(x)[1])
    y.2 <- matrix(y.1, ncol = inputDimension, byrow = TRUE)

    if(set == -1){
      dominant <- apply(x >= y.2, 1, sum) == inputDimension
    }else{
      dominant <- apply(x <= y.2, 1, sum) == inputDimension
    }

    return(dominant)
  }



#-----------------------------------------------------------------------------------------------------------------#
#		              Function which compute the exact bounds
#-----------------------------------------------------------------------------------------------------------------#
  
  Volume.bounds <- function(S, set){ 

   if(set == 1){
     S <- 1 - S
   }
 
 
    if(is.null(dim(S))){
      if(set == 1){
        return(1 - prod(S))
      }
      if(set == -1){
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

   RES.VOL <- dominated_hypervolume(1 - t(S), rep(1, inputDimension))

    if(set == 1){
      RES.VOL <- 1 - RES.VOL
    }

    return(RES.VOL) 

  }

#-----------------------------------------------------------------------------------------------------------------
#				                     Gives the frontier of a set
#-----------------------------------------------------------------------------------------------------------------

  Frontier <- function(S, set){
    if(is.null(dim(S)) |( dim(S)[1] == 1)){
      return(S)
    }
    R <- NULL
    if(set == 1){
      S <- 1 - S
    }
    while(!is.null(dim(S))){
      aa  <- apply(S, MARGIN = 1, prod)
      temp <- S[which.max(aa), ]
      R <- rbind(R, temp)
      S <- S[-which.max(aa), ]
      ss <- is.dominant(S, temp, inputDimension, -1)
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


#-----------------------------------------------------------------------------------------------------------------
#
#
#
#     				Monte Carlo method under monotonicity constraints
#
#
#-----------------------------------------------------------------------------------------------------------------
  monteCarloMonotone <- function(N.calls){ 


    NN    <- 0        #Number of call to f
    N.tot <- 0
    res   <- NULL

    Z.safe <- NULL
    Z.fail <- NULL

    X <- NULL
    Y <- NULL

    is.Call <- 0    
    while(NN < N.calls){
      if(silent == FALSE){
        if(N.tot%%100 == 0){print(NN);flush.console();}
      }

      U <- runif(inputDimension)
      X <- rbind(X, U)
      N.tot <- N.tot + 1

      if( is.null(Z.safe)& is.null(Z.fail) ){
        t.u <- G(U)
        NN  <- NN + 1

        is.Call[NN] <- N.tot
      }

      if(is.null(Z.safe)&( !is.null(Z.fail)) ){
        ttf <- is.dominant(Z.fail, U, inputDimension, set = -1)
        if( sum(ttf) == 0 ){
          t.u <- G(U)
          NN  <- NN + 1


          is.Call[NN] <- N.tot          
        }else{
          t.u <- -1
        }
      }
      
      if(!is.null(Z.safe) & is.null(Z.fail) ){       
        tts <- is.dominant(Z.safe, U, inputDimension, set = 1)
        if(sum(tts) == 0){
          t.u <- G(U)
          NN  <- NN + 1
          is.Call[NN] <- N.tot          
        }else{
          t.u <- 1
        }
      }
    
      if((!is.null(Z.safe)) &( !is.null(Z.fail)) ){      
        ttf <- is.dominant(Z.fail, U, inputDimension, set = -1)       
        tts <- is.dominant(Z.safe, U, inputDimension, set = 1)
        if( (sum(tts) == 0) & (sum(ttf) == 0) ){
          t.u <- G(U)
          NN  <- NN + 1
          is.Call[NN] <- N.tot
        }
        if( (sum(tts) == 0)& (sum(ttf) != 0) ){
          t.u <- -1
        }
        if( (sum(tts) != 0)& (sum(ttf) == 0) ){
          t.u  <- 1
        }
        
      }

      Y <- c(Y, t.u)

      if(t.u <= 0){
        Z.fail <- rbind(Z.fail, U)
        res    <- c(res, 1)
      }else{        
        Z.safe <- rbind(Z.safe, U)
        res    <- c(res, 0)
      }

    }

    I <- 1:N.tot

    alpha <- 0.05
 
    cum.res <- cumsum(res)

    estimation_MC <- cum.res/I                                        #Monte Carlo Estimator
    Var_MC        <- (estimation_MC)*(1 - estimation_MC)/I            #Variance of the estimator
    IC.inf        <- estimation_MC - qnorm(1 - alpha/2)*sqrt(Var_MC)  #Confidence Interval
    IC.sup        <- estimation_MC + qnorm(1 - alpha/2)*sqrt(Var_MC)  #Confidence Interval
    CV_MC         <- 100*sqrt(Var_MC)/estimation_MC

    if(is.null(Z.fail)){
      Um <- 0
    }else{
      ZF <- Frontier(Z.fail, -1)        
      Um <- Volume.bounds(ZF, -1)    
    }

    if(is.null(Z.safe)){
      UM <-1
    }else{ 
      ZS <- Frontier(Z.safe, 1) 
      UM <- Volume.bounds(ZS, 1)  
    }
    
    return(list(cbind(IC.inf, IC.inf, estimation_MC, CV_MC, Var_MC)[is.Call, ], Um, UM, N.tot))
  }


#-----------------------------------------------------------------------------------------------------------------
#	               	              Maximum Likelihood estimator
#-----------------------------------------------------------------------------------------------------------------

    # p.k = (p.k^-, p.k^+)
    # p : unknow
    # signature[k] = 1 si H(Y[k,]) <= 0, 0 otherwise

    log.likehood <- function(p , p.k, signature){
      gamma <- (p - p.k[,1])/(p.k[,2] - p.k[,1])
      u     <- (gamma^signature)*((1 - gamma)^(1 - signature))
      return(prod(u))
    }

#-----------------------------------------------------------------------------------------------------------------
#	   Convert an integer into a binary
#-----------------------------------------------------------------------------------------------------------------
    as.binary <- function (x) { 
      base <- 2;
      r <- numeric(inputDimension)
      for (i in inputDimension:1){ 
        r[i] <- x%%base 
	  x  <- x%/%base
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
        ttf1 <- is.dominant(Z.fail, Y.temp2, inputDimension, -1)
        if( (sum(tts1) == 0 ) & ( sum(ttf1) == 0) ){
          Y <- rbind(Y,Y.temp2)
          CP1 <- CP1 + 1
        } 
      }
      return(Y)
    }

#-------------------------------------------
#
# 			mrmEstimation
#
#-------------------------------------------

  mrmEstimation <- function(N.calls, H){

   
    V <- list()

    V <- Intersect(inputDimension, H)    

    list.set <- 0
    for(i in 1:length(V)){
      list.set[i] <- V[[i]][[2]]
    }

    u.dep   <- list()
    u.other <- list()

    u.dep[[1]] <- V[[max(which(list.set == -1))]][[1]]
    u.dep[[2]] <- V[[max(which(list.set == -1))]][[2]]

    u.other[[1]] <- V[[max(which(list.set == 1))]][[1]]
    u.other[[2]] <- V[[max(which(list.set == 1))]][[2]]

    Z.fail <- t(as.matrix(u.dep[[1]]))           #current failure points usefull in the non-dominated set 
    Z.safe <- t(as.matrix(u.other[[1]]))         #current safety points usefull in the non-dominated set 

    cp     <- length(V)  
    
    um <- 0
    uM <- 1

    Um  <- 0
    UM  <- 1

    eps   <- 1e-7
    alpha <- 0.05

    SIGN   <- 0
    ICinf  <- 0
    ICsup  <- 0
    VAR    <- 0
    CV.MLE <- 0
    MLE    <- 0

    p.hat <- 0

    X <- NULL
    Y <- NULL

    um <- prod(V[[cp]][[1]])
    uM <- 1 - prod(1 - V[[cp]][[1]])

    j  <- 1
    Um <- um
    UM <- uM
    W  <- 1:(2^(inputDimension) - 1)

    while(cp < N.calls){
      if(silent == FALSE){
        print(paste("Current number of runs =",cp));flush.console()
      }
      uu <-  Sim.non.dominated.space (1, Z.safe, Z.fail, W)

      H.u <- H(uu)
      SIGN[j] <- (1-sign(H.u))/2

      X <- rbind(X, uu)
      Y <- c(Y, H.u)

      if(H.u > 0){
        Z.safe.old <- rbind(uu, Z.safe)
        ss     <- is.dominant(Z.safe, uu, inputDimension, -1)
        Z.safe <- Z.safe[which(ss == FALSE), ]
        Z.safe <- rbind(uu, Z.safe)

        vol <- Volume.bounds(Z.safe, 1)

         Um[j+1] <- Um[j]
         if(vol >= UM[j]){
           UM[j + 1] <- UM[j]
         }else{
           UM[j + 1] <- vol
         }
         CC <- ifelse(cp == N.calls, 1, 0)

      }else{
        Z.fail.old <- rbind(uu, Z.fail)

        ff     <- is.dominant(Z.fail, uu, inputDimension, 1)
        Z.fail <- Z.fail[which(ff == FALSE), ]
        Z.fail <- rbind(uu, Z.fail)

        vol <- Volume.bounds(Z.fail, -1)

        UM[j+1] <- UM[j]
        if(vol <= Um[j]){
           Um[j+1] <- Um[j]
        }else{
          Um[j+1] <- vol
        }

      }

      cp <- cp + 1

      MLE.test <- optimize(f = log.likehood,
                               interval = c(Um[j],UM[j]),
                               maximum = TRUE,
                               signature = SIGN,
                               p.k = cbind(Um[1:j], UM[1:j])
                               )   
      MLE[j] <- as.numeric(MLE.test[1])

      VAR <- sum( 1/((MLE - Um[1:j])*(UM[1:j]- MLE)))

      bn  <- 1/VAR
      an  <- eps*VAR^(5/2)/abs( sum( 1/((MLE + eps - Um[1:j])*(UM[1:j] - MLE - eps))) - sum( 1/((MLE - Um[1:j])*(UM[1:j]- MLE)))  )

      ICinf[j]  <- MLE[j] - qnorm(1 - alpha/2)/sqrt(VAR - alpha/an)
      ICsup[j]  <- MLE[j] + qnorm(1 - alpha/2)/sqrt(VAR + alpha/an) 

      CV.MLE[j] <- 100/(sqrt(VAR)*MLE[j])   
      
      p.hat[j] <- mean(Um[1:j] + (UM[1:j] - Um[1-j])*SIGN[1:j]  )
      j <- j + 1 

    } #end of "while cp < N.calls"

    RR <- list( cbind(Um[1:(j-1)], UM[1:(j-1)], MLE[1:j-1], ICinf[1:(j-1)] , ICsup[1:(j-1)], CV.MLE[1:(j-1)], p.hat[1:(j-1)]), X, Y)
    return(RR)
  }

#-----------------------------------------------------------------------------------------------------------------
#
#
#
#
# 					               FIN	METHODE 1.1
#
#
#
#-----------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------
#
#
#
#
# 					            	FIN DE LA PARTIE 2
#
#
#
#-----------------------------------------------------------------------------------------------------------------  
  if(Method == "MRM"){
    RESULT <- mrmEstimation(N.calls, G)
  }


  if(Method == "MC"){
    RESULT <- monteCarloMonotone(N.calls)
  }

  return(RESULT)

}
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#						Fin de MRM	
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------


