#---------------------------------------------------------------------------------#
#--------------------------------Algorithm FORM :--------------------------------#
#---------------------------------------------------------------------------------#
#' @export
FORM <- function(f, u.dep, inputDist, N.calls, eps = 1e-7, 
                 Method = "HLRF", IS = FALSE, q = 0.5, copula = "unif"){
  
  # Few error messages : 
  if (mode(f) != "function") stop("f is not a function")
  if( mode(inputDist) != "list") stop("inputDist is not a list")
  if (missing(inputDist)) stop("inputDist is missing")
  if (missing(N.calls)) stop("N.calls is missing")
  
  ndim <- length(inputDist)
  
  #-------------------------------------------------------------------------------#
  #--------------- Create the list whith paramters of the input ------------------#
  transformtionToInputSpace <- function(inputDist){
    
    InputDist <- list()
    InputDist <- inputDist
    
    for(i in 1:ndim){
      
      nparam <- length(inputDist[[i]][[2]])
      
      for(j in 1:nparam){
        InputDist[[i]]$q <- paste("q",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$p <- paste("p",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$d <- paste("d",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$r <- paste("r",InputDist[[i]][[1]], sep = "");
      }
    }
    
    return(InputDist)
  }
  
  InputDist <- transformtionToInputSpace(inputDist)
  
  #------ Transformation to the Gaussian space ------#

  Gtransf <- function(X){
    ndim <- length(InputDist)
    XU   <- pnorm(X)
    
    for(i in 1:ndim){
      XU[i] <- do.call(InputDist[[i]]$q, c(list(XU[i, drop = FALSE]), InputDist[[i]][[2]]))
    }
    return(f(XU))
  }

  #------ Gradient------#
  GRADfct <- function(func,X){
    GG <- 0
    df <- 0
    RES <- list()
    dfunc <- 0
    RES[[1]] <- 0
    
    f.X <- func(X)
    RES.2 <- X
    RES.3 <- f.X
    
    for (i in 1:ndim) {
      dX <- X
      dX[i] <- dX[i] + 1e-6 
      f.dX <- func(dX)
      dfunc[i] <- (f.dX-f.X)/1e-6
      RES.2 <- rbind(RES.2, dX)
      RES.3 <- rbind(RES.3, f.dX)
    }
    RES[[2]] <- RES.2
    RES[[3]] <- RES.3
    RES[[1]] <- dfunc

    return(RES)
  }
  #----------------------------------------------------------------------#

  x <- NULL                 # contains  the input design
  y <- NULL                 # contains f(x)
  dy <- NULL                # contains all the derivatives of f(x)

  List.res <- list()				   # List of the output
  u.new    <- 0
  cp       <- 0   			       # Numbers of calls to f
  fact.imp <- 0					       # Importance factor
  N.FORM   <- floor(N.calls*q) # Number of calls reserved to FORM
  
  u.ans <- Gtransf(u.dep)
  cp <- cp + 1
  
  x <- u.dep
  y <- u.ans

  res     <- GRADfct(Gtransf, u.dep)
  cp <- cp + ndim + 1

  x <- rbind(x,res[[2]])
  y <- rbind(y,res[[3]])
  dy <- rbind(dy,t(res[[1]]))
  
  g <- res[[1]]
  
  err   <- 1					       #distance of two tested points

  #----------------------------------Algorithme HLRF :----------------- #
  if(Method == "HLRF"){
    
    a <- g/(sqrt(g%*%g))
     
    while( (err > eps) & (cp <= (N.FORM - (2 + ndim)) ) ){ # modif BIS
      
      u.new <- (t(u.dep)%*%a)*a - u.ans*a/sqrt(g%*%g)
      err   <- sqrt(t(u.dep-u.new)%*%(u.dep-u.new))
      u.dep <- u.new
      
      res   <- GRADfct(Gtransf, u.dep)
      cp <- cp + ndim + 1

      x <- rbind(x,res[[2]])
      y <- rbind(y,res[[3]])   
      dy <- rbind(dy,t(res[[1]]))
      
      u.ans <- res[[3]][1]
      g <- res[[1]]   

      a <- ifelse( g == 0,c(0,0), g/(sqrt(g%*%g)) )
    }
    
    B 	 <- -a%*%u.dep 				 #indices of reliability
    P        <- pnorm(-B)   				 #failure probability
    fact.imp <- a^2				       
    List.res <- list(P, B, cp, u.dep, fact.imp, x, y, dy);
    
    names(List.res) <- c("pf", "beta", "compt.f", "design.point", "fact.imp", "x", "y", "dy")
    if(IS == FALSE){return(List.res)}
    }
  #-----------------------------End of the algorithm HLRF :--------------------------#
  
  #-----------------------------Algorithm Abdo-Rackwitz :--------------------------#
  if(Method == "AR"){
    
    lambda <- 2*(u.ans - (t(g)%*%u.dep))/(g%*%g)
    
    while( (err > eps) & (cp <= (N.FORM - (2 + ndim)) ) ){ # modif BIS
        
      u.new <- -lambda*g/2
      err   <- sqrt(t(u.dep-u.new)%*%(u.dep-u.new))
      u.dep <- u.new
      
      res   <- GRADfct(Gtransf, u.dep)
      cp <- cp + ndim + 1

      x <- rbind(x, res[[2]])
      y <- rbind(y, res[[3]])
      dy <- rbind(dy,t(res[[1]]))
      
      u.ans <- res[[3]][1]
      g <- res[[1]]   
      
      lambda <- 2*(u.ans - (t(g)%*%u.dep))/(g%*%g)
      if(cp >= N.FORM){break}
    }
    
    a <- g/(sqrt(g%*%g))
    B <- -a%*%u.dep  			
    P <- pnorm(-B)  				

    fact.imp <- a^2				
    List.res <- list(P, B, cp, u.dep, fact.imp, x, y, dy);
    
    names(List.res) <- c("pf", "beta", "compt.f", "design.point", "fact.imp", "x", "y", "dy")
    if(IS == FALSE){return(List.res)}
    }
  #------------------------End of the algorithm of Abdo-Rackwitz----------------------#
  
  #---------------------------Importance Sampling :---------------------------------#
  if(IS == TRUE){
    alpha  <- 0.05				
    N.IS   <- N.calls - cp			           #Number of calls reserved to the importance sampling
    norm.u <- u.dep%*%u.dep
    P.temp <- 0
    P 	 <- 0
    VAR 	 <- 0
    s	 <- 0
    Z 	 <- numeric(N.IS)

    v.temp <- matrix(rnorm(N.IS*ndim,0,1),ncol=ndim)
    v.temp <- v.temp  +  matrix(rep(u.dep,N.IS),ncol=ndim,byrow=TRUE)     #One simulate points around u.dep

    zz <- v.temp -  matrix(rep(u.dep,N.IS), ncol = ndim, byrow = TRUE)
    s  <- exp( u.dep%*%u.dep/2 - as.numeric(v.temp%*%u.dep) )

    test.v.temp <- apply(v.temp, MARGIN = 1, Gtransf)
    x <- rbind(x, v.temp)
    y <- rbind(y, matrix(test.v.temp,ncol=1))
    indic.Def <- 1*(test.v.temp <=0) 
    P.temp <- indic.Def*s

    VAR <- sum(P.temp^2) 
    cp  <- cp + N.IS					            #Numbers of call to the failure fonction inside the "sapply"
    P   <- sum(P.temp)/N.IS    			            #failure probability
    B   <- -qnorm(P)  	        			            #indices of reliability
    
    VAR <- (sum(P.temp^2)/N.IS - P^2)/(N.IS-1)              #standard error of the estimator
    
    IC.inf     <- P - qnorm(1-alpha/2)*sqrt(VAR) 		#low bound of confidence interval
    IC.sup     <- P + qnorm(1-alpha/2)*sqrt(VAR)		#high bound of confidence interval

    Inter.conf <- c(IC.inf,IC.sup)
    fact.imp   <- a^2					          
    List.res   <- list(P, B, cp, u.dep, fact.imp, VAR, Inter.conf, x, y, dy)
    
    names(List.res) <-c("pf", "beta", "compt.f", "design.point", "fact.imp", "variance", "conf", "x", "y", "dy");
    return(List.res)
  } 
  #----------------------------End Importance Sampling----------------------------#
  
  }
#---------------------------------------------------------------------------------#
#----------------------------- End of algorithm  ---------------------------------#
#------------------------------------FORM-----------------------------------------#
#---------------------------------------------------------------------------------#
