
###############################################################################
# ParMLE                                                                      #
#                                                                             #
# Maximum likelihood estimation of item parameters for cognitive diagnostic   #
# models.                                                                     #
#                                                                             #
# Input:                                                                      #
# (1) Y: a matrix of binary responses (1=correct, 0=incorrect). Rows          #
#               represent persons and columns represent items.                #
# (2) Q: the Q-matrix of the test. Rows represent items and columns represent #
#        attributes.                                                          #
# (3) alpha: examinee attribute profiles. Rows represent persons and columns  #
#            represent attributes.                                            #
# (4) model: currently has three options, "DINA", "DINO", "NIDA", "GNIDA",    #
#            and "RRUM".                                                      #
#                                                                             #
# Output:                                                                     #
#     For DINA, DINO, and NIDA models:                                        #
# (1) slip: a vector of slip parameters                                       #
# (2) guess: a vector of guessing parameters                                  #
# (3) se.slip: a vector of standard error for slip parameters.                #
# (4) se.guess: a vector of standard error for guessing parameters.           #
#     For GNIDA:                                                              #
# (1) slip: a matrix (# items by # attributes) of slip parameters.            #
# (2) guess: a matrix of guessing parameters                                  #
# (3) se.slip: a matrix of standard error for slip parameters.                #
# (4) se.guess: a matrix of standard error for guessing parameters.           #
#     For RRUM:                                                               #
# (1) pi: a vector of pi parameters for each item                             #
# (2) r: a matrix (# items by # attributes) of r parameters                   #
# (3) se.pi: a vector of standard errors for pi parameters.                   #
# (4) se.r: a matrix of standard errors for r parameters.                     #
#                                                                             #
###############################################################################


ParMLE <- function(Y, Q, alpha, model=c("DINA", "DINO", "NIDA", "GNIDA", "RRUM"))
{
  #####
  # 1 #
  ##### Check dimension consistency and convert data to the right formats 
  
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  check <- NULL
  check <- CheckInput(Y, Q)  
  if (!is.null(check)) return(warning(check))
  
  model <- match.arg(model)
  
  #####
  # 2 #
  ##### Estimation
  
  nitem <- dim(Y)[2]
  nperson <- dim(Y)[1]
  natt <- dim(Q)[2]
  
  if (model == "DINA")
  {
    slip <- se.slip <- guess <- se.guess <- matrix(NA, nitem, 1)
    
    for (i in 1:nitem)
    {
      ita <- NULL
      for (j in 1:nperson)
      {
        ita[j] <- prod(alpha[j, ] ^ Q[i, ])
      }
      
      slip[i] <- sum((1 - Y[ , i]) * ita) / sum(ita)
      se.slip[i] <- sqrt(slip[i] * (1 - slip[i]) / sum(ita))
      guess[i] <- sum(Y[ , i] * (1 - ita)) / sum(1 - ita)   
      se.guess[i] <- sqrt(guess[i] * (1 - guess[i]) / sum(1 - ita))
    }        
        
  } else if (model == "DINO")
  {
    slip <- se.slip <- guess <- se.guess <- matrix(NA, nitem, 1)
    
    for (i in 1:nitem)
    {
      omega <- NULL
      for (j in 1:nperson)
      {
        omega[j] <- 1 - prod((1 - alpha[j, ]) ^ Q[i, ])
      }
            
      slip[i] <- sum((1 - Y[ , i]) * omega) / sum(omega)
      se.slip[i] <- sqrt(slip[i] * (1 - slip[i]) / sum(omega))
      guess[i] <- sum(Y[ , i] * (1 - omega)) / sum(1 - omega)   
      se.guess[i] <- sqrt(guess[i] * (1 - guess[i]) / sum(1 - omega))
    }        
        
  } else if (model == "NIDA")
  {
    slip <- se.slip <- guess <- se.guess <- rep(NA, natt)
  
    alpha_mesh <- matrix(t(alpha), nrow=nitem*nperson, ncol=natt, byrow=T)
    Q_mesh <- matrix(apply(t(Q), 2, function(x) rep(x, nperson)), nrow=nitem*nperson, ncol=natt, byrow=T)
    Y_mesh <- matrix(Y, nrow=nitem*nperson, ncol=1)
    
    derlike_NIDA <- function(par, alpha_mesh, Q_mesh, Y_mesh)
    {      
      dls <- dlg <- rep(0, natt)
      
      s <- matrix(par[1:natt], nrow=nitem*nperson, ncol=natt, byrow=T)
      g <- matrix(par[(natt+1):(2*natt)], nrow=nitem*nperson, ncol=natt, byrow=T)
      
      for (k in 1:natt)
      { 
        P <- apply((1 - s) ^ (alpha_mesh * Q_mesh) * g ^ ((1 - alpha_mesh) * Q_mesh), 1, prod)
        dls[k] <- sum((P - Y_mesh) * alpha_mesh[, k] * Q_mesh[, k] / (1 - P)) / (1 - par[k])
        dlg[k] <- sum((Y_mesh - P) * (1 - alpha_mesh[, k]) * Q_mesh[, k] / (1 - P)) / par[natt + k]        
      }     
      return(c(dls, dlg))
    }
    
    # Initial values for optimization
    
    p0 <- rep(0.3, 2 * natt)
    for (k in 1:natt)
    {
    	Q_temp <- rep(0, natt)
    	Q_temp[k] <- 1
    	if (sum(apply(Q, 1, function(x) all(x==Q_temp))) > 0)
    	{
    		p0[k] <- sum(alpha_mesh[,k]==1 & apply(Q_mesh, 1, function(x) all(x==Q_temp)) & Y_mesh==0
    		) / sum(alpha_mesh[,k]==1 & apply(Q_mesh, 1, function(x) all(x==Q_temp)))
    		
    		p0[natt + k] <- sum(alpha_mesh[,k]==0 & apply(Q_mesh, 1, function(x) all(x==Q_temp)) & Y_mesh==1
    		) / sum(alpha_mesh[,k]==0 & apply(Q_mesh, 1, function(x) all(x==Q_temp)))
		 }
    }
    p0[p0 < 0.05 | p0 > 0.5] <- 0.3
    
    # Optimization
    
    ans <- dfsane(par=p0, fn=derlike_NIDA, alpha_mesh=alpha_mesh, Q_mesh=Q_mesh, Y_mesh=Y_mesh, control=list(maxit=300, trace=FALSE))
    slip <- ans$par[1:natt]
    guess <- ans$par[(natt + 1):(2 * natt)]
    
    slip[slip < 0] <- 0; slip[slip > 1] <- 1
    guess[guess < 0] <- 0; guess[guess > 1] <- 1
    
    for (k in 1:natt)
    {
      se.slip[k] <- sqrt(slip[k] * (1 - slip[k]) / (sum(alpha[, k]) * nitem))
      se.guess[k] <- sqrt(guess[k] * (1 - guess[k]) / ((nperson - sum(alpha[,k])) * nitem))
    } 

  } else if (model == "GNIDA")
  {
    slip <- se.slip <- guess <- se.guess <- matrix(0, nitem, natt)
    
    derlike_GNIDA <- function(par, alpha, Q_j, Y_j)     # Q_j and Y_j are vectors
    {              
      ind <- which(Q_j == 1)
      K1 <- sum(Q_j)
      Q_mesh <- matrix(rep(Q_j[ind], nperson), nrow=nperson, ncol=K1, byrow=T)
      Y_mesh <- matrix(Y_j, nrow=nperson, ncol=K1)
      s <- matrix(par[1:K1], nrow=nperson, ncol=K1, byrow=T)
      g <- matrix(par[(K1+1):(2*K1)], nrow=nperson, ncol=K1, byrow=T)      
	  P <- apply(((1 - s) ^ (alpha[ ,ind] * Q_mesh) * g ^ ((1 - alpha[ ,ind]) * Q_mesh)), 1, prod)
	  P_mesh <- matrix(P, nrow=nperson, ncol=K1)
      dls <- colSums((P_mesh - Y_mesh) * alpha[ ,ind] * Q_mesh / (1 - P_mesh) / (1 - s))
      dlg <- colSums((Y_mesh - P_mesh) * (1 - alpha[, ind]) * Q_mesh / (1 - P_mesh) / g)
      return(c(dls, dlg))
    }
  
    
    for (j in 1:nitem)
    {       	
      # Initial values for optimization    
	  p0 <- rep(0.3, 2 * natt)
	  if (sum(Q[j,] == 1) == 1)
	  {
	  	k <- which(Q[j,] == 1)
   	  	p0[k] <- sum(alpha[,k]==1 & Y[,j]==0) / sum(alpha[,k]==1)
   	  	p0[natt + k] <- sum(alpha[,k]==0 & Y[,j]==1) / sum(alpha[,k]==0)
      }
      p0[p0 < 0.05 | p0 > 0.5] <- 0.3      
    
      # Optimization    	
      ind <- which(Q[j,] == 1)
      K1 <- sum(Q[j,])
      ans <- dfsane(par=c(p0[ind], p0[natt+ind]), fn=derlike_GNIDA, alpha=alpha, Q_j=Q[j,], Y_j=Y[,j], control=list(maxit=300, trace=FALSE))
      slip[j, ind] <- ans$par[1:K1]
      guess[j, ind] <- ans$par[(K1 + 1):(2 * K1)]
    }
    
    slip[slip < 0] <- 0; slip[slip > 1] <- 1
    guess[guess < 0] <- 0; guess[guess > 1] <- 1
    
    for (j in 1:nitem)
    {
      for (k in 1:natt)
      {
        se.slip[j, k] <- sqrt(slip[j, k] * (1 - slip[j, k]) / sum(alpha[, k]))
        se.guess[j, k] <- sqrt(guess[j, k] * (1 - guess[j, k]) / (nperson - sum(alpha[, k])))
      } 
    }   
        
  } else if ( model == "RRUM") 
  {
  	pi <- se.pi <- rep(NA, nitem)
	r <- se.r <- matrix(0, nitem, natt)
	
  	derlike_RRUM <- function(par, alpha, Q_j, Y_j)
	{
	  dlpi <- dlr <- 0
	  ind <- which(Q_j == 1)
	  K1 <- sum(Q_j)
	  P <- rep(par[1], nperson) * apply(matrix(rep(par[2:(K1+1)], nperson), nperson, K1, byrow=TRUE) ^ 
	  			((1 - alpha[, ind]) * matrix(rep(Q_j[ind], nperson), nperson, K1, byrow=TRUE)), 1, prod)
	  dlpi <- dlpi + (Y_j - P) / (1 - P)
	  dlr <- dlr + (1 - alpha[, ind]) * matrix(rep(Q_j[ind], nperson), nperson, K1, byrow=TRUE) * 
	  		(Y_j - P) / (1 - P)
	  return(apply(cbind(dlpi, dlr), 2, sum))
	}
	
	pi0 <- rep(0.8, nitem)
	r0 <- matrix(0.8, nitem, natt) * Q
	for (j in 1:nitem)
	{
	  p0 <- c(pi0[j], r0[j, which(Q[j, ] == 1)])
	  ans <- dfsane(par=p0, fn=derlike_RRUM, alpha=alpha, Q_j=Q[j,], Y_j=Y[,j], control=list(maxit=300, trace=FALSE))
	  pi.temp <- ans$par[1]
	  r.temp <- ans$par[2:(sum(Q[j,]) + 1)]
	  pi[j] <- pi.temp
	  r[j,which(Q[j,] == 1)] <- r.temp
	}

    pi[pi < 0] <- 0; pi[pi > 1] <- 1
    r[r < 0] <- 0; r[r > 1] <- 1
    
    for (j in 1:nitem)
    {
      se.pi[j] <- sqrt(pi[j] * (1 - pi[j]) / sum(apply(alpha, 1, function(x) prod(x ^ Q[j, ])) == 1))
      for (k in 1:natt) if (Q[j, k] == 1) se.r[j, k] <- sqrt(r[j, k] * (1 - r[j, k]) / (nperson - sum(alpha[, k])))
    }   
      
  } else
  {
    return(warning("Model specification is not valid."))
  }
   
  if (model %in% c("DINA", "DINO", "NIDA", "GNIDA")) output <- list(slip=slip, guess=guess, se.slip=se.slip, se.guess=se.guess, model=model, Q=Q, Y=Y, alpha=alpha)
  if (model == "RRUM") output <- list(pi=pi, r=r, se.pi=se.pi, se.r=se.r, model=model, Q=Q, Y=Y, alpha=alpha)

  class(output) <- "ParMLE"
  return(output)
  
}

