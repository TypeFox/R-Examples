#' environ --- conducts environ analysis 
#' INPUT = network object
#' OUTPUT = input and/or output environs
#' 
#' M. Lau July 2011 | DEH edited Feb 2013
#' ---------------------------------------------------

enaEnviron <- function(x,input=TRUE,output=TRUE,type='unit',err.tol=1e-10,balance.override=FALSE){
                                        #check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}
                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }
                                        #Assume 'rc' orientation of flows
                                        #Don't transpose flows for calculations, they will be transposed in enaFlow
                                        #calculate enaFlow with RC input
  user.orient <- get.orient()
  set.orient('rc')
  Flow <- enaFlow(x)
  set.orient(user.orient)
                                        #Transpose for calculations in Patten school
  Flow$G <- t(Flow$G)
  Flow$GP <- t(Flow$GP)
  Flow$N <- t(Flow$N)
  Flow$NP <- t(Flow$NP)
                                        #Unit environ calculations
  if (input){
                                        #Input perspective
    EP <- list()
    for (i in (1:nrow(Flow$NP))){
      dNP <- diag(Flow$NP[i,]) #diagonalized N matrix
      EP[[i]] <- dNP %*% Flow$GP         # calculate internal environ flows
      EP[[i]] <- EP[[i]] - dNP      # place negative environ throughflows on the principle diagonal
      EP[[i]] <- cbind(EP[[i]],apply((-EP[[i]]),1,sum)) #attach z column
      EP[[i]] <- rbind(EP[[i]],c(apply((-EP[[i]]),2,sum))) #attach y row
      EP[[i]][nrow(EP[[i]]),ncol(EP[[i]])] <- 0 #add zero to bottom right corner to complete matrix
      EP[[i]][abs(EP[[i]]) < err.tol] <- 0 #ignore numerical error
                                        #add labels to matrices
      labels <- c(rownames(Flow$GP))
      colnames(EP[[i]]) = c(labels, 'z')
      rownames(EP[[i]]) = c(labels, 'y')
    }
                                        #add environ names
    names(EP) <- labels
  }
  if (output){
                                        #Output perspective
    E <- list()
    for (i in (1:nrow(Flow$N))){
      dN <- diag(Flow$N[,i]) #diagonalized N matrix
      E[[i]] <- Flow$G %*% dN
      E[[i]] <- E[[i]] - dN
      E[[i]] <- cbind(E[[i]], apply((-E[[i]]), 1, sum))
      E[[i]] <- rbind(E[[i]], c(apply((-E[[i]]), 2, sum)))
      E[[i]][nrow(E[[i]]), ncol(E[[i]])] <- 0 
      E[[i]][abs(E[[i]]) < err.tol] <- 0 
                                        #add labels to matrices
      labels <- c(rownames(Flow$G))
      colnames(E[[i]]) = c(labels, 'z')
      rownames(E[[i]]) = c(labels, 'y')
    }
                                        #add environ names
    names(E) <- labels
  }
                                        #Realized environ calculations
  if (type == 'realized'){
                                        #Input perspective
    if (input){
      for (i in (1:nrow(Flow$N))){
        EP[[i]] <- EP[[i]]*unpack(x)$y[i] #Construct realized environ
      }
    }
                                        #Output perspective
    if (output){
      for (i in (1:nrow(Flow$N))){
        E[[i]] <- E[[i]]*unpack(x)$z[i]
      }
    }
  }
                                        #Wrap-up output into list  
  if (input & output){
    out <- list('input' = EP,'output' = E)
  } else if (input & output == FALSE){
    out <- EP
  } else if (input == FALSE & output){
    out <- E
  }      
  
  if (type != 'unit' && type!= 'realized'){
    print('WARNING: Invalid input in type, input ignored')
  }
                                        #re-orient matrices
  if (user.orient == 'rc'){
    for (i in 1:length(out)){
      for (j in 1:length(out[[i]])){
        out[[i]][[j]] <- t(out[[i]][[j]])
      }
    }
  }else{}
                                        #output
  return(out)
}

