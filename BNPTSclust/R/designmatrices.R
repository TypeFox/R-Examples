designmatrices <-
function(level,trend,seasonality,deg,T,n,fun){
    
    # Function that generates the design matrices of the clustering
    # algorithm based on the parameters that the user wants to consider,
    # i.e. level, polinomial trend and/or seasonal components. It also 
    # returns the number of parameters that are considered and not
    # considered for clustering. Since this function is for internal use,
    # its arguments are taken directly from the clustering functions.
    #
    # IN:
    #
    # level       <- Variable that indicates if the level of the time 
    #                series will be considered for clustering. If 
    #                level = 0, then it is omitted. If level = 1, then it 
    #                is taken into account.
    # trend       <- Variable that indicates if the polinomial trend of 
    #                the model will be considered for clustering. If 
    #                trend = 0, then it is omitted. If trend = 1, then it
    #                is taken into account.
    # seasonality <- Variable that indicates if the seasonal components 
    #                of the model will be considered for clustering. 
    #                If seasonality = 0, then they are omitted. If 
    #                seasonality = 1, then they are taken into account.
    # deg         <- Degree of the polinomial trend of the model.
    # T           <- Number of periods of the time series.
    # n           <- Number of time series.
    # fun         <- Clustering function being used. 
    #
    # OUT:
    #
    # Z <- Design matrix of the parameters not considered for clustering.
    # X <- Design matrix of the parameters considered for clustering.
    # p <- Number of parameters not considered for clustering.
    # d <- Number of parameters considered for clustering.
    
    if(fun == "tseriesca"){
      M <- matrix(0,T,1+deg)      # Matrix with all components.
      M[,1] <- 1                     # Level components.
      for(i in 1:deg){               # Trend components.
        M[,i+1] <- seq(T)^i
      }
      
      if(level == 0 & trend == 0){
        p <- 1+deg
        d <- 0
        Z <- M
        return(list(p=p,d=d,Z=Z))
      }
      
      if(level == 1 & trend == 0){
        p <- deg
        d <- 1
        Z <- M[,(2:(deg+1))]
        X <- as.matrix(M[,1])
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 0 & trend == 1){
        p <- 1
        d <- deg
        Z <- as.matrix(M[,1])
        X <- M[,(2:(deg+1))]
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 1 & trend == 1){
        p <- 0
        d <- 1+deg
        X <- M
        return(list(p=p,d=d,X=X))
      }
      
    }
    
    if(fun == "tseriescm"){
      
      M <- matrix(0,T,1+deg+11)      # Matrix with all components.
      M[,1] <- 1                     # Level components.
      for(i in 1:deg){               # Trend components.
        M[,i+1] <- seq(T)^i
      }
      # Seasonal components
      num <- floor(T/12)                     # Number of years present in the data
      
      if (num < 1){                          # If the number of months in the data is less than 12, the design matrix is filled this way
        X2 <- diag(1,(T-1))
        X2 <- cbind(X2,matrix(0,(T-1),1))
        X <- rbind(X2,matrix(0,1,T))  
      }else{
        X21 <- rbind(diag(1,11),matrix(0,1,11))  # Matrix that contains the indicator functions for the 11 months and one row of zeros to avoid singularity problems in the design matrix 
        X2 <- X21
        resid <- T %% 12                         # Number of the year (num+1) present in the data
        
        if (num >= 2){
          for (i in 2:num){
            X2 <- rbind(X2,X21)   
          }  
        }
      }
        
      M[,((deg+2):(1+deg+11))] <- rbind(X2,X21[0:resid,])          
        
      if(level == 0 & trend == 0 & seasonality == 0){
        p <- 1+deg+11
        d <- 0
        Z <- M
        return(list(p=p,d=d,Z=Z))
      }
      
      if(level == 0 & trend == 0 & seasonality == 1){
        p <- 1+deg
        d <- 11
        Z <- M[,(1:(deg+1))]
        X <- M[,((deg+2):(1+deg+11))]
        return(list(p=p,d=d,Z=Z,X=X))        
      }
      
      if(level == 0 & trend == 1 & seasonality == 0){
        p <- 1+11
        d <- deg
        Z <- cbind(M[,1],M[,(deg+2):(1+deg+11)])
        X <- M[,(2:(deg+1))]
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 1 & trend == 0 & seasonality == 0){
        p <- deg+11
        d <- 1
        Z <- M[,(2:(1+deg+11))]
        X <- as.matrix(M[,1])
        return(list(p=p,d=d,Z=Z,X=X))        
      }
      
      if(level == 1 & trend == 1 & seasonality == 0){
        p <- 11
        d <- 1+deg
        Z <- M[,(deg+2):(1+deg+11)]
        X <- M[,(1:(deg+1))]
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 1 & trend == 0 & seasonality == 1){
        p <- deg
        d <- 1+11
        Z <- M[,(2:(deg+1))]
        X <- cbind(M[,1],M[,((deg+2):(1+deg+11))])
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 0 & trend == 1 & seasonality == 1){
        p <- 1
        d <- deg+11
        Z <- as.matrix(M[,1])
        X <- M[,(2:(1+deg+11))]
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 1 & trend == 1 & seasonality == 1){
        p <- 0
        d <- 1+deg+11
        X <- M
        return(list(p=p,d=d,X=X))
      }
       
    }
    
    if(fun == "tseriescq"){
      
      M <- matrix(0,T,1+deg+3)      # Matrix with all components.
      M[,1] <- 1                     # Level components.
      for(i in 1:deg){               # Trend components.
        M[,i+1] <- seq(T)^i
      }
      # Seasonal components
      num <- floor(T/4)                     # Number of years present in the data
      
      if (num < 1){                          # If the number of months in the data is less than 12, the design matrix is filled this way
        X2 <- diag(1,(T-1))
        X2 <- cbind(X2,matrix(0,(T-1),1))
        X <- rbind(X2,matrix(0,1,T))  
      }else{
        X21 <- rbind(diag(1,3),matrix(0,1,3))  # Matrix that contains the indicator functions for the 11 months and one row of zeros to avoid singularity problems in the design matrix 
        X2 <- X21
        resid <- T %% 4                         # Number of the year (num+1) present in the data
        
        if (num >= 2){
          for (i in 2:num){
            X2 <- rbind(X2,X21)   
          }  
        }
      }
      
      M[,((deg+2):(1+deg+3))] <- rbind(X2,X21[0:resid,])          
      
      if(level == 0 & trend == 0 & seasonality == 0){
        p <- 1+deg+3
        d <- 0
        Z <- M
        return(list(p=p,d=d,Z=Z))
      }
      
      if(level == 0 & trend == 0 & seasonality == 1){
        p <- 1+deg
        d <- 3
        Z <- M[,(1:(deg+1))]
        X <- M[,((deg+2):(1+deg+3))]
        return(list(p=p,d=d,Z=Z,X=X))        
      }
      
      if(level == 0 & trend == 1 & seasonality == 0){
        p <- 1+3
        d <- deg
        Z <- cbind(M[,1],M[,(deg+2):(1+deg+3)])
        X <- M[,(2:(deg+1))]
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 1 & trend == 0 & seasonality == 0){
        p <- deg+3
        d <- 1
        Z <- M[,(2:(1+deg+3))]
        X <- as.matrix(M[,1])
        return(list(p=p,d=d,Z=Z,X=X))        
      }
      
      if(level == 1 & trend == 1 & seasonality == 0){
        p <- 3
        d <- 1+deg
        Z <- M[,(deg+2):(1+deg+3)]
        X <- M[,(1:(deg+1))]
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 1 & trend == 0 & seasonality == 1){
        p <- deg
        d <- 1+3
        Z <- M[,(2:(deg+1))]
        X <- cbind(M[,1],M[,((deg+2):(1+deg+3))])
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 0 & trend == 1 & seasonality == 1){
        p <- 1
        d <- deg+3
        Z <- as.matrix(M[,1])
        X <- M[,(2:(1+deg+3))]
        return(list(p=p,d=d,Z=Z,X=X))
      }
      
      if(level == 1 & trend == 1 & seasonality == 1){
        p <- 0
        d <- 1+deg+3
        X <- M
        return(list(p=p,d=d,X=X))
      }
      
    }
  }
