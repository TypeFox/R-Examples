lambdasearch <-
  function(L=NULL,
           U=NULL,
           y=NULL,
           Eigenobject=NULL,
           tol=NULL,
           noisy=FALSE,
  				 eigtrunc=NULL){
    
    n <- nrow(y)  
    if(is.null(tol)){
      tol <- 10^-3 * n
    } else {
      stopifnot(is.vector(tol),
                length(tol)==1,
                is.numeric(tol),
                tol>0)    
    }
   
  # get upper bound starting value
    if(is.null(U)){
    U <- n
    while(sum(Eigenobject$values / (Eigenobject$values + U)) < 1){
      U <- U-1    
     }
    } else {
      stopifnot(is.vector(U),
                length(U)==1,
                is.numeric(U),
                U>0)
    }
    
  # get lower bound starting value
    if(is.null(L)){
      q <- which.min(abs(Eigenobject$values - (max(Eigenobject$values)/1000)))    
      
      #L <- 0
      L = .Machine$double.eps  #CJH: to avoid Inf in next statement
      
      while(sum(Eigenobject$values / (Eigenobject$values + L)) > q){
        L <- L+.05    
      }
    }  else {
      stopifnot(is.vector(L),
                length(L)==1,
                is.numeric(L),
                L>=0)
    }
    
    # create new search values    
    X1 <- L + (.381966)*(U-L)
    X2 <- U - (.381966)*(U-L)
    
    # starting LOO losses
    S1 <- looloss(lambda=X1,y=y,Eigenobject=Eigenobject, eigtrunc=eigtrunc)
    S2 <- looloss(lambda=X2,y=y,Eigenobject=Eigenobject, eigtrunc=eigtrunc)
  
    if(noisy){cat("L:",L,"X1:",X1,"X2:",X2,"U:",U,"S1:",S1,"S2:",S2,"\n") }
    
    while(abs(S1-S2)>tol){ # terminate if difference between S1 and S2 less than tolerance
    
     # update steps and use caching
      if(S1 < S2){
        U  <- X2
        X2 <- X1
        X1 <- L + (.381966)*(U-L)
        S2 <- S1
        S1 <- looloss(lambda=X1,y=y,Eigenobject=Eigenobject, eigtrunc=eigtrunc)
        
       } else { #S2 < S1
        L  <- X1
        X1 <- X2
        X2 <- U - (.381966)*(U-L)
        S1 <- S2
        S2 <- looloss(lambda=X2,y=y,Eigenobject=Eigenobject,eigtrunc=eigtrunc)
       }

      if(noisy){cat("L:",L,"X1:",X1,"X2:",X2,"U:",U,"S1:",S1,"S2:",S2,"\n") }
    }
    out <- ifelse(S1<S2,X1,X2)
    if(noisy){cat("Lambda:",out,"\n")}  
    return(invisible(out))
}


#lambdasearch(L=NULL,U=NULL,y=y,Eigenobject=Eigenobject)

