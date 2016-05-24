## Not for public use.
characterize <- function(S1,S2){
  if(is.null(ncol(S1))){  ## simple point process without marks
    n.mark <- 0
    N1 <- length(S1);N2 <- length(S2)
    T1 <- S1;T2 <- S2
  }else{
    n.mark <- ncol(S1)-1  # dimension of marks
    N1 <- nrow(S1); N2 <- nrow(S2)
    T1 <- S1$time;T2 <- S2$time
  }
  return(list(T1=T1,T2=T2,N1=N1,N2=N2,n.mark=n.mark))
}


## Heaviside step function
hev<-function(t){ifelse(t>0,1,0)}
