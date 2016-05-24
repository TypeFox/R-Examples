WeightsRootsClenshawCurtis <- function(n){
    # n <- 3      
    wn <- rep(1,n)
    xn <- rep(0,n)
    
    if (n == 1){
        xn[1] <- 0.0;
        wn[1] <- 2.0;    
    } else {
      xn <- cos ((n-(1:n))*pi/(n-1));
#         for (i in (1:n)){
#             xn[i] <- cos ((n-i)*pi/(n-1));        
#         }
        for (i in (1:n)){
            theta = (i-1)*pi/(n-1)    
            for (j in (1:floor((n-1)/2))){
                if ((2*j) == (n-1)){b <- 1.0}
                else {b <- 2.0}            
                wn[i] <- wn[i] - b*cos(2.0*j*theta)/(4*j*j-1)
            }
        }   
        wn[1]     <- wn[1]/2;#(n-1);
        wn[2:n-1] <- 2.0*wn[2:n-1]/(n-1);
        wn[n]     <- wn[n]/(n-1);
    }
    res <- list(Xn=xn,Wn=wn)
    return(res)
}