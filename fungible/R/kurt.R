#######################################################
# algorithm G2 in D. N. Joanes and C. A. Gill (1998), 
# Comparing measures of sample skewness and kurtosis. The Statistician, 
# 47, 183-189.

kurt<-function(x){
 
       kt<-function(xx){
            n <- length(xx)
            mn <- mean(xx)
            dif.x <- xx - mn
            m2 <- sum(dif.x^2)/n
            m4 <- sum(dif.x^4)/n
            b2 <- (m4/m2^2)       
            g2 <- ((n + 1) * (n - 1) * (b2 - (3 * (n - 1))/(n + 1)))/((n -2) * (n - 3))
            g2
       }     
                
           if(ncol(x)==1 || is.null(dim(x)))
             return(kt(x))
          else
           return(apply(x,2,kt))                                                      
 }         
