
np.est <- function(data=data, h.seq=NULL, newt=NULL, 
                   estimator="NW", kernel = "quadratic")
{

if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data) != 2)  stop("data must have 2 columns: y and t")

if ( (!is.null(h.seq)) && (sum(is.na(h.seq))  != 0) ) stop ("h.seq must be numeric")
if ( (!is.null(h.seq)) && (any(h.seq<=0)) ) stop ("h.seq must contain one ore more positive values")

if ( (!is.null(newt)) && (sum(is.na(newt))  != 0) ) stop ("newt must be numeric")
if ( (!is.null(newt)) && (any(newt<=0)) ) stop ("newt must contain one ore more positive values") 

if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")

if ((kernel != "quadratic") & (kernel != "Epanechnikov") & (kernel != "triweight") & (kernel != "gaussian") & (kernel != "uniform"))  stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")


kernel.function <- get(kernel)
n <- nrow(data)  
y <- data[,1]
t <- data[,2] 

if (is.null(newt)) newt <- t

if (is.null(h.seq)) h.seq <- np.cv(data=data, estimator=estimator, kernel=kernel)$h.opt[2,1]

num.band <- length(h.seq)

Yhat <- matrix(0, length(newt), num.band)



if ((sum(y)==Inf) | (sum(y)=="NaN")) {Yhat[,]<-0/0; return(Yhat)}

else {          
  
  Ymat <- matrix(rep(y, num.band), nrow = num.band, byrow = T)
    
  if (estimator=="NW") 
    
      for (i in 1:length(newt)) {
        
          diff <- t-newt[i]                              
          Zmat <- matrix(rep(diff, num.band), nrow = num.band, byrow = T)
          Umat <- Zmat/h.seq
         
          Kmat <- kernel.function(Umat)
          Kmat[(Kmat<0)] <- 0 

          S0 <- apply(Kmat, 1, sum)
          Kmat <- Kmat/S0
          Yhat[i,] <- apply(Ymat * Kmat, 1, sum)
    
      } # for
                                               
  else if (estimator=="LLP")
      
      for (i in 1:length(newt)) {
	
          diff <- t-newt[i]                              
          Zmat <- matrix(rep(diff, num.band), nrow = num.band, byrow = T)
          Umat <- Zmat/h.seq

          Kmat <- kernel.function(Umat)
          Kmat[(Kmat<0)] <- 0

          S0 <- apply(Kmat, 1, sum)
				  S1 <- apply(Kmat*Zmat, 1, sum)
				  S2 <- apply(Kmat*Zmat^2, 1, sum)

          Kmat <- Kmat * (S2 - Zmat*S1)/(S0*S2 - S1^2)
          Yhat[i,] <- apply(Ymat * Kmat, 1, sum)
                        
# The following is done because, when the bandwidth is "small", S1 and S2 are equal to 0, and then the corresponding estimate would be NaN.
# Nevertheless, if the bandwidth is "small" but x_i=t_k, then the corresponding estimate would be mean(y_k).

          for (j in 1:num.band){
          
              aux <- rep(0,n)      
              for (k in 1:n) if (newt[i]==t[k]) aux[k] <- 1
              num.t.newt <- sum(aux)
              if ( (num.t.newt != 0)  & ((Yhat[i,j]==Inf) | (Yhat[i,j]=="NaN"))) Yhat[i,j] <- sum(y*aux)/num.t.newt
              
          } # for j              
 
    } # for i
  
}
 
return (Yhat)              

}
