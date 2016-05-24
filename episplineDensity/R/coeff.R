coeff <- function (x, params) 
{
#      x: numeric data item (length 1)
# params: list of epiparameters
#
# Evaluate coefficients at x using params.
#
if (missing (params))
    stop ("Please specify epiparameters")

N <- params$Ndiscr
m0 <- params$m0
mN <- params$mN
ord <- params$order

Delta <- (mN - m0)/N
c.out <- numeric (N * (ord + 2) +1) # return value

if ((x - m0) %% Delta == 0) {   # if x falls on segment endpoint
    k <- round((x-m0)/Delta)+1; 
    c.out[k] <- 1;
}
else {
   k <- ceiling((x-m0)/Delta);
   xk <- x - ((k-1)*Delta + m0);
   c.out[N+1+(k-1)*(ord+1)+1] <- 1; 
   c.out[N+1+(k-1)*(ord+1)+2] <- xk; 
   c.out[N+1+(k-1)*(ord+1)+3] <- xk^2;
}
return (c.out)

}
