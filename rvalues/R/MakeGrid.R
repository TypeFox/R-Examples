MakeGrid <- function(nunits, type = "log", ngrid = NULL, lower = 1/nunits,
                     upper = 1 - lower)  {
  if((lower <=0) | (upper >= 1) | (lower >= upper)) {
     stop("Both lower and upper must be numbers strictly between zero and one. 
           Also, lower must be greater than upper.")
  }                   
  switch(type,
         uniform={
            if(is.null(ngrid)) {
              if(nunits <= 1000) {
                ngrid <- nunits
              }
              else{
                ## Size of grid grows according to
                ## f(x) = 1000 + 25*log(x - 1000)*(x - 1000)^(1/4)
                
                nn <- nunits - 1000
                ngrid <- 1000 + ceiling(25*exp(.25*log(nn))*log(nn))
              }
            }
             grid <- seq(lower,upper,length=ngrid)
         },
         log={
             if(is.null(ngrid)) {
               if(nunits <= 1000) {
                 ngrid <- nunits
               }
               else{
                 ## Size of grid grows according to
                 ## f(x) = 1000 + 25*log(x - 1000)*(x - 1000)^(1/4)
                 
                 nn <- nunits - 1000
                 ngrid <- 1000 + ceiling(25*exp(.25*log(nn))*log(nn))
               } 
             }
            
             grid <- exp( seq(log(lower),log(upper), length=ngrid ) )
         },
         log.symmetric={
            if(is.null(ngrid)) {
               if(nunits <= 1000) {
                 ngrid <- nunits
               }
               else{
                 ## Size of grid grows according to
                 ## f(x) = 1000 + 25*log(x - 1000)*(x - 1000)^(1/4)
                
                 nn <- nunits - 1000
                 ngrid <- 1000 + ceiling(25*exp(.25*log(nn))*log(nn))
               }
             }
             ### ignore the value for upper in this case (due to the symmetry requirement)
             grid1 <- exp( seq(log(lower),log(1/2), length=ceiling((ngrid-1)/2) + 1 ) )
             grid2 <- rev(1 - grid1[-length(grid1)])
           
             # if ngrid is odd, length(grid) = ngrid
             # if ngrid is even, length(grid) = ngrid + 1
             grid <- c(grid1,grid2)
        },
  )
  return(grid)
}
