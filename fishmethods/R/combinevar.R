combinevar <- function(xbar = NULL,s_squared = NULL,n = NULL) {
   if(length(xbar)!=length(s_squared)|length(xbar)!=length(n)|
      length(s_squared)!=length(n)) stop ("Vector lengths are different.")
   sum_of_squares <- sum((n-1)*s_squared + n*xbar^2)
   grand_mean <- sum(n*xbar)/sum(n)
   combined_var <- (sum_of_squares - sum(n)*grand_mean^2)/(sum(n)-1)
   return(c(grand_mean,combined_var))
}
