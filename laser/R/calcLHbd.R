`calcLHbd` <-
function(x, r, a)
{
   if (!is.numeric(x)) stop("object x not of class 'numeric'")
   x <- rev(sort(x))
   N <- length(x)+1
   x <- c(0, x)
   LH <- ( sum(log(1:(N-1))) + ((N-2)*log(r))
        + (r*sum(x[3:N]))
        +(N*log(1-a)) - 2 * sum(log(exp(r * x[2:N])-a)));

   return(LH);
}

