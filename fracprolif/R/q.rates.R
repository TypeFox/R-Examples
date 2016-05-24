# See Powell, 1956
setClass("qrate",
  representation(stderr = "numeric", n="numeric", df="numeric"),
  contains = "numeric")

q.rates <- function(dist, est)
{
    q.rates.inner <- function(dist, params)
    {
        params   <- as.list(params)
        Q        <- params$Q
        params$Q <- NULL

        g <- function(t, v)  { params <- params;
                               params$x <- t;
                               exp(-v * t) *
                                 do.call(paste("d",dist, sep=''), params) }
        f <- function(v) {integrate(g, 0, Inf, v)$value - 0.5}
           # May need to adjust max value of integration range from Inf to 100
        r <- uniroot(f, lower=0.0001, upper=1)$root
        new("qrate", c(r,  Q*r/(1-Q) ) )
    }

   v <- q.rates.inner(dist, coef(est))
   
   # Now find stderr
   n <- est@details$counts[1]
   
   gd <- grad(function(x) {q.rates.inner("emg", x)[1]}, coef(est))[1:3]
   
   gq <- grad(function(x) {q.rates.inner("emg", x)[2]}, coef(est))
   
   attr(v, "stderr") <- c(
       sqrt(t(gd) %*% vcov(est)[1:3,1:3] %*% gd / n),
       sqrt(t(gq) %*% vcov(est) %*% gq / n)
   )
   
   attr(v, "n")  <- n
   attr(v, "df") <- n-2
   
   attr(v, "names") <- c("d", "q")
   
   v
}

setMethod("show","qrate", 
  function(object)
  {
     cat("         d          q\n")
     cat(object, "\n")
  }
)

setMethod("summary","qrate", 
  function(object)
  { 
     cat("         d          q\n")
     cat(object)
     cat("\n\n")
     cat("       se(d)        se(q)\n")
     cat(object@stderr)
     cat("\n\n")
     cat("n  = ")
     cat(object@n)
     cat("\n")
     cat("df = ")
     cat(object@df)
  }
)
