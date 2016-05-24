"simm.levy" <- function (date = 1:500, mu = 2, l0 = 1, x0 = c(0, 0), id = "A1",
   burst = id, typeII = TRUE)
{
   if (typeII)
       class(date) <- c("POSIX", "POSIXct")
   n <- length(date)
   dt <- c(diff(unclass(date)))
   if (all(dt - dt[1] > 1e-07))
       stop("the time lag between relocations should be constant")
   ang <- runif(n - 2, -pi, pi)
   v = dt * (l0 * (runif(n - 1)^(1/(1 - mu))))
   ang = cumsum(c(runif(1, 0, 2 * pi), ang))
   si = c(x0[2], x0[2]+cumsum(v * sin(ang)))
   co = c(x0[1], x0[1]+cumsum(v * cos(ang)))
   res <- as.ltraj(data.frame(co, si), date, id, burst, typeII = typeII)
   return(res)
}
