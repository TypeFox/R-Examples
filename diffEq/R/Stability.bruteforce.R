## =============================================================================
## Stability by the brute force method
## =============================================================================
stability.bruteforce <- function (
           Rez = seq(-2,2,by=0.02), 
           Imz = seq(-2,2,by=0.02), 
           func = function (z) return(abs(1 + h*z) <=1),
           fill = "grey", 
           cex = 1.5,
           add = FALSE,
           ...)    { 
  h <- 1
  if (!add) 
    plot(0, type = "n", xlim = range(Rez), ylim = range(Imz),# asp=TRUE,
       xlab = "Re(z)", ylab = "Im(z)", ...)  
  gg  <- expand.grid(Rez,Imz)
  zz <- complex( real = gg[,1], imaginary = gg[,2])

  for(i in 1:nrow(gg)){
    if(func(zz[i]))
    points(gg[i,1],gg[i,2],pch = 18, col = fill, cex= cex)
  }
  abline(h=0)
  abline(v=0)
}


stabfuncExpRK <- function(z, order = 1) {
      h  <- 1
      ss <- 1
      for (p in 1: order) ss <- ss + (h*z)^p / factorial(p)
      return (abs(ss) <= 1) 
}

stabfuncImpRK <- function(x, ...) {   # as in Hairer...
  if (is.list(x)) x <- x$ID
  else if (!is.character(x)) stop ("'x' should be a list of type rkMethod or a name of an implicit RK")

  ## Radau IIA order 5
  if (x == "irk5r")
  fun  <-  function (z) return(abs((1 + 2*z/5 + z^2/20) /
                            (1 - 3*z/5 + 3*z^2/20 - z^3/60)) <= 1)  
  ## Hammer - Hollingsworth coefficients , order 4
  else if (x == "irk4hh")
  fun  <-  function (z) return(abs((1 + z/2 + z^2/12) /
                                   (1 - z/2 + z^2/12)) <= 1)  
  ## Kuntzmann and Butcher order 6
  else if (x == "irk6kb")
  fun  <-  function (z) return(abs((1 + 2*z/3 + z^2/5 + z^3/30 + z^4/360) /
                            (1 - z/3 + z^2/30)) <= 1)  
  ## Lobatto order 4
  else if (x == "irk4l")
  fun  <-  function (z) return(abs((1 + 3*z/4 + z^2/4 + z^3/24) /
                            (1 - z/4)) <= 1)  

  ## Lobatto order 6      
  else if (x == "irk6l")

  return(fun)
 
}

