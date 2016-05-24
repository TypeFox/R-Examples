
## "wild" function , global minimum at about -15.81515
fw <- function(x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80

# optimization with fewer function evaluations compared to SANN
res1 <- optim_LS(50, fw,lower = -50, upper=50, line_length = 10000)

# Upper and lower bounds should be considered carefully
res2 <- optim_LS(50, fw, lower=-Inf,upper=Inf,line_length = 10000)

# Only integer values allowed
res_int <- optim_LS(50, fw, allowed_values = seq(-50,50,by=1))

\dontrun{ 
  #plot of the function and solutions
  require(graphics)
  plot(fw, -50, 50, n = 1000, main = "Minimizing 'wild function'")
  points(-15.81515, fw(-15.81515), pch = 16, col = "red", cex = 1)
  points(res1$par, res1$ofv, pch = 16, col = "green", cex = 1)
  points(res_int$par, res_int$ofv, pch = 16, col = "blue", cex = 1)
} 

# Rosenbrock Banana function
# f(x, y) = (a-x)^2 + b*(y-x^2)^2
# global minimum at (x, y)=(a, a^2), where f(x, y)=0. 
# Usually a = 1 and b = 100 so x=1 and y=1
\dontrun{ 
  fr <- function(x,a=1,b=100) {   
    x1 <- x[1]
    x2 <- x[2]
    b*(x2 - x1*x1)^2 + (a - x1)^2
  }
  
  res3 <- optim_LS(c(-1.2,1), fr,lower = -5, upper = 5, line_length = 1000)

  # plot the surface
  x <- seq(-50, 50, length= 30)
  y <- x
  f <- function(x,y){apply(cbind(x,y),1,fr)}
  z <- outer(x, y, f)
  persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue", ticktype="detailed") -> res
  points(trans3d(1, 1, 0, pmat = res), col = 2, pch = 16,cex=2)
  points(trans3d(res3$par[1], res3$par[1], res3$ofv, pmat = res), col = "green", pch = 16,cex=2)
}

# box constraints
flb <- function(x){
  p <- length(x) 
  sum(c(1, rep(4, p-1)) * (x - c(1, x[-p])^2)^2) 
}

## 25-dimensional box constrained
\dontrun{ 
optim(rep(3, 25), flb,lower = rep(2, 25), upper = rep(4, 25),method = "L-BFGS-B") 
}
res_box <- optim_LS(rep(3, 25), flb,lower = rep(2, 25), upper = rep(4, 25),line_length = 1000) 

# one-dimensional function
\dontrun{ 
  f <- function(x)  abs(x)+cos(x)
  res5 <- optim_LS(-20,f,lower=-20, upper=20, line_length = 500)
  
  curve(f, -20, 20)
  abline(v = res5$par, lty = 4,col="green")
}  

# one-dimensional function
f <- function(x)  (x^2+x)*cos(x) # -10 < x < 10
res_max <- optim_LS(0,f,lower=-10, upper=10,maximize=TRUE,line_length = 1000) 

\dontrun{ 
  res_min <- optim_LS(0,f,lower=-10, upper=10, line_length = 1000) 
  
  curve(f, -10, 10)
  abline(v = res_min$par, lty = 4,col="green")
  abline(v = res_max$par, lty = 4,col="red")
}


# two-dimensional Rastrigin function
#It has a global minimum at f(x) = f(0) = 0.
\dontrun{ 
  Rastrigin <- function(x1, x2){
    20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
  }
  
  
  x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
  z <- outer(x1, x2, Rastrigin)
  
  res6 <- optim_LS(c(-4,4),function(x) Rastrigin(x[1], x[2]),
                   lower=-5.12, upper=5.12, line_length = 1000)
  
  # color scale
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <-
    colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  # Generate the desired number of colors from this palette
  nbcol <- 100
  color <- jet.colors(nbcol)
  # Compute the z-value at the facet centres
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  # Recode facet z-values into color indices
  facetcol <- cut(zfacet, nbcol)
  persp(x1, x2, z, col = color[facetcol], phi = 30, theta = 30)
  filled.contour(x1, x2, z, color.palette = jet.colors)
}


## Parallel computation  
## works better when each evaluation takes longer
## here we have added extra time to the computations
## just to show that it works
\dontrun{ 
  res7 <- optim_LS(c(-4,4),function(x){Sys.sleep(0.01); Rastrigin(x[1], x[2])},
                   lower=-5.12, upper=5.12, line_length = 200)
  res8 <- optim_LS(c(-4,4),function(x){Sys.sleep(0.01); Rastrigin(x[1], x[2])},
                   lower=-5.12, upper=5.12, line_length = 200, parallel = TRUE)
  res9 <- optim_LS(c(-4,4),function(x){Sys.sleep(0.01); Rastrigin(x[1], x[2])},
                   lower=-5.12, upper=5.12, line_length = 200, parallel = TRUE, 
                   parallel_type = "snow")
}
