## =============================================================================
## Example 5 from Shampine et al.
## 
##   Falkner-Skan BVPs are discussed in T. Cebeci and H.B. Keller, 
##   Shooting and parallel shooting methods for solving the Falkner-Skan
##   boundary-layer equation, J. Comp. Phy., 7 (1971) 289-300.  This is 
##   the positive wall shear case for which the parameter beta is known 
##   and the problem is to be solved for a range of the parameter.  This
##   is the hardest case of the table in the paper.
##
## =============================================================================

require(bvpSolve)
falkner <- function(x, y, beta)  
  list( c(y[2],
          y[3],
          -y[1]*y[3] - beta*(1 - y[2]^2))
             )

beta     <- 0.5
infinity <- 6

Sol <- bvptwp(func = falkner, x = seq(0,infinity, by=0.1), parms=beta,
              yini = c(f=0, df=0, df2=NA), yend = c(NA, 1, NA) )
              
plot(Sol)
Sol[1,"df2"]  #  0.92768

# solving for a range of the parameter

bet <- seq(0,2,by=0.2)
solbet <- NULL

for (beta in bet) {
  Sol <- bvptwp(func = falkner, x = seq(0,infinity, by=0.1), parms=beta,
              yini = c(f=0, df=0, df2=NA), yend = c(NA, 1, NA))
              
  solbet<-cbind(solbet,Sol[,2])
}

matplot(solbet, type = "l")
legend("topleft", col = 1:length(bet), lty = 1,
     legend = bet, title = "beta")
