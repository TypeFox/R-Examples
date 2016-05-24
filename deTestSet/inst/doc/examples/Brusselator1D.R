## =============================================================================
##
## The brusselator problem in 1D  
##
##		A ----> U           (1)
## 	  B + U ----> R + V   (2)
##		V + 2 U ----> 3 U   (3)
##		U ----> S.          (4)
##
## =============================================================================

require(ReacTran)
require(deTestSet)

#-----------------------------
# model function
#-----------------------------

brusselator1D <- function(t, y, parms) {

   U <- y[1:Nx]
   V <- y[(Nx+1):(2*Nx)]

   dU <- A + U^2*V - (B+1)*U +                  # reaction
         tran.1D (C=U, C.up = 1,  C.down = 1,   # transport
                       D = Dx, dx=Grid)$dC      
                       
   dV <- B*U - U^2*V +                          # reaction
         tran.1D (C=V, C.up = 3,  C.down = 3,
                       D = Dx, dx=Grid)$dC     

   list(c(dU,dV))                                    
}

#-----------------------------
# boxes and the grid
#-----------------------------

Nx   <- 200    # number of boxes
Grid <- setup.grid.1D(x.up = 0, x.down = 1, N = Nx)

#-----------------------------
# Parameters, initial values 
#-----------------------------

A  <- 1
B  <- 3
Dx <- 0.02       # diffusion coefficient

# initial conditions
uini <- 1 + sin(2*pi*Grid$x.mid)
vini <- rep(x = 3, times = Nx)

yini <- c(uini, vini)

#-----------------------------
# solve model 
#-----------------------------

times <- seq(from = 0, to = 10, by = 0.01)

print(system.time(
  out <- ode.1D(y = yini, parms = NULL, func = brusselator1D, nspec = 2,
     dimens = Nx, times = times)
))

# works only for deSolve version 1.10
print(system.time(
  out <- ode.1D(y = yini, parms = NULL, func = brusselator1D, nspec = 2,
     dimens = Nx, times = times, method = mebdfi, restructure = TRUE)
))


#-----------------------------
# plot output
#-----------------------------
par(mfrow = c(2, 2))
image(out, main = c("U","V"), ylab = "distance",
  grid = Grid$x.mid, add.contour = TRUE, mfrow = NULL)
par(mar=c(1,1,1,1))
image(out, main = c("U","V"), ylab = "distance", method = "persp",
  grid = Grid$x.mid, border = NA, theta = 30, mfrow = NULL)

mtext(side = 3, outer = TRUE, line = -1.2,"1-D Brusselator",
      cex = 1.5, font = 2)

