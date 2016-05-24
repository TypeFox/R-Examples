## =============================================================================
##
## The medical AKZO problem 
##                     
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 400
##
## =============================================================================
#
# use ReacTran to solve this..
require(ReacTran)
require(deTestSet)

# -------------------------------------------------------
# problem formulation
# -------------------------------------------------------

# The grid properties
Len <- 10 
dx  <- 0.05
x   <- seq (from = dx/2, to = Len, by = dx)
N   <- length(x)

# initial condition of state variables
yini <- c(rep(0,N),rep(1,N))

# model parameters
k   <- 100      # reaction rate
D   <- 1        # diffusion coefficient

# derivative function
medAKZO <- function(t,y,parms) {

    # two state variable vectors
    u <- y[1:N]
    v <- y[-(1:N)]
    
    # boundary condition for u  
    phi <- if (t <= 5) 2 else 0

    # rate of change: only u is transported
    du <- tran.1D(C = u, D = D, C.up = phi, dx = dx)$dC - k*u*v
    dv <- - k*u*v

    # return the rate of changes, concatenated, as a list
    list(c(du, dv))
}

# -------------------------------------------------------
# solve the model
# -------------------------------------------------------

# time sequence
times <- seq(from = 0, to = 20, by = 0.01)

print(system.time(
  out  <- ode.1D(func = medAKZO, times = times, parms = NULL,
                 names = c("u", "v"),  y = yini, nspec = 2)
))
print(system.time(
  out2  <- ode.1D(func = medAKZO, times = times, parms = NULL,
                 y = yini, nspec = 2, method = mebdfi,  restructure = TRUE)
))
### TAKES TOOO LONG
#print(system.time(
#  out3  <- ode.1D(func = medAKZO, times = times, parms = NULL,
#                 y = yini, nspec = 2, method = gamd,  restructure = TRUE)
#))

# -------------------------------------------------------
# plot
# -------------------------------------------------------

# remove time column
Out <- out[,-1]
image( out, which = c("u", "v"), zlim=c(0,2), mfrow = c(2,2))
emptyplot()
colorlegend(zlim = c(0, 2), digit = 2, posx = c(0.5, 0.53))

plot(out, which = 1, lwd = 2, mfrow = NULL, col = "darkblue", main = "u")
for ( i in 1:10) lines(times, out[,10*i], col = "lightblue")
mtext(outer = TRUE,side = 3,"Chemical AKZO",line = -1.5, cex = 1.5)
