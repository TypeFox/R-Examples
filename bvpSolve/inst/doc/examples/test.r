## ================================================================     comments
## fluid injection problem
## Original definition:
## f'''- R[(f')^2 -f*f''] + A = 0
## h'' + R*f*h' + 1 = 0
## O'' + P*f*O' = 0
## A is unknown, P = 0.7*R
##
## rewritten as:
## df=f1                   #f
## df1=f2                  #f1=f'
## df2=R(f1^2-f*f2)-A      #f2=f''
## dh=h1
## dh1= -Rfh1-1
## dO=O1
## dO1 = O2
## dO2 = -P*f*O1
## dA = 0                  # the constant to be estimated
## ================================================================ end comments


# load the package with the solver
require(bvpSolve)

# the derivative function
fluid<-function(t, y, parms, R) {
 P    = 0.7*R
 with(as.list(y),  {
   df  <- f1
   df1 <- f2
   df2 <- R*(f1^2-f*f2)-A
   dh  <- h1
   dh1 <- -R*f*h1-1
   dO  <- O1
   dO1 <- -P*f*O1
   dA  <- 0

   return(list( c(df,df1,df2,dh,dh1,dO,dO1,dA)) )
 })
}

# fixed points
times  <- seq(0, 1, by = 0.1)      # -> fixpnt

# boundary values
yini   <- c(f=0, f1=0, f2=NA, h=0, h1=NA, O=0, O1=NA, A=NA)
yend   <- c(f=1, f1=0, f2=NA, h=0, h1=NA, O=1, O1=NA, A=NA)

# Solving the model, using twpbvpC this works
Soltwp <- bvptwp(func = fluid, x = times, parms = NULL, R = 200,
                nmax = 195, cond = TRUE, #verbose=TRUE,
                yini = yini, yend = yend)
diagnostics(Soltwp)

# Solving the model, using twpbvpC + very large value of R + conditioning
Soltwp2 <- bvptwp(func = fluid, x = times, parms = NULL, R=5000,
                nmax = 1100, cond = TRUE, verbose = TRUE,    #
                yini = yini, yend = yend, xguess = Soltwp[,1],
                yguess = t(Soltwp[,-1]))
diagnostics(Soltwp2)

# plot the results
plot(Soltwp2, type = "l", lwd = 2)

