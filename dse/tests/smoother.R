
require("dse")
Sys.info()
DSEversion()
 
fuzz <- 1e-6
digits <- 18
all.ok <- TRUE  

test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

###################################################

# test with input and having output dim < state dim.

###################################################

if(is.R()) data("eg1.DSE.data.diff", package="dse")
model <- TSmodel(toSSChol(estVARXls(eg1.DSE.data.diff))) 

model0 <- model
model0$G <- NULL
simdata0 <- simulate(model0, rng=test.rng) 

z  <- smoother(model0, simdata0, compiled=TRUE)
zz <- smoother(model0, simdata0, compiled=FALSE)

#tfplot(simdata0$state, state(z, smooth=TRUE), state(zz, smooth=TRUE), graphs.per.page=3)


# using simulated data gives a true state for comparison.
simdata <- simulate(model, input= inputData(eg1.DSE.data.diff), rng=test.rng) 

z  <- smoother(model, simdata, compiled=TRUE)
zz <- smoother(model, simdata, compiled=FALSE)


error <- max(abs((state(zz, smooth=TRUE) - zz$smooth$state)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

tfplot(state(zz), simdata$state, graphs.per.page=3)
#tfplot(state(z, smoother=TRUE) - state(zz, smoother=TRUE), graphs.per.page=3)
tfplot(state(z, smoother=TRUE),  state(zz, smoother=TRUE), graphs.per.page=3)

# plot smoother agains true state
#tfplot(state(z, smoother=TRUE), simdata$state, graphs.per.page=3)

# plot smoother agains true state
#tfplot(state(zz, smoother=TRUE), simdata$state, graphs.per.page=3)

#tfplot(state(z, smoother=TRUE), simdata$state, graphs.per.page=3)
#tfplot(simdata$state, state(z, smoother=TRUE), state(zz, smoother=TRUE), graphs.per.page=3)

#tfplot(simdata$state, state(zz, filter=TRUE), state(zz, smoother=TRUE), graphs.per.page=3)

# compare fortran and S versions
error <- max(abs((state(z, smoother=TRUE) - state(zz, smoother=TRUE))))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(z$smooth$track - zz$smooth$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(z$filter$track - zz$filter$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }


######################################

# test output dim exceeds state dim.

######################################

Hloadings <- t(matrix(c(
    8.8,   5.2,
   23.8, -12.6,
    5.2,  -2.0,
   36.8,  16.9,
   -2.8,  31.0,
    2.6,  47.6), 2,6))

ss.ar1 <- SS(F=array(c(.5, .4, .3, .2),c(2,2)),  
		H=Hloadings,     
		Q=array(c(1.0, 2.0),c(2,2)),  
		R=diag(1,6) 
		)

simdata2 <- simulate(ss.ar1, rng=test.rng)

z  <- smoother(ss.ar1, simdata2, compiled = TRUE)
zz <- smoother(ss.ar1, simdata2, compiled = FALSE)

#tfplot(state(zz, smooth=TRUE), state(z, smooth=TRUE), simdata2$state, graphs.per.page=3)
#tfplot(state(zz, smooth=TRUE) - state(z, smooth=TRUE), graphs.per.page=3)

error <- max(abs((state(zz, filter=TRUE) - zz$filter$state)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }


error <- max(abs((state(zz, filter=TRUE) - state(z, filter=TRUE))))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs((state(zz, smooth=TRUE) - zz$smooth$state)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs((state(z, smoother=TRUE) - state(zz, smoother=TRUE))))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }


error <- max(abs((z$smooth$track) - zz$smooth$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }


tfplot(simdata2$state, state(zz, smoother=TRUE), state(zz, filter=TRUE))
tfplot(simdata2$state, state(z,  smoother=TRUE), state(z,  filter=TRUE))
tfplot(simdata2$state, state(z,  smoother=TRUE), state(zz, smoother=TRUE))

######################################

# test "big k" (which is numerically sensitive).

######################################

#  Starting P0  ("big k") symmetric with off diagonal element smaller than diag.
P0 <- matrix(1e6,4,4) 
diag(P0 )<- 1e7

mod4 <-  SS(F=t(matrix(c(
    		  0.8, 0.04,  0.2, 0,
    		  0.2,  0.5,    0, -0.3,
    		    1,    0,    0, -0.2,
    		    0,	  1,    0,  0   ), c(4,4))),
	       H=cbind(Hloadings, matrix(0,6,2)),	
	       Q=diag(c(1, 1, 0, 0),4),  
	       R=diag(1,6),
	       z0=c(10, 20, 30,40),
	       P0=P0
	       )

z  <- simulate(SS(F=t(matrix(c(
    		  0.8, 0.04,  0.2, 0,
    		  0.2,  0.5,    0, -0.3,
    		    1,    0,    0, -0.2,
    		    0,	  1,    0,  0   ), c(4,4))),
	       H=cbind(Hloadings, matrix(0,6,2)),	
	       Q=diag(c(1, 1, 0, 0),4),  
	       R=diag(1,6),
	       z0=c(10, 20, 30,40),
	       P0=diag(c(10, 10, 10, 10)) ),
               rng=test.rng)  

state.sim  <- z$state  # for comparison below
y.sim  <- outputData(z) # simulated indicators


error <- max(abs(l(mod4, TSdata(output=y.sim), return.state=TRUE)$filter$state -
                 l(mod4, TSdata(output=y.sim), return.state=TRUE,
	                                           compile=FALSE)$filter$state))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }


zz  <- smoother(l(mod4, TSdata(output=y.sim)))   
zzz <- smoother(l(mod4, TSdata(output=y.sim)), compiled=FALSE)

tfplot(state.sim,  state(zz))
tfplot(state.sim,  state(zzz))
tfplot(state.sim,  state(zz,  smoother=TRUE))
tfplot(state.sim,  state(zzz, smoother=TRUE))

error <- max(abs(state(zz, filter=TRUE) - state(zzz, filter=TRUE)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(state(zz, smoother=TRUE) - state(zzz, smoother=TRUE)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(zz$filter$track - zzz$filter$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(zz$smooth$track - zzz$smooth$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }



if (! all.ok) stop("some tests FAILED")

