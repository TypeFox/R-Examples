########################################################################
## DEMO: ecological model applications have common properties,
##       that belong together:
##       - equations
##       - constant parameters
##       - initial values of the state variables
##       - time steps
##       - an adequate solver, and
##       - utility functions
########################################################################

# if (dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))

if (require(odesolve)) {

	library(odesolve)
	library(proto)

	## object creation from scratch (ex-nihilo)
	lv <- proto(expr = {
	    equations <- function(t, x, p) {
	    dx1.dt <-   p["k1"] * x[1] - p["k2"] * x[1] * x[2]
	    dx2.dt <- - p["k3"] * x[2] + p["k2"] * x[1] * x[2]
	    list(c(dx1.dt, dx2.dt))
	    }
	    # vectors of parameters, timesteps and initial values
	    parms  <- c(k1=0.2, k2=0.2, k3=0.2)
	    times  <- 1:100
	    init   <- c(prey=0.5, predator=1)
	    # two methods
	    solve  <- function(.) {
		  # must use .$with(equations) instead of .$equations
		  # as latter can only be used in a call
		  equations <- .$with(equations)
		  res <- as.data.frame(rk4(.$init, .$times, equations, .$parms))
		}
	    plot <- function(.) {
	      res <- .$solve()
	      
	      graphics::plot(res$time, res$predator, col = "red", 
		 type = "b", pch = 20)
	      graphics::plot(res$time, res$prey, col = "green",
		 type = "b", pch = 20)
	    }
	})    

	## the created object is fully functional
	par(mfrow=c(2,1))
	res <- lv$solve() 
	plot(res)
	lv$plot()  

	## derive a child object with a different ODE solver
	lv.lsoda <- lv$proto(
		      solve <- function(.) {
			equations <- .$with(equations)
			res <- lsoda(.$init, .$times, equations, .$parms)
			as.data.frame(res)
		      }
		    )

	## test this child object
	lv.lsoda$plot()

	### derive two scenarios

	## scenario 1 with initial values, that are in equilibrium
	sc1 <- lv.lsoda$proto(init   <- c(prey=1, predator=1))

	## scenario 2, similar to scenario 1, but with different parameters
	sc2 <- sc1$proto(parms  <- c(k1=0.3, k2=0.2, k3=0.2))

	## solve, plot, compare these scenarios
	par(mfrow=c(2,2))
	sc1$plot()
	sc2$plot()
			      
	## and as a last example, we modify parameters of the parent
	lv$times <- seq(0, 100, 0.1)

	sc1$main <- "Scenario 1"
	sc2$main <- "Scenario 2"
	lv$plot <- function(.) {
	      res <- .$solve()
	      
	      graphics::plot(res$time, res$predator, col = "red", 
		 type = "l", pch = 20, main=.$main)
	      graphics::plot(res$time, res$prey, col = "green",
		 type = "l", pch = 20, main=.$main)
	    }

	par(mfrow=c(2,2))

	sc1$plot()
	sc2$plot()

	## fragile base object: the folowing would give an error
	# lv$plot()

	## but it works if we add the required main title
	lv$main <- "Reference"
	lv$plot()

	## show the object relationships
	g <- graph.proto()
	plot(g)

	## Conclusion:
	##   - the responsibility is due to the user, but
	##   - prototypes allow oop without overhead

}

par(opar)

