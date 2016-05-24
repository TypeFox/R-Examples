odesolv <- function(bwtlist, ystart=diag(rep(1,norder)),
                    h0=width/100, hmin=width*1e-10, hmax=width*0.5,
                    EPS=1e-4, MAXSTP=1000)
{
#  Solve L u = 0,
#  L being an order M homogeneous linear differential operator,
#     (Lu)(t) = w_1(t) u(t) + w_2(t) Du(t) + ...
#                  w_m(t) D^{m-1} u(t) + D^m u(t) = 0
#  for function u and its derivatives up to order m - 1.
#  Each such solution is determined by the values of u and its
#  m - 1 derivatives at time t = 0.  These initial conditions are
#  contained in the columns of the matrix YSTART, which has exactly
#  m rows.  The number of solutions computed by ODESOLV is equal to
#  the number of columnsof YSTART.  In order for the solutions to be
#  linearly independent functions, the columns of YSTART must be
#  linearly independent.  This means that the maximum number of
#  linearly independent solutions is m.
#  The solution for each value of t is a matrix, y(t).
#  Any column of y(t) contains u, Du, ... , D^{m-1} u at argument
#  value t, for the corresponding set of starting values for these
#  first m derivatives. It is the job of this function to estimate
#  these values, and ODESOLV will choose a set of values of TNOW at
#  which these can be estimated with satisfactory accuracy.
#  ODESOLV uses the Runge-Kutta method, which is a good general
#  purpose numerical method.  But it does not work well for stiff
#  systems, and it can fail for poor choices of initial conditions
#  as well as other problems.

#  Arguments:
#  BWTLIST ... a list containing m = 1 weightfunctions.
#             The weight functions w_1, ... , w_m are the functions with
#             indices 2, 3, ..., m+1.
#  YSTART ... initial values for Y.  This is a matrix with M rows,
#             were M is the order of the operator L.  Any column of M
#             specifies intial values for derivatives 0, 1, ... M-1.  Each
#             column must specify a unique set of initial conditions.  A
#             frequent choice is the identity matrix of order M, and this is
#             the default.
#  H0     ... initial trial step size
#  HMIN   ... minimum step size
#  HMAX   ... maximum step size
#  EPS    ... error tolerance
#  MAXSTP ... maximum number of Runge Kutta steps permitted. If the equation
#             is difficult to solve, this may have to be increased.

#  Returns:
#  TP     ... vector of TRUE values used
#  YP     ... m by m by length(TP) array of Y-values generated for
#             values in TP.

#  Note that ODESOLV calls function DERIVS in order to evaluate the
#  differential operator.  Also, it works by redefining the order m
#  linear differential equation as a linear system of m first order
#  differential equations.  DERIVS evaluates the right side of this
#  linear system.

#  Last modified 30 October 2005

  	MAXWARN <- 10

#  determine the order of the system m
  	
   norder <- length(bwtlist)

#  determine the range of values over which the equation is solved

   bfdPar1  <- bwtlist[[1]]
   wbasis   <- bfdPar1$fd$basis
  	rangeval <- wbasis$rangeval
  	tbeg     <- rangeval[1]
  	tend     <- rangeval[2]
  	width    <- tend - tbeg
  	tnow     <- tbeg

  	h        <- min(c(h0, hmax))
  	tp       <- tnow

#  set up the starting values

  	ystartd  <- dim(ystart)
  	if (ystartd[1] != norder) stop("YSTART has incorrect dimensions")
  	n        <- ystartd[2]
  	yp       <- c(ystart)
  	index    <- abs(yp) > 1e-10
  	y        <- ystart

#  Since ODESOLVE is slow, progress is displayed

  	cat("Progress:  each dot is 10% of interval\n")

#  initialize the solution

  	tdel  <- (tend - tbeg)/10
  	tdot  <- tdel
  	iwarn <- 0

#  solve the equation using a maximum of MAXSTP steps

  	for (nstp in 1:MAXSTP) {
		#  evaluate the right side at the current value
    	dydt  <- derivs(tnow, y, bwtlist)
    	yscal <- c(abs(y) + abs(h*dydt) + 1e-30)[index]
    	if (nstp > 1) {
      		tp <- c(tp,tnow)
      		yp <- c(yp,c(y))
    	}
    	if (tnow >= tdot) {
      		cat(".")
      		tdot <- tdot + tdel
    	}
    	if ((tnow+h-tend)*(tnow+h-tbeg) > 0) h <- tend-tnow
		#  take a Runge-Kutta step to the next value
    	result <-  rkqs(y, dydt, tnow, h, bwtlist, yscal, index, EPS)
    	tnow   <- result[[1]]
    	y      <- result[[2]]
    	h      <- result[[3]]
    	hnext  <- result[[4]]
		#  test to see if interval has been covered, and return
		#  if it has.
    	if ((tnow-tend)*(tend-tbeg) >= 0) {
      		cat(".")
      		tp <- c(tp,tnow)
      		yp <- c(yp,c(y))
      		yp <- array(yp,c(norder,n,length(tp)))
      		return(list(tp, yp))
    	}
		#  test if the step is too small.
    	if (abs(hnext) < hmin) {
			warning("Stepsize smaller than minimum")
			hnext <- hmin
			iwarn <- iwarn + 1
			if (iwarn >= MAXWARN) stop("Too many warnings.")
		}
		#  test if the step is too large.
    	h <- min(c(hnext, hmax))
  	}
  	warning("Too many steps.")
}

#  ---------------------------------------------------------------

rkqs <- function(y, dydt, tnow, htry, bwtlist, yscal, index, EPS)
{
#  Take a single step using the Runge-Kutta procedure to
#  control for accuracy.  The function returns the new
#  value of t, the value of the solution at t, the step
#  size, and a proposal for the next step size.
	h <- htry
	#  check the accuracy of the step
   	result <- rkck(y, dydt, tnow, h, bwtlist)
   	ytemp  <- result[[1]]
   	yerr   <- c(result[[2]])[index]
   	errmax <- max(abs(yerr/yscal))/EPS
	#  modify the step size if ERRMAX is too large
   	while (errmax > 1) {
        	htemp <- 0.9*h*(errmax^(-0.25))
        	h     <- max(c(abs(htemp),0.1*abs(h)))
        	tnew  <- tnow + h
        	if (tnew == tnow) stop("stepsize underflow in rkqs")
        	result <- rkck(y, dydt, tnow, h, bwtlist)
        	ytemp  <- result[[1]]
        	yerr   <- result[[2]]
        	errmax <- max(abs(yerr/yscal))/EPS
 	}
	# set up the proposed next step
   	if (errmax > 1.89e-4) {
        	hnext <- 0.9*h*(errmax^(-0.2))
   	} else {
        	hnext <- 5.*h
  	}
   	tnow <- tnow + h
   	y <- ytemp
  	return( list(tnow, y, h, hnext) )
}

#  ------------------------------------------------------------------

rkck <- function(y, dydt, tnow, h, bwtlist)
{
#  Take a single Runge-Kutta step.
#  Return the solution, and an estimate of its error.
      C1  <-   37/378
      C3  <-  250/621
      C4  <-  125/594
      C6  <-  512/1771
      DC5 <- -277/14336
      DC1 <-  C1 - 2825/27648
      DC3 <-  C3 - 18575/48384
      DC4 <-  C4 - 13525/55296
      DC6 <-  C6 - 0.25

      ytemp <- y + h*0.2*dydt
      ak2   <- derivs(tnow+0.2*h,   ytemp, bwtlist)
      ytemp <- y + h*(0.075*dydt  + 0.225*ak2)
      ak3   <- derivs(tnow+0.3*h,   ytemp, bwtlist)
      ytemp <- y + h*(0.3*dydt    - 0.9*ak2 + 1.2*ak3)
      ak4   <- derivs(tnow+0.6*h,   ytemp, bwtlist)
      ytemp <- y + h*(-11*dydt/54 + 2.5*ak2 + (-70*ak3+35*ak4)/27)
      ak5   <- derivs(tnow+h,       ytemp, bwtlist)
      ytemp <- y + h*(1631*dydt/55296  + 575*ak2/512 + 575*ak3/13824 +
                      44275*ak4/110592 + 253*ak5/4096)
      ak6   <- derivs(tnow+0.875*h, ytemp, bwtlist)
      yout  <- y + h*(C1*dydt + C3*ak3 + C4*ak4 + C6*ak6)
      yerr  <- h*(DC1*dydt    + DC3*ak3 + DC4*ak4 + DC5*ak5 + DC6*ak6)
      return (list( yout, yerr ) )
}
