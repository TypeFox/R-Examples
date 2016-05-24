# The ISOpureR package is copyright (c) 2014 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### FUNCTION: ISOpure.model_optimize.cg_code.rminimize.R ########################################################################
#
# Minimize a differentiable multivariate function.
# 
# This function is a conjugate-gradient search with interpolation/extrapolation
# by Carl Edward Rasmussen.  A description of the Matlab code can be found at
# http://learning.eng.cam.ac.uk/carl/code/minimize/ (accessed Jan. 21, 2014).
# This is a implementation in R.  
# 
# See also the explanantion of the interpolation/extrapolation in the pdf document
# cd_code_details.pdf in the ISOpureR directory. 
#
# Function call: ISOpure.model_optimize.cg_code.rminimize <- function(X, f, df, run_length, ...)

### INPUT #########################################################################################
# X: the starting point is given by "X" (D by 1)
#    X must be either a scalar or a column vector/matrix, not a row matrix
#
# f: the name of the function to be minimized, returning a scalar
#
# df: the name of the function which returns the vector of partial derivatives of f wrt X,
#     where again the partial derivatives must be in scalar or column vector/matrix form
#
# run_length: gives the length of the run: if it is positive, it gives the maximum number 
#     of line searches, if negative its absolute gives the maximum allowed number of function 
#     evaluations. You can (optionally) give "run_length" a second component, which will indicate the
#     reduction in function value to be expected in the first line-search (defaults
#     to 1.0). 
#
#     REVISIT: this R code was only used with a positive, single component for run_length, so 
#     would need to test this.
#     
# ...: parameters to be passed on to the function f.

### OUTPUT ########################################################################################
# The function return a list with three components:
#    X: the found solution X
#    fX: a vector of function values fX indicating the progress made
#        TO DO: check whether the last entry is the most recent!
#    i: the number of iterations (line searches or function evaluations, depending on the sign 
#       of "run_length") used.

### NOTES #########################################################################################
# 
# The function returns when either its length is up, or if no further progress
# can be made (ie, we are at a (local) minimum, or so close that due to
# numerical problems, we cannot get any closer). NOTE: If the function
# terminates within a few iterations, it could be an indication that the
# function values and derivatives are not consistent (ie, there may be a bug in
# the implementation of your "f" function).
#
# The Polack-Ribiere flavour of conjugate gradients is used to compute search
# directions, and a line search using quadratic and cubic polynomial
# approximations and the Wolfe-Powell stopping criteria is used together with
# the slope ratio method for guessing initial step sizes. Additionally a bunch
# of checks are made to make sure that exploration is taking place and that
# extrapolation will not be unboundedly large.
#
# See also: checkgrad 
#
# Copyright (C) 2001 - 2006 by Carl Edward Rasmussen (2006-09-08).

ISOpure.model_optimize.cg_code.rminimize <- function(X, f, df, run_length, ...) {
	INT <- 0.1;  # don't reevaluate within 0.1 of the limit of the current bracket
	EXT <- 3.0;  # extrapolate maximum 3 times the current step-size
	MAX <- 20;   # max 20 function evaluations per line search
	RATIO <- 10; # maximum allowed slope ratio
	
	# SIG and RHO are the constants controlling the Wolfe-
	# Powell conditions. SIG is the maximum allowed absolute ratio between
	# previous and new slopes (derivatives in the search direction), thus setting
	# SIG to low (positive) values forces higher precision in the line-searches.
	# RHO is the minimum allowed fraction of the expected (from the slope at the
	# initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
	# Tuning of SIG (depending on the nature of the function to be optimized) may
	# speed up the minimization; it is probably not worth playing much with RHO.
	SIG <- 0.1;
	RHO <- SIG/2;

	# The code falls naturally into 3 parts, after the initial line search is
	# started in the direction of steepest descent. 1) we first enter a while loop
	# which uses point 1 (p1) and (p2) to compute an extrapolation (p3), until we
	# have extrapolated far enough (Wolfe-Powell conditions). 2) if necessary, we
	# enter the second loop which takes p2, p3 and p4 chooses the subinterval
	# containing a (local) minimum, and interpolates it, unil an acceptable point
	# is found (Wolfe-Powell conditions). Note, that points are always maintained
	# in order p0 <= p1 <= p2 < p3 < p4. 3) compute a new search direction using
	# conjugate gradients (Polack-Ribiere flavour), or revert to steepest if there
	# was a problem in the previous line-search. Return the best value so far, if
	# two consecutive line-searches fail, or whenever we run out of function
	# evaluations or line-searches. During extrapolation, the "f" function may fail
	# either with an error or returning Nan or Inf, and minimize should handle this
	# gracefully.

	if (length(run_length) == 2) {
		red <- run_length[2];
		run_length <- run_length[1];
		}
	else {
		red = 1;
		}
	
	if (run_length > 0) {
		S = 'Linesearch';
		}
	else{
		S = 'Function evaluation';
		}
	
	# zero the run length counter
	i <- 0;
	# no previous line search has failed
	ls_failed <- 0;
	
	# get function value and gradient
	f0 <- f(X,...);
	df0 <- df(X,...);
	fX <- f0;

	# count epochs?!
	i <- i + (run_length < 0);
	# initial search direction (steepest) and slope
	s <- -df0;
	# conjugate transpose fix by Francis, to match the Matlab multiplication
	d0 <- Conj(t(-s)) %*% s;  
	# initial step is red/(|s|+1)                  
	x3 <- as.numeric(red / (1 - d0)); 
	
	# while not finished
	while (i < abs(run_length)) {
		# count iterations?!
		i <- i + (run_length > 0);
		# make a copy of current values
		X0 <- X;
		F0 <- f0;
		dF0 <- df0;
		if (run_length > 0) {
			M <- MAX;
			}
		else {
			M <- min(c(MAX, - run_length - i));
			}

		# keep extrapolating as long as necessary
		while(TRUE) {
			x2 <- 0;
			f2 <- f0;
			d2 <- d0;
			f3 <- f0;
			df3 <- df0;
			success <- 0;
			while( !(success) && (M > 0)) {
			# count epochs?!
				M <- M - 1;
				i <- i + (run_length < 0);
				x3 <- as.numeric(x3);  
				f3 <- f(as.vector(X + x3 * s), ...); 
				df3 <- df(as.vector(X + x3 * s), ...); 
				# catch any error which occured in f
				if (is.nan(f3) || is.infinite(f3) || any(is.nan(df3)) || any(is.infinite(df3))) {
					# bisect and try again
					x3 <- (x2 + x3)/2; 
					}
				else {
					success <- 1;
					}
				}
			# keep best values   
			if (ISOpure.util.matlab_less_than(f3,F0)){
				X0 <- X + x3 * s;
				F0 <- f3;
				dF0 <- df3;
				}
			# new slope; conjugate transpose fix
			d3 <- Conj(t(df3)) %*% s;                
			if (ISOpure.util.matlab_greater_than(d3, SIG*d0) || ISOpure.util.matlab_greater_than(f3, f0+x3*RHO*d0) || M == 0){
				break;
				}
			# move point 2 to point 1
			x1 <- x2; 
			f1 <- f2; 
			d1 <- d2; 
			# move point 3 to point 2
			x2 <- x3; 
			f2 <- f3; 
			d2 <- d3; 
			# make cubic extrapolation                              
			A <- 6*(f1-f2)+3*(d2+d1)*(x2-x1); 
			B <- 3*(f2-f1)-(2*d1+d2)*(x2-x1); 
			# num. error possible, ok!
			x3 <- x1-d1*(x2-x1)^2/(B+sqrt(B*B-A*d1*(x2-x1))); 
			if (is.nan(x3) || is.infinite(x3) || ISOpure.util.matlab_less_than(x3, 0)) {
				x3 <- x2 * EXT;    
				}
			# new point beyond extrapolation limit?
			else if (ISOpure.util.matlab_greater_than(x3, x2 * EXT)) {
				x3 <- x2 * EXT;
				}
			# new point too close to previous point?
			else if (ISOpure.util.matlab_less_than(x3, x2 + INT*(x2-x1))){
				x3 <- x2 + INT*(x2-x1); 
				}
			}
		# end extrapolation
		
		# keep interpolating 
		while ( (ISOpure.util.matlab_greater_than(abs(d3), -SIG*d0) || ISOpure.util.matlab_greater_than(f3,f0+x3*RHO*d0)) && M > 0){
			# choose subinterval
			if ( ISOpure.util.matlab_greater_than(d3, 0) || ISOpure.util.matlab_greater_than(f3, f0+x3*RHO*d0)){
				# move point 3 to point 4
				x4 <- x3; 
				f4 <- f3;
				d4 <- d3;
				}
			else {
				# move point 3 to point 2
				x2 <- x3; 
				f2 <- f3;
				d2 <- d3;                                
				}
			if (ISOpure.util.matlab_greater_than(f4, f0)) {
				# quadratic interpolation
				x3 <- x2 - (0.5*d2*(x4-x2)^2)/(f4-f2-d2*(x4-x2)); 
				}
			else {
				# cubic interpolation
				A <- 6*(f2-f4)/(x4-x2)+3*(d4+d2); 
				B <- 3*(f4-f2)-(2*d2+d4)*(x4-x2); 
				x3 <- x2+(sqrt(B*B-A*d2*(x4-x2)^2)-B)/A; 
				}
			# if we had a numerical problem then bisect
			if (is.nan(x3) || is.infinite(x3)) {
				x3 <- (x2+x4)/2; 
				}
			# don't accept too close
			x3 <- max(min(x3, x4 - INT * (x4 - x2)), x2 + INT * (x4 - x2)); 
			x3 <- as.numeric(x3); 
			f3 <- f(as.vector(X + x3 * s), ...);
			df3 <- df(as.vector(X + x3 * s), ...);
			# keep best values
			if (ISOpure.util.matlab_less_than(f3, F0)) {  
				X0 <- X + x3 * s;
				F0 <- f3;
				dF0 <- df3;
				}
			# count epochs?!
			M <- M - 1;
			i <- i + (run_length < 0);
			# new slope; conjugate transpose fix 
			d3 <- Conj(t(df3)) %*% s;                    
			}
		# end interpolation
		# if line search succeeded
		if (ISOpure.util.matlab_less_than(abs(d3),-SIG*d0) && ISOpure.util.matlab_less_than(f3, f0+x3*RHO*d0)){ 
			# update variables
			X <- X + x3 * s;
			f0 <- f3;
			# conjugate transpose fix
			fX <- t(c(Conj(t(fX)), f0));            
			m1 <- ((t(df3) %*% df3)-(t(df0) %*% df3));
			m2 <- (t(df0) %*% df0);
			# Polack-Ribiere CG direction
			s <- as.numeric(m1/m2)*s - df3;
			# swap derivatives 
			df0 <- df3;
			d3 <- d0;
			d0 <- Conj(t(df0)) %*% s; 
			# new slope must be negative
			if (ISOpure.util.matlab_greater_than(max(d0), 0)) {
				# otherwise use steepest direction
				s <- -df0;                
				d0 <- Conj(t(-s)) %*% s;
				}
			# slope ratio but max RATIO
			x3 <- x3 * min( RATIO, d3/(d0 - .Machine$double.xmin)); 
			# this line search did not fail
			ls_failed <- 0;
		}
		else { 
			# restore best point so far
			X <- X0;
			f0 <- F0;
			df0 <- dF0;
			# line search failed twice in a row
			# or we ran out of time, so we give up
			if (ls_failed || i > abs(run_length)) {
				break;
				}
			# try steepest    
			s <- -df0;
			d0 <- Conj(t(-s)) %*% s;
			x3 <- 1/(1-d0); 
			# this line search failed
			ls_failed <- 1;
			}
		}
		return(list(as.vector(X),as.vector(fX),i))
	}
