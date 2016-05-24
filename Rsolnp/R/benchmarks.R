#################################################################################
##
##   R package Rsolnp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009-2013
##   This file is part of the R package Rsolnp.
##
##   The R package Rsolnp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package Rsolnp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

benchmarkids <- function()
{
return(c("Powell", "Wright4", "Wright9", "Alkylation", "Entropy", "Box", "RosenSuzuki",
				"Electron", "Permutation"))
}


benchmark <- function( id = "Powell")
{
  if( !any(benchmarkids() == id[ 1L ]) )
    stop( "invalid benchmark id" )
 	ans = switch(id,
			Powell = .powell(),
			Wright4 = .wright4(),
			Wright9 = .wright9(),
			Alkylation = .alkylation(),
			Entropy = .entropy(),
			Box = .box(),
			RosenSuzuki = .rosensuzuki(),
			Electron = .electron(),
			Permutation = .permutation())
	return(ans)
}

.powell = function()
{
	.fn1 = function(x)
	{
		exp(x[1]*x[2]*x[3]*x[4]*x[5])
	}

	.eqn1 = function(x){
		z1=x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]
		z2=x[2]*x[3]-5*x[4]*x[5]
		z3=x[1]*x[1]*x[1]+x[2]*x[2]*x[2]
		return(c(z1,z2,z3))
	}

	.x0 = c(-2, 2, 2, -1, -1)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = c(10,0,-1), control=ctrl)
	minos = list()
	minos$fn = 0.05394985
	minos$pars = c(-1.717144, 1.595710, 1.827245, 0.763643, 0.763643)
	minos$nfun = 524
	minos$iter = 12
	minos$elapsed = 0.2184

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("Powell's exponential problem is a function of five variables with
three nonlinear equality constraints on the variables.")

	return(bt)
}

.wright4 = function()
{
	.fn1 = function(x)
	{
		(x[1]-1)^2+(x[1]-x[2])^2+(x[2]-x[3])^3+(x[3]-x[4])^4+(x[4]-x[5])^4
	}

	.eqn1 = function(x){
		z1=x[1]+x[2]*x[2]+x[3]*x[3]*x[3]
		z2=x[2]-x[3]*x[3]+x[4]
		z3=x[1]*x[5]
		return(c(z1,z2,z3))
	}

	.x0 = c(1, 1, 1, 1, 1)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = c(2+3*sqrt(2),-2+2*sqrt(2),2), control=ctrl)
	minos = list()
	minos$fn = 0.02931083
	minos$pars = c(1.116635, 1.220442, 1.537785, 1.972769, 1.791096)
	minos$nfun = 560
	minos$iter = 9
	minos$elapsed = 0.249

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("Wright's fourth problem is a function of five variables with
three non linear equality constraints on the variables. This popular
test problem has several local solutions and taken from Wright (1976).")
	return(bt)
}

.wright9 = function()
{
	.fn1 = function(x)
	{
		10*x[1]*x[4]-6*x[3]*x[2]*x[2]+x[2]*(x[1]*x[1]*x[1])+
				9*sin(x[5]-x[3])+x[5]^4*x[4]*x[4]*x[2]*x[2]*x[2]
	}

	.ineqn1 = function(x){
		z1=x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]
		z2=x[1]*x[1]*x[3]-x[4]*x[5]
		z3=x[2]*x[2]*x[4]+10*x[1]*x[5]
		return(c(z1,z2,z3))
	}
	ineqLB = c(-100, -2, 5)
	ineqUB = c(20, 100, 100)
	.x0 = c(1, 1, 1, 1, 1)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, ineqfun = .ineqn1, ineqLB = ineqLB, ineqUB = ineqUB, control=ctrl)
	minos = list()
	minos$fn = -210.4078
	minos$pars = c(-0.08145219, 3.69237756, 2.48741102,  0.37713392, 0.17398257)
	minos$nfun = 794
	minos$iter = 11
	minos$elapsed = 0.281

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("Wright's ninth problem is a function of five variables with three
non linear inequality constraints on the variables. This popular test
problem has several local solutions and taken from Wright (1976).")
	return(bt)
}

.alkylation = function()
{
	.fn1 = function(x)
	{
		-0.63*x[4]*x[7]+50.4*x[1]+3.5*x[2]+x[3]+33.6*x[5]
	}

	.eqn1 = function(x){
		z1=98*x[3]-0.1*x[4]*x[6]*x[9]-x[3]*x[6]
		z2=1000*x[2]+100*x[5]-100*x[1]*x[8]
		z3=122*x[4]-100*x[1]-100*x[5]
		return(c(z1,z2,z3))
	}
	.ineqn1 = function(x){
		z1=(1.12*x[1]+0.13167*x[1]*x[8]-0.00667*x[1]*x[8]*x[8])/x[4]
		z2=(1.098*x[8]-0.038*x[8]*x[8]+0.325*x[6]+57.25)/x[7]
		z3=(-0.222*x[10]+35.82)/x[9]
		z4=(3*x[7]-133)/x[10]
		return(c(z1,z2,z3,z4))
	}
	ineqLB = c(0.99,0.99,0.9,0.99)
	ineqUB = c(100/99,100/99,10/9,100/99)
	eqB = c(0,0,0)
	LB = c(0,0,0,10,0,85,10,3,1,145)
	UB = c(20,16,120,50,20,93,95,12,4,162)
	.x0 = c(17.45,12,110,30,19.74,89.2,92.8,8,3.6,155)
	ctrl = list(rho = 0, trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = eqB, ineqfun = .ineqn1, ineqLB = ineqLB,
			ineqUB = ineqUB, LB = LB, UB = UB, control = ctrl)
	minos = list()
	minos$fn = -172.642
	minos$pars = c(16.996427, 16.000000, 57.685751, 30.324940, 20.000000, 90.565147, 95.000000, 10.590461, 1.561636, 153.535354)
	minos$nfun = 2587
	minos$iter = 13
	minos$elapsed = 0.811

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("The Alkylation problem models a simplified alkylation process. It
is a function of ten variables with four non linear inequality and
three non linear equality constraints as well as variable bounds.
The problem is taken from Locke and Westerberg (1980).")
	return(bt)
}

.entropy = function()
{
	.fn1 = function(x)
	{
		m = length(x)
		f = 0
		for(i in 1:m){
			f = f-log(x[i])
		}
		ans = f-log(.vnorm(x-1) + 0.1)
		ans
	}

	.eqn1 = function(x){
		sum(x)
	}
	eqB = 10
	LB = rep(0,10)
	UB = rep(1000,10)
	.x0 = runif(10, 0, 1000)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = eqB, LB = LB, UB = UB, control=ctrl)
	minos = list()
	minos$fn = 0.1854782
	minos$pars = c(2.2801555, 0.8577605, 0.8577605, 0.8577605, 0.8577605, 0.8577605,
			0.8577605, 0.8577605, 0.8577605, 0.8577605)
	minos$nfun = 886
	minos$iter = 4
	minos$elapsed = 0.296

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("The Entropy problem is non convex in n variables with one linear
equality constraint and variable positivity bounds.")
	return(bt)
}

.box = function()
{
	.fn1 = function(x)
	{
		-x[1]*x[2]*x[3]
	}

	.eqn1 = function(x){
		4*x[1]*x[2]+2*x[2]*x[3]+2*x[3]*x[1]
	}

	eqB = 100
	LB = rep(1, 3)
	UB = rep(10, 3)

	.x0 = c(1.1, 1.1, 9)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = eqB, LB = LB, UB = UB, control=ctrl)
	minos = list()
	minos$fn = -48.11252
	minos$pars = c(2.886751, 2.886751, 5.773503)
	minos$nfun = 394
	minos$iter = 9
	minos$elapsed = 0.156

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("The box problem is a function of three variables with one non
linear equality constraint and variable bounds.")
	return(bt)
}

.rosensuzuki = function()
{
	.fn1 = function(x)
	{
		x[1]*x[1]+x[2]*x[2]+2*x[3]*x[3]+x[4]*x[4]-5*x[1]-5*x[2]-21*x[3]+7*x[4]
	}

	.ineqn1 = function(x){
		z1=8-x[1]*x[1]-x[2]*x[2]-x[3]*x[3]-x[4]*x[4]-x[1]+x[2]-x[3]+x[4]
		z2=10-x[1]*x[1]-2*x[2]*x[2]-x[3]*x[3]-2*x[4]*x[4]+x[1]+x[4]
		z3=5-2*x[1]*x[1]-x[2]*x[2]-x[3]*x[3]-2*x[1]+x[2]+x[4]
		return(c(z1,z2,z3))
	}
	ineqLB = rep(0, 3)
	ineqUB = rep(1000, 3)
	.x0 = c(1, 1, 1, 1)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, ineqfun = .ineqn1, ineqLB = ineqLB, ineqUB = ineqUB, control=ctrl)
	minos = list()
	minos$fn = -44
	minos$pars = c(2.502771e-07, 9.999997e-01, 2.000000e+00, -1.000000e+00)
	minos$nfun = 527
	minos$iter = 12
	minos$elapsed = 0.203

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("The Rosen-Suzuki problem is a function of four variables with
three nonlinear inequality constraints on the variables. It is taken
from Problem 43 of Hock and Schittkowski (1981).")
	return(bt)
}



#----------------------------------------------------------------------------------
# Some Problems in Global Optimization
#----------------------------------------------------------------------------------


# Distribution of Electrons on a Sphere
# Given n electrons, find the equilibrium state distribution (of minimal Coulomb potential)
# of the electrons positioned on a conducting sphere. This model is from the COPS benchmarking suite.
# See http://www-unix.mcs.anl.gov/~more/cops/.

.electron = function()
{
	gofn = function(dat, n)
	{

		x = dat[1:n]
		y = dat[(n+1):(2*n)]
		z = dat[(2*n+1):(3*n)]
		ii = matrix(1:n, ncol = n, nrow = n, byrow = TRUE)
		jj = matrix(1:n, ncol = n, nrow = n)
		ij = which(ii<jj, arr.ind = TRUE)
		i = ij[,1]
		j = ij[,2]
		#  Coulomb potential
		potential = sum(1.0/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2))
		potential
	}

	goeqfn = function(dat, n)
	{
		x = dat[1:n]
		y = dat[(n+1):(2*n)]
		z = dat[(2*n+1):(3*n)]
		apply(cbind(x^2, y^2, z^2), 1, "sum")
	}

	n = 25
	LB = rep(-1, 3*n)
	UB = rep(1, 3*n)
	eqB = rep(1, n)
	ans = gosolnp(pars  = NULL, fixed = NULL, fun = gofn, eqfun = goeqfn, eqB = eqB, LB = LB, UB = UB,
			control = list(), distr = rep(1, length(LB)), distr.opt = list(outer.iter = 10, trace = 1),
			n.restarts = 2, n.sim = 20000, rseed = 443, n = 25)

	conopt = list()
	conopt$fn  = 243.813
	conopt$iter = 33
	conopt$nfun = NA
	conopt$elapsed = 0.041
	conopt$pars = c(-0.0117133872042326,	0.627138691757704,	-0.471025867741051,	-0.164419761338935,	-0.0315460712487934,
			-0.12718981058582,	-0.540049624346613,	0.600346770449059,	0.29796281847713,	-0.740960572770077,	0.972512148478245,
			-0.870858895858346,	0.84178885636396,	-0.182471994739506,	0.603293664844919,	0.0834172554171806,	0.51317309921937,
			0.260639996237799,	-0.0972877803105543,	-0.979882559381314,	-0.64809471648373,	-0.722351411610064,	0.847184430059647,
			0.514683899757428,	-0.574607272207711, 0.114609211815613,	-0.748168886860133,	-0.379763612890494,	0.743271243797936,
			0.846784034469756,	0.220425955966718,	0.839147591392778,	-0.613810163641104,	-0.499794531840362,	-0.199680552951248,
			-0.105141937435843,	-0.434753057357539,	-0.127562956191463,	-0.895691740627038,	0.574257349984438,	-0.967631920158332,
			-0.0243647398873149,	0.959445715727407,	-0.517406241891138,	0.197677956191858,	0.503867654081605,	0.286619450971711,
			0.522509289901788,	0.474911361600545,	-0.768915699421274, -0.993341595387613,	-0.21665728244143,	-0.796187308516752,	-0.648470508368975,
			0.531000606737778,	0.967075565826839,	-0.0646353084836963,	0.512660548728168,	-0.813279524362711,	0.641174786133487,
			0.207762590603952,	-0.229335044471304,	-0.524518077390237,	-0.405512363446903,	0.55341236880543,	0.23818066376044,
			0.857939234263016,	-0.107381147942665,	0.850191665834437,	0.0274516931378701,	0.571043453369507,	-0.629331175510656,
			-0.0962423162171412,	-0.713834491989006,	0.280348229724225)
	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			conopt =  rbind(round(conopt$fn, 5L),
					round(conopt$iter, 0L),
					round(0, 0L),
					round(conopt$nfun, 0L),
					round(conopt$elapsed, 3L),
					matrix(round(conopt$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	colnames(bt) = c("solnp", "conopt")
	attr(bt, "description") = paste("The equilibrium state distribution (of minimal Coulomb potential)\n of the electrons positioned on a conducting sphere.")
	return(bt)
}

# Permutation Problem -- Unique Solution f(x) = 0 and x(i) = i

.permutation = function()
{
	.perm = function(x, n, b){
		F = 0
		for(k in 1:n){
			S = 0
			for(i in 1:n){
				S = S + ( ( (i^k) + b ) * (( x[i]/i )^k -1))
			}
			F = F + S^2
		}
		F
	}

	ans = gosolnp(pars  = NULL, fixed = NULL, fun = .perm, eqfun = NULL, eqB = NULL, LB = rep(-4, 4), UB = rep(4, 4),
			control = list(outer.iter = 25, trace = 1, tol = 1e-9), distr = rep(1, 4), distr.opt = list(),
			n.restarts = 6, n.sim = 20000, rseed = 99, n = 4, b =0.5)


	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			actual =  rbind(0,
					NA,
					round(0, 0L),
					NA,
					NA,
					matrix(1:4, ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	colnames(bt) = c("solnp", "expected")
	attr(bt, "description") = paste("Permutation Problem PERM(4,0.5).")

}