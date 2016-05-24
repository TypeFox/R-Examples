# File:   csem.R
# Author: Zhenghua Nie 
# Date:   Mon 26 Dec 2011 23:54:22 EST
#
#
#
# Copyright (C) 2011 Zhenghua Nie. All Rights Reserved.
# This code is published under GNU GENERAL PUBLIC LICENSE.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License,  or
# (at your option) any later version.
#      
# This program is distributed WITHOUT ANY WARRANTY. See the
# GNU General Public License for more details.
#           
# If you do not have a copy of the GNU General Public License,  
# write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


# The following function is a wrapper to compute the objective function and its gradient.
# If hessian=TRUE,  csem will return Hessian computed by the numerical method,  but
# is flexible to return Hessian computed by the analytical solution.
CompiledObjective <- function(par, model.description, gradient=TRUE, hessian=FALSE, objective=c("objectiveML", "objectiveGLS", "objectiveFIML", "objectivelogLik"), ...)
{
		if(missing(objective)) objective <- "objectiveML"
		objective <- match.arg(objective)

		res <- csem(model=model.description,  start=par, objective=objective,  opt.flag=0,  gradient=gradient,  opts=list("hessian"=hessian, "check.analyticals"=FALSE), ...)
		ret <- list();
		ret$f <- res$minimum
		ret$parameters <- res$estimate
		ret$C <- res$C
		ret$A <- res$A
		ret$P <- res$P
		ret$gradient <- res$gradient
		ret$hessian <- res$hessian

		return(ret)
}

msemCompiledObjective <- function(par, model.description, gradient=TRUE, hessian=FALSE, objective=c("objectiveML", "objectiveGLS", "objectiveFIML"), ...)
{
		if(missing(objective)) objective <- "objectiveML"
		objective <- match.arg(objective)

		res <- cmsem(model=model.description,  start=par, objective=objective,  opt.flag=0, gradient=gradient,   opts=list("hessian"=hessian, "check.analyticals"=FALSE), ...)
		AA <- PP <- CC <- vector(model.description$G,  mode="list")
		indAP <- 1
		indC <- 1
		ff <- numeric(model.description$G)
		for(g in 1:model.description$G)
		{
				m <- model.description$m[g]
				n <- model.description$n[g]
				AA[[g]] <- matrix(res$A[indAP:(indAP+m*m-1)], m, m);
				PP[[g]] <- matrix(res$P[indAP:(indAP+m*m-1)], m, m);
				indAP <- indAP + m*m;
				CC[[g]] <- matrix(res$C[indC:(indC+n*n-1)], n, n);
				indC <- indC + n*n;
				ff[g] <- as.numeric(res$f[g])
		}
		ret <- list();
		ret$f <- res$minimum
		ret$parameters <- res$estimate
		ret$C <- CC
		ret$A <- AA
		ret$P <- PP
		ret$ff <- ff
		ret$gradient <- res$gradient
		ret$hessian <- res$hessian

		return(ret)
}

# The wrapper function for solving optimization problems. Please note that the objective function is written in C/C++,  we need to know the name.
CompiledSolve <- function(model.description, start, objective=c("objectiveML", "objectiveGLS", "objectiveFIML", "objectivelogLik"),  gradient=TRUE, typsize=rep(1.0, length(start)), debug=FALSE, maxiter=100,...)
{
		if(missing(objective)) objective <- "objectiveML"
		objective <- match.arg(objective)

		stepmax=max(1000.0 * sqrt(sum((start/typsize)^2)),  1000.0)

		res <- csem(model=model.description, start, opt.flag=1, typsize=typsize,objective=objective,  
								gradient=gradient, 
								opts=list("iterlim"=maxiter, "print.level"=if(debug) 2 else 0,
													"hessian"=TRUE, "check.analyticals"=FALSE, "stepmax"=stepmax), ...)

		return(res)
}

# The wrapper function for solving optimization problems. Please note that the objective function is written in C/C++,  we need to know the name.
msemCompiledSolve <- function(model.description, start, objective=c("objectiveML", "objectiveGLS", "objectiveFIML"),  
															gradient=TRUE, 
															typsize=rep(1.0, length(start)), debug=FALSE, maxiter=100,gradtol=1e-6, ...)
{
		if(missing(objective)) objective <- "objectiveML"
		objective <- match.arg(objective)

		stepmax=max(1000.0 * sqrt(sum((start/typsize)^2)),  1000.0)

		res <- cmsem(model=model.description, start, opt.flag=1, typsize=typsize,objective=objective,  
								 gradient=gradient, 
								opts=list("iterlim"=maxiter, "print.level"=if(debug) 2 else 0,"gradtol"=gradtol, 
													"hessian"=TRUE, "check.analyticals"=FALSE, "stepmax"=stepmax), ...)

#reoraginize the matrix A,  P,  C 
		AA <- PP <- CC <- vector(model.description$G,  mode="list")
		indAP <- 1
		indC <- 1
		ff <- numeric(model.description$G)
		for(g in 1:model.description$G)
		{
				m <- as.integer(model.description$m[g])
				n <- as.integer(model.description$n[g])
				AA[[g]] <- matrix(as.numeric(res$A)[indAP:(indAP+m*m-1)], m, m);
				PP[[g]] <- matrix(as.numeric(res$P)[indAP:(indAP+m*m-1)], m, m);
				indAP <- indAP + m*m;
				CC[[g]] <- matrix(as.numeric(res$C)[indC:(indC+n*n-1)], n, n);
				indC <- indC + n*n;
				ff[g] <- as.numeric(res$f[g])
		}

		ret <- list();
		ret$minimum <- res$minimum
		ret$estimate <- res$estimate
		ret$gradient <- res$gradient
		ret$hessian <- res$hessian
		ret$code <- res$code
		ret$iterations <- res$iterations
		ret$C <- CC
		ret$A <- AA
		ret$P <- PP
		ret$ff <- ff

		return(ret)
}

print.f <- function(input)
{
		print(input);   # call R function "print" 
}

#optimze:0 we only compute the objective function,  gradients or hessian and return them.
# 
csem <- function(model=NULL, start=NULL,opt.flag=1,  typsize=rep(1, model$t), objective=c("objectiveML", "objectiveGLS", "objectiveFIML", "objectivelogLik", "test_objective"),  
								 gradient=TRUE, 
								 opts=list("hessian"=1, "fscale"=1, "gradtol"=1e-6, "steptol"=1e-6, "stepmax"=max(1000 * sqrt(sum((start/typsize)^2)),  1000), "iterlim"=100, 
													 "ndigit"=12,"print.level"=0, "check.analyticals"=1), 
								 csem.environment = new.env(), ...){

		environment(print.f) <- csem.environment; 
 ## Write wrappers around user-defined functions to pass additional
  ## arguments
  print.f.wrapper <- function(x){ print.f(x,...) }

		if(missing(model)) stop("Must provide the model.")
		if(missing(objective)) objective <- "objectiveML"
		objective <- match.arg(objective)
		if(missing(typsize) || is.null(typsize)) typsize <- rep(1, model$t)
		if(missing(start)) start <- rep(0.10, model$t)
		if(length(opts$print.level)==0) 
				print.level <- 0
		else 
				print.level <- as.integer(opts$print.level)
		if(print.level < 0 || print.level > 2) stop("'print.level' must be in {0, 1, 2}")


		## the following is for generating gradient.
		if(objective != "objectivelogLik")
		{
				arrows.1.seq <- model$ram[model$ram[, 1]==1 & model$ram[, 4]!=0,  4] 
				arrows.2.seq <- model$ram[model$ram[, 1]==2 & model$ram[, 4]!=0,  4]
		}

		# this function is modfied from ipoptr developed by Jelmer Ypma (http://www.ucl.ac.uk/~uctpjyy/ipoptr.html).
		# Please reference the license of ipoptr.
		get.option.types <- function(opts) {
				# define types of nlm options,  we should add all options here.
				nlm.option.types <- list(
																 "fscale"="numeric", 
																 "gradtol"="numeric", 
																 "steptol"="numeric", 
																 "stepmax"="numeric", 
																 "hessian"="integer", 
																 "iterlim"="integer", 
																 "ndigit"="integer", 
																 "print.level"="integer", 
																 "check.analyticals"="integer"
																 )


				# initialize list with options sorted by type
				converted.opts <- list( "integer"=list(), "string"=list(), "numeric"=list() )

				is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

				# check if we have at least 1 element in the list, otherwise the 
				# loop runs from 1 to down 0 and we get errors
				if ( length( opts ) > 0 ) {

						# loop over all options and give them the correct type
						for ( i in 1:length( opts ) ) {
								tmp.type <- nlm.option.types[[match( names(opts)[i], names(nlm.option.types) )]]
								if ( is.null( tmp.type ) ) {
										# determine type
										if ( is.character(opts[[i]]) ) {
												tmp.type <- "string"
										} else if ( is.wholenumber(opts[[i]]) ) {
												tmp.type <- "integer"
										} else {
												tmp.type <- "numeric"
										}
										cat( paste( "Warning: ", names(opts)[i], " is not a recognized option, we try to pass it to nlm as ", tmp.type, "\n" ) )
								}

								if ( tmp.type=="string" ) {
										converted.opts$string[[ names(opts)[i] ]] <- as.character(opts[[i]])
								} else if ( tmp.type=="integer" ) {
										converted.opts$integer[[ names(opts)[i] ]] <- as.integer(opts[[i]])
								} else if ( tmp.type=="numeric" ) {
										converted.opts$numeric[[ names(opts)[i] ]] <- as.numeric(opts[[i]])
								} else {
										stop(paste("Type of option ", names(opts)[i], " not recognized"))
								}
						}
				}

				return ( converted.opts )
		}

		if(objective != "objectivelogLik")
		{
				ret <- list( 
										"objective" = objective, 
										"gradient" = as.integer(gradient), 
										"opt.flg" = as.integer(opt.flag), 
										"start" = start, 
										"options" = get.option.types(opts), 
										"data" = model$data, 
										"pattern.number" = model$pattern.number, 
										"valid.data.patterns" = model$valid.data.patterns, 
										"S" = model$S, 
										"logdetS" = as.numeric(model$logdetS), 
										"invS" = model$invS, 
										"N" = as.integer(model$N), 
										"m" = as.integer(model$m), 
										"n" = as.integer(model$n), 
										"t" = as.integer(model$t), 
										"fixed" = model$fixed, 
										"ram" = model$ram, 
										"sel.free" = model$sel.free, 
										"arrows.1" = model$arrows.1, 
										"arrows.1.free" = model$arrows.1.free, 
										"one.head" = model$one.head, 
										"arrows.2t" = model$arrows.2t, 
										"arrows.2" = model$arrows.2, 
										"arrows.2.free" = model$arrows.2.free, 
										"unique.free.1" = model$unique.free.1, 
										"unique.free.2" = model$unique.free.2, 
										"J" = model$J, 
										"correct" = model$correct, 
										"param.names" = model$param.names, 
										"var.names" = model$var.names, 
										"one.free" = model$one.free, 
										"two.free" = model$two.free, 
										"raw" = as.integer(model$raw), 
										"arrows.1.seq" = arrows.1.seq, 
										"arrows.2.seq" = arrows.2.seq, 
										"typsize" = typsize, 
										"print.f" = print.f.wrapper, 
										"csem.environment"=csem.environment)
				attr(ret, "class") <- "csem"
		}
		else
		{
				ret <- list(
										"objective" = objective, 
										"gradient" = as.integer(gradient), 
										"opt.flg" = as.integer(opt.flag), 
										"start" = start, 
										"t" = length(start), 
										"options" = get.option.types(opts), 
										"data" = model$data, 
										"pattern.number" = model$pattern.number, 
										"valid.data.patterns" = model$valid.data.patterns, 
										"tri" = model$tri, 
										"posn.intercept" = model$posn.intercept, 
										"typsize" = typsize, 
										"print.f" = print.f.wrapper, 
										"csem.environment"=csem.environment
										)
		}
		# add the current call to the list
		# ret$call <- match.call()

		solution <- .Call("csemSolve", ret)

		# ret$environment <- NULL
		# ret$solution <- solution

		ret <- solution   #this is for simplifing the interface.
		# add solution variables to object
		#ret$status <- solution$status

		return(ret)
}

cmsem <- function(model=NULL, start=NULL,opt.flag=1,  typsize=rep(1, model$t), objective=c("objectiveML", "objectiveGLS", "objectiveFIML", "test_objective"),  
									gradient=TRUE, 
									opts=list("hessian"=1, "fscale"=1, "gradtol"=1e-6, "steptol"=1e-6, "stepmax"=max(1000 * sqrt(sum((start/typsize)^2)),  1000), "iterlim"=100, 
														"ndigit"=12,"print.level"=0, "check.analyticals"=1), 
									csem.environment = new.env(), ...){

		environment(print.f) <- csem.environment; 
 ## Write wrappers around user-defined functions to pass additional
  ## arguments
  print.f.wrapper <- function(x){ print.f(x,...) }

		if(missing(model)) stop("Must provide the model.")
		if(missing(objective)) objective <- "objectiveML"
		objective <- match.arg(objective)
		if(missing(typsize) || is.null(typsize)) typsize <- rep(1, model$t)
		if(missing(start)) start <- rep(0.10, model$t)
		if(length(opts$print.level)==0) 
				print.level <- 0
		else 
				print.level <- as.integer(opts$print.level)
		if(print.level < 0 || print.level > 2) stop("'print.level' must be in {0, 1, 2}")


		## the following is for generating gradient.
		G <- model$G
		arrows.1.seq <- arrows.2.seq <- vector(G, mode="list")
		for(g in 1:G)
		{
				arrows.1.seq[[g]] <- model$ram[[g]][model$ram[[g]][, 1]==1 & model$ram[[g]][, 4]!=0,  4] 
				arrows.2.seq[[g]] <- model$ram[[g]][model$ram[[g]][, 1]==2 & model$ram[[g]][, 4]!=0,  4]
		}

		get.option.types <- function(opts) {
				# define types of nlm options,  we should add all options here.
				nlm.option.types <- list(
																 "fscale"="numeric", 
																 "gradtol"="numeric", 
																 "steptol"="numeric", 
																 "stepmax"="numeric", 
																 "hessian"="integer", 
																 "iterlim"="integer", 
																 "ndigit"="integer", 
																 "print.level"="integer", 
																 "check.analyticals"="integer"
																 )


				# initialize list with options sorted by type
				converted.opts <- list( "integer"=list(), "string"=list(), "numeric"=list() )

				is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

				# check if we have at least 1 element in the list, otherwise the 
				# loop runs from 1 to down 0 and we get errors
				if ( length( opts ) > 0 ) {

						# loop over all options and give them the correct type
						for ( i in 1:length( opts ) ) {
								tmp.type <- nlm.option.types[[match( names(opts)[i], names(nlm.option.types) )]]
								if ( is.null( tmp.type ) ) {
										# determine type
										if ( is.character(opts[[i]]) ) {
												tmp.type <- "string"
										} else if ( is.wholenumber(opts[[i]]) ) {
												tmp.type <- "integer"
										} else {
												tmp.type <- "numeric"
										}
										cat( paste( "Warning: ", names(opts)[i], " is not a recognized option, we try to pass it to nlm as ", tmp.type, "\n" ) )
								}

								if ( tmp.type=="string" ) {
										converted.opts$string[[ names(opts)[i] ]] <- as.character(opts[[i]])
								} else if ( tmp.type=="integer" ) {
										converted.opts$integer[[ names(opts)[i] ]] <- as.integer(opts[[i]])
								} else if ( tmp.type=="numeric" ) {
										converted.opts$numeric[[ names(opts)[i] ]] <- as.numeric(opts[[i]])
								} else {
										stop(paste("Type of option ", names(opts)[i], " not recognized"))
								}
						}
				}

				return ( converted.opts )
		}

		ret <- list( 
								"objective" = objective, 
								"gradient" = as.integer(gradient), 
								"opt.flg" = as.integer(opt.flag), 
								"start" = start, 
								"options" = get.option.types(opts), 
								"G" = as.integer(model$G), 
								"data" = model$data, 
								"pattern.number" = model$pattern.number, 
								"valid.data.patterns" = model$valid.data.patterns, 
								"S" = model$S, 
								"logdetS" = model$logdetS, 
								"invS" = model$invS, 
								"N" = model$N, 
								"m" = model$m, 
								"n" = model$n, 
								"t" = as.integer(model$t), 
								"fixed" = model$fixed, 
								"ram" = model$ram, 
								"sel.free" = model$sel.free, 
								"arrows.1" = model$arrows.1, 
								"arrows.1.free" = model$arrows.1.free, 
								"one.head" = model$one.head, 
								"arrows.2t" = model$arrows.2t, 
								"arrows.2" = model$arrows.2, 
								"arrows.2.free" = model$arrows.2.free, 
								"unique.free.1" = model$unique.free.1, 
								"unique.free.2" = model$unique.free.2, 
								"J" = model$J, 
								"correct" = model$correct, 
								"param.names" = model$param.names, 
								"var.names" = model$var.names, 
								"one.free" = model$one.free, 
								"two.free" = model$two.free, 
								"raw" = as.integer(model$raw), 
								"arrows.1.seq" = arrows.1.seq, 
								"arrows.2.seq" = arrows.2.seq, 
								"typsize" = typsize, 
								"print.f" = print.f.wrapper, 
								"csem.environment"=csem.environment)
		attr(ret, "class") <- "cmsem"

		# add the current call to the list
		# ret$call <- match.call()

		solution <- .Call("cmsemSolve", ret)

		# ret$environment <- NULL
		# ret$solution <- solution

		ret <- solution   #this is for simplifing the interface.
		# add solution variables to object
		#ret$status <- solution$status

		return(ret)
}

