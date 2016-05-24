#############################################################################
#   Copyright (c) 2012 Christophe Dutang                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################
### GNE computation in GNE
###
###         R functions
### 

GNE <- function(approach = 
	c("non smooth", "fixed point", "minimization", "constrained equation"), 
	method = "default", xinit, control=list(), ...)
{
	approach <- match.arg(approach, c("non smooth", "fixed point", "minimization", "constrained equation"))
	
	if(approach == "non smooth")
		res <- GNE.nseq(xinit, method=method, control=control, ...)

	if(approach == "constr equation")
		res <- GNE.ceq(xinit, method=method, control=control, ...)
	
	if(approach == "fixed point")
		res <- GNE.fpeq(xinit, method=method, control.outer=control, ...)
	
	if(approach == "minimization")
		res <- GNE.minpb(xinit, method=method, control.outer=control, ...)
	
	c(res, list(approach=approach))
}

GNE.nseq <- function(init, dimx, dimlam, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	compl, gcompla, gcomplb, argcompl, 
	dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint, 
	method="default", control=list(), silent=TRUE, ...)
{
	if(method == "default") method <- "Newton"
	
	argtest1 <- testargfunSSR(init, dimx, dimlam, grobj, arggrobj, constr, argconstr,  grconstr, arggrconstr, 
						 compl, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint)
	
	#basic tests for funSSR
	test.try <- try( funSSR(init, dimx, dimlam, grobj, arggrobj, constr, argconstr,  
							grconstr, arggrconstr, compl, argcompl, dimmu, joint, argjoint,
							grjoint, arggrjoint), silent=silent )
	
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate Phi(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Phi(init) has infinite or NaN values.", fvec=NA) )

	argtest2 <- testargjacSSR(init, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
						 heconstr, argheconstr, gcompla, gcomplb, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint,
						 hejoint, arghejoint)	
	
	#basic tests for jacSSR
	test.try <- try( jacSSR(init, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
					  heconstr, argheconstr, gcompla, gcomplb, argcompl, dimmu, joint, argjoint,
					  grjoint, arggrjoint, hejoint, arghejoint), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evaluate Jac Phi(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Jac Phi(init) has infinite or NaN values.", fvec=NA) )
	
	#wrapped functions
	myfunSSR <- function(x, argfun, argjac)
		evalwitharglist(funSSR, x, argfun) 
	
	myjacSSR <- function(x, argfun, argjac)
		evalwitharglist(jacSSR, x, argjac)
	
	arg1 <- list(argtest1$dimx, argtest1$dimlam, argtest1$grobj, argtest1$arggrobj, 
				 argtest1$constr, argtest1$argconstr, argtest1$grconstr, argtest1$arggrconstr, 
				 argtest1$compl, argtest1$argcompl, argtest1$dimmu, 
				 argtest1$joint, argtest1$argjoint, argtest1$grjoint, argtest1$arggrjoint)
	arg2 <- list(argtest2$dimx, argtest2$dimlam, argtest2$heobj, argtest2$argheobj, 
				 argtest2$constr, argtest2$argconstr, argtest2$grconstr, argtest2$arggrconstr, 
				 argtest2$heconstr, argtest2$argheconstr, argtest2$gcompla, argtest2$gcomplb, argtest2$argcompl, 
				 argtest2$dimmu, argtest2$joint, argtest2$argjoint, argtest2$grjoint, 
				 argtest2$arggrjoint, argtest2$hejoint, argtest2$arghejoint)	
	if(!silent)
		print("init completed.")
	
	res <- nseq(init, myfunSSR, myjacSSR, argfun=arg1, argjac=arg2, method=method, 
				control=control, silent=silent, ...)	
	class(res) <- "GNE"
	if(!silent)
		print("computation completed.")
	
	res
}

GNE.ceq <- function(init, dimx, dimlam, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint, 
	method="PR", control=list(), silent=TRUE, ...)
{
	if(method == "default") method <- "PR"
	
	argtest1 <- testargfunCER(init, dimx, dimlam, grobj, arggrobj, constr, argconstr, 
							  grconstr, arggrconstr, dimmu, joint, argjoint, grjoint, arggrjoint)
	
	#basic tests for funCER
	test.try <- try( funCER(init, dimx, dimlam, grobj, arggrobj, constr, argconstr, 
				grconstr, arggrconstr, dimmu, joint, argjoint, grjoint, arggrjoint), silent=silent )

	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate H(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="H(init) has infinite or NaN values.", fvec=NA) )
	
	argtest2 <- testargjacCER(init, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, 
							  arggrconstr, heconstr, argheconstr, dimmu, joint, argjoint, grjoint, 
							  arggrjoint, hejoint, arghejoint)
	
	#basic tests for jacCER
	test.try <- try( jacCER(init, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, 
							arggrconstr, heconstr, argheconstr, dimmu, joint, argjoint,
							grjoint, arggrjoint, hejoint, arghejoint), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evaluate Jac H(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Jac H(init) has infinite or NaN values.", fvec=NA) )
	
	#wrapped functions
	myfunCER <- function(x, argfun, argjac)
		evalwitharglist(funCER, x, argfun) 
	
	myjacCER <- function(x, argfun, argjac)
		evalwitharglist(jacCER, x, argjac)
	
	arg1 <- list(dimx = argtest1$dimx, dimlam = argtest1$dimlam, grobj = argtest1$grobj, 
		arggrobj = argtest1$arggrobj, constr = argtest1$constr, argconstr = argtest1$argconstr, 
		grconstr = argtest1$grconstr, arggrconstr = argtest1$arggrconstr, dimmu = argtest1$dimmu, 
		joint = argtest1$joint, argjoint = argtest1$argjoint, grjoint = argtest1$grjoint, 
		arggrjoint = argtest1$arggrjoint)
	arg2 <- list(argtest2$dimx, argtest2$dimlam, argtest2$heobj, argtest2$argheobj, 
				 argtest2$constr, argtest2$argconstr, argtest2$grconstr, argtest2$arggrconstr, 
				 argtest2$heconstr, argtest2$argheconstr, argtest2$dimmu, argtest2$joint, 
				 argtest2$argjoint, argtest2$grjoint, argtest2$arggrjoint,
				 argtest2$hejoint, argtest2$arghejoint)	
	if(!silent)
		print("init completed.")
		
	res <- ceq(init, dimx, dimlam, Hfinal=myfunCER, jacHfinal=myjacCER, 
			   argfun=arg1, argjac=arg2, method=method, control=control, 
			   silent=silent, ...)	
	class(res) <- "GNE"
	if(!silent)
		print("computation completed.")
	
	res
}


GNE.fpeq <- function(init, dimx, obj, argobj, grobj, arggrobj, 
	heobj, argheobj, joint, argjoint, jacjoint, argjacjoint, 
	method = "default", problem = c("NIR", "VIR"), 
	merit = c("NI", "VI", "FP"), order.method=1, control.outer=list(), 
	control.inner=list(), silent=TRUE, param=list(), stepfunc, argstep, ...)
{

	if(method == "default") 
		method <- "MPE"
	method <- match.arg(method, c("pure", "vH", "UR", "RRE", "MPE", "SqRRE", "SqMPE"))
	problem <- match.arg(problem)
	merit <- match.arg(merit)
	
	if(method %in% c("vH") && merit == "FP")
		stop("incompatible method and merit arguments.")
	
	if(order.method > 3 || order.method < 1)
		stop("wrong order argument.")
	if(problem == "NIR")
	{
		arggap <- testarggapNIR(init, dimx, obj, argobj)
		argfp <- testargfpNIR(init, dimx, obj, argobj, joint, argjoint,  
					 grobj, arggrobj, jacjoint, argjacjoint)
		
		yfun <- function(x)
			fpNIR(x, argfp$dimx, argfp$obj, argfp$argobj, argfp$joint, 
				  argfp$argjoint, argfp$grobj, argfp$arggrobj, argfp$jacjoint,
				  argfp$argjacjoint, echo=!silent, control=control.inner, 
				  param=param, ...)

		if(merit == "FP")
			merit <- NULL
		else if(merit == "NI")
		{
			merit <- function(x, y=NULL)
			{
				if(is.null(y))
				{
					y <- yfun(x)
					res <- gapNIR(x, y$par, arggap$dimx, arggap$obj, arggap$argobj, 
								  param=param, echo=!silent)
					y$iter <- ifelse(!is.na(y$counts[2]), y$counts[2], y$counts[1])
					list(value=res, counts=y$counts+c(1,0), iter=y$iter)
				}else
					list(value=gapNIR(x, y, arggap$dimx, arggap$obj, arggap$argobj, 
								  param=param, echo=!silent),
						 counts=c(1,0), iter=0)
			}
		
		}else
			stop("wrong merit function")
	}
	if(problem == "VIR")
	{
		arggap <- testarggapVIR(init, dimx, grobj, arggrobj)
		argfp <- testargfpVIR(init, dimx, obj, argobj, joint, argjoint,  
							  grobj, arggrobj, jacjoint, argjacjoint)
		
		yfun <- function(x)
			fpVIR(x, argfp$dimx, argfp$obj, argfp$argobj, argfp$joint, 
				  argfp$argjoint, argfp$grobj, argfp$arggrobj, argfp$jacjoint,
				  argfp$argjacjoint, echo=!silent, control=control.inner, 
				  param=param, ...)
		if(merit == "FP")
			merit <- NULL
		else if(merit == "VI")
		{
			merit <- function(x, y=NULL)
			{
				if(is.null(y))
				{
					y <- yfun(x)
					res <- gapVIR(x, y$par, arggap$dimx, arggap$grobj, arggap$arggrobj, 
								  param=param, echo=!silent)
					y$iter <- ifelse(!is.na(y$counts[2]), y$counts[2], y$counts[1])
					list(value=res, counts=y$counts+c(1,0), iter=y$iter)
				}else
					list(value=gapVIR(x, y, arggap$dimx, arggap$grobj, arggap$arggrobj, 
								  param=param, echo=!silent),
						 counts=c(1,0), iter=0)
			}
			
		}else
			stop("wrong merit function")
		
	}	
		
	if(!silent)
	{
		print("init completed.")
		cat("control parameters for fpeq.\n")
		print(control.inner)
		print(control.outer)
	}
	res <- fpeq(init, fn=yfun, merit=merit, method=method, control=control.outer, stepfunc=stepfunc,
				argstep=argstep, silent=silent, order.method=order.method, ...)
	class(res) <- "GNE"
	if(!silent)
		print("computation completed.")
	res
}


GNE.minpb <- function(init, dimx, obj, argobj, grobj, arggrobj, 
	heobj, argheobj, joint, argjoint, jacjoint, argjacjoint, 
	method="default", problem = c("NIR", "VIR"), control.outer=list(), 
	control.inner=list(), silent=TRUE, param=list(), optim.type=c("free","constr"), ...)
{
	if(method == "default")
		method <- "BFGS"
	method <- match.arg(method, c("BB","BFGS", "CG"))
	problem <- match.arg(problem)
	optim.type <- match.arg(optim.type)
	if(length(param) != 0)
		if(any(!names(param) %in% c("alpha", "beta")))
			stop("argument param should be a named list with alpha and optionnaly beta")
	
	if(problem == "NIR")
	{
		arggap <- testarggapNIR(init, dimx, obj, argobj)
		arggradgap <- testarggradgapNIR(init, dimx, grobj, arggrobj)
		argfp <- testargfpNIR(init, dimx, obj, argobj, joint, argjoint,  
							  grobj, arggrobj, jacjoint, argjacjoint)
		
		yfun <- function(x, param=param)
			fpNIR(x, argfp$dimx, argfp$obj, argfp$argobj, argfp$joint, 
			  argfp$argjoint, argfp$grobj, argfp$arggrobj, argfp$jacjoint,
			  argfp$argjacjoint, echo=!silent, control=control.inner, 
			  param=param, ...)
		
		if(length(param) <= 1)
		{
			merit <- function(x, y=NULL)
			{
				if(is.null(y))
				{
					y <- yfun(x)
					res <- gapNIR(x, y$par, arggap$dimx, arggap$obj, arggap$argobj, 
								  param=param, echo=!silent)
					y$iter <- ifelse(!is.na(y$counts[2]), y$counts[2], y$counts[1])
					list(value=res, counts=y$counts+c(1,0), iter=y$iter)
				}else
					list(value=gapNIR(x, y, arggap$dimx, arggap$obj, arggap$argobj, 
								  param=param, echo=!silent), counts=c(1,0), iter=0)
			}
			gradmerit <- function(x, y=NULL)
			{
				if(is.null(y))
				{
					y <- yfun(x)
					res <- gradxgapNIR(x, y$par, arggradgap$dimx, arggradgap$grobj, 
									   arggradgap$arggrobj, param=param, echo=!silent)
					y$iter <- ifelse(!is.na(y$counts[2]), y$counts[2], y$counts[1])
					list(value=res, counts=y$counts+c(0,1), iter=y$iter)
				}else
					list(value=gradxgapNIR(x, y, arggradgap$dimx, arggradgap$grobj, 
						arggradgap$arggrobj, param=param, echo=!silent), counts=c(0,1), iter=0)
			}
		}else if(length(param) == 2)
		{
			merit <- function(x, y=NULL)
			{
				yalpha <- yfun(x, param["alpha"])
				Valpha <- gapNIR(x, yalpha$par, arggap$dimx, arggap$obj, arggap$argobj, 
								  param=param["alpha"], echo=!silent)
				yalpha$iter <- ifelse(!is.na(yalpha$counts[2]), yalpha$counts[2], yalpha$counts[1])
				
				ybeta <- yfun(x, param["beta"])
				Vbeta <- gapNIR(x, ybeta$par, arggap$dimx, arggap$obj, arggap$argobj, 
								 param=param["beta"], echo=!silent)
				ybeta$iter <- ifelse(!is.na(ybeta$counts[2]), ybeta$counts[2], ybeta$counts[1])
				
				list(value=(Valpha-Vbeta)^2, counts=yalpha$counts+ybeta$counts+c(1,0), iter=yalpha$iter+ybeta$iter)
			}
			gradmerit <- function(x, y=NULL)
			{
				yalpha <- yfun(x, param["alpha"])
				Valpha <- gapNIR(x, yalpha$par, arggap$dimx, arggap$obj, arggap$argobj, 
								 param=param["alpha"], echo=!silent)
				grValpha <- gradxgapNIR(x, yalpha$par, arggradgap$dimx, arggradgap$grobj, 
									   arggradgap$arggrobj, param=param["alpha"], echo=!silent)
				yalpha$iter <- ifelse(!is.na(yalpha$counts[2]), yalpha$counts[2], yalpha$counts[1])

				ybeta <- yfun(x, param["beta"])
				Vbeta <- gapNIR(x, ybeta$par, arggap$dimx, arggap$obj, arggap$argobj, 
								param=param["beta"], echo=!silent)
				grVbeta <- gradxgapNIR(x, ybeta$par, arggradgap$dimx, arggradgap$grobj, 
										arggradgap$arggrobj, param=param["beta"], echo=!silent)
				ybeta$iter <- ifelse(!is.na(ybeta$counts[2]), ybeta$counts[2], ybeta$counts[1])
				
				list(value=2*(grValpha-grVbeta)*(Valpha-Vbeta), counts=yalpha$counts+ybeta$counts+c(1,1), iter=yalpha$iter+ybeta$iter)
			}
			
		}else
			stop("wrong argument param.")
	}
	if(problem == "VIR")
	{
		arggap <- testarggapVIR(init, dimx, grobj, arggrobj)
		arggradgap <- testarggradxgapVIR(init, dimx, grobj, arggrobj, heobj, argheobj)
		argfp <- testargfpVIR(init, dimx, obj, argobj, joint, argjoint,  
							  grobj, arggrobj, jacjoint, argjacjoint)

		yfun <- function(x, param=param)
			fpVIR(x, argfp$dimx, argfp$obj, argfp$argobj, argfp$joint, 
				  argfp$argjoint, argfp$grobj, argfp$arggrobj, argfp$jacjoint,
				  argfp$argjacjoint, echo=!silent, control=control.inner, 
				  param=param, ...)
		
		if(length(param) <= 1)
		{
			
			merit <- function(x, y=NULL)
			{
				if(is.null(y))
				{
					y <- yfun(x)
					res <- gapVIR(x, y$par, arggap$dimx, arggap$grobj, arggap$arggrobj, 
								  param=param, echo=!silent)
					y$iter <- ifelse(!is.na(y$counts[2]), y$counts[2], y$counts[1])
					list(value=res, counts=y$counts+c(1,0), iter=y$iter)
				}else
					list(value=gapVIR(x, y, arggap$dimx, arggap$grobj, arggap$arggrobj, 
								  param=param, echo=!silent), counts=c(1,0), iter=0)
			}
			gradmerit <- function(x, y=NULL)
			{
				if(is.null(y))
				{
					y <- yfun(x)
					res <- gradxgapVIR(x, y$par, arggradgap$dimx, arggradgap$grobj, 
									   arggradgap$arggrobj, arggradgap$heobj, 
									   arggradgap$argheobj, param=param, echo=!silent)
					y$iter <- ifelse(!is.na(y$counts[2]), y$counts[2], y$counts[1])
					list(value=res, counts=y$counts+c(0,1), iter=y$iter)
				}else
					list(value=gradxgapVIR(x, y, arggradgap$dimx, arggradgap$grobj, 
									   arggradgap$arggrobj, arggradgap$heobj, 
									   arggradgap$argheobj, param=param, echo=!silent), counts=c(0,1), iter=0)
			}
		}else if(length(param) == 2)
		{
			merit <- function(x, y=NULL)
			{
				yalpha <- yfun(x, param["alpha"])
				Valpha <- gapVIR(x, yalpha$par, arggap$dimx, arggap$grobj, arggap$arggrobj, 
							  param=param["alpha"], echo=!silent)
				yalpha$iter <- ifelse(!is.na(yalpha$counts[2]), yalpha$counts[2], yalpha$counts[1])
				
				ybeta <- yfun(x, param["beta"])
				Vbeta <- gapVIR(x, ybeta$par, arggap$dimx, arggap$grobj, arggap$arggrobj, 
								 param=param["beta"], echo=!silent)
				ybeta$iter <- ifelse(!is.na(ybeta$counts[2]), ybeta$counts[2], ybeta$counts[1])
				
				list(value=(Valpha-Vbeta)^2, counts=yalpha$counts+ybeta$counts+c(1,1), iter=yalpha$iter+ybeta$iter)
			}
			gradmerit <- function(x, y=NULL)
			{
				yalpha <- yfun(x, param["alpha"])
				Valpha <- gapVIR(x, yalpha$par, arggap$dimx, arggap$grobj, arggap$arggrobj, 
								 param=param["alpha"], echo=!silent)
				grValpha <- gradxgapVIR(x, yalpha$par, arggradgap$dimx, arggradgap$grobj, 
								   arggradgap$arggrobj, arggradgap$heobj, 
								   arggradgap$argheobj, param=param["alpha"], echo=!silent)
				
				yalpha$iter <- ifelse(!is.na(yalpha$counts[2]), yalpha$counts[2], yalpha$counts[1])

				ybeta <- yfun(x, param["beta"])
				Vbeta <- gapVIR(x, ybeta$par, arggap$dimx, arggap$grobj, arggap$arggrobj, 
								param=param["beta"], echo=!silent)
				grVbeta <- gradxgapVIR(x, ybeta$par, arggradgap$dimx, arggradgap$grobj, 
										arggradgap$arggrobj, arggradgap$heobj, 
										arggradgap$argheobj, param=param["beta"], echo=!silent)
				
				ybeta$iter <- ifelse(!is.na(ybeta$counts[2]), ybeta$counts[2], ybeta$counts[1])

				list(value=2*(grValpha-grVbeta)*(Valpha-Vbeta), counts=yalpha$counts+ybeta$counts+c(1,1), iter=yalpha$iter+ybeta$iter)
			}
			
		}else
			stop("wrong argument param.")	
	}
	

	
	if(optim.type == "constr")
		res <- minpb(init, merit, gr=gradmerit, hin=joint, arghin=argjoint, 
			hin.jac=jacjoint, arghin.jac=argjacjoint, method=method, 
			control=control.outer, silent=silent, ...)
	if(optim.type == "free")
		res <- minpb(init, merit, gr=gradmerit, method=method, 
				 control=control.outer, silent=silent, ...)
	
	class(res) <- "GNE"
	res
}





#print function
print.GNE <- function(x, ...)
{
	if (!inherits(x, "GNE"))
		stop("Use only with 'GNE' objects")
	cat("GNE:", x$par, "\nwith optimal norm", x$value, "\n")
	cat("after ", x$iter, "iterations with exit code", x$code, ".\n")
	cat("Output message:", x$message, "\n")
	if(!is.null(x$counts))	
		cat("Function/grad/hessian calls:", x$counts, "\n")
	if(!is.null(x$outer.counts))	
		cat("Outer Function/grad/hessian calls:", x$outer.counts, "\n")	
	if(!is.null(x$inner.counts))	
		cat("Inner Function/grad/hessian calls:", x$inner.counts, "\n")	
	if(!is.null(x$fvec))
		cat("Optimal (vector) value:", x$fvec, "\n")
}


#summary function
summary.GNE <- function(object, ...)
{
	structure(object, class = c("summary.GNE", class(object)))	
}


#print function
print.summary.GNE <- function(x, ...)
{	
	if (!inherits(x, "GNE"))
		stop("Use only with 'GNE' objects")
	cat("GNE:", x$par, "\nwith optimal norm", x$value, "\n")
	cat("after ", x$iter, "iterations with exit code", x$code, ".\n")
	if(!is.null(x$playerdead))
		cat("dead players? ", x$playerdead, ".\n")
}

