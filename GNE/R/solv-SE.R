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
### SE computation in GNE
###
###         R functions
### 

#Stackelberg equilibrium computation


SE.nseq <- function(leaders, init, dimx, dimlam, 
	obj, argobj, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	compl, gcompla, gcomplb, argcompl, 
	dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint, 
	method.follower="default", method.leader="default", 
	control.follower=list(), control.leader=list(), 
	maxit.follower=10, silent=TRUE, 
	simpleconstr=FALSE, ...)
{
	if(method.follower == "default") method.follower <- "Newton"
	if(method.leader == "default") method.leader <- "BFGS"	
	
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
	
	nbplay <- argtest1$nplayer

	argtest3 <- testarggapNIR(init[1:nbplay], dimx, obj, argobj)

	if(!is.numeric(leaders) || length(leaders) > nbplay-1)
		stop("wrong leaders argument.")
	if(any(!leaders %in% 1:nbplay))
		stop("wrong leaders argument.")
	followers <- (1:nbplay)[!(1:nbplay %in% leaders)]
	
	dimx <- argtest1$dimx
	n <- sum(dimx)
	nfoll <- sum(dimx[followers])
	nlead <- sum(dimx[leaders])
	dimlam <- argtest1$dimlam
	m <- sum(dimlam)
	
	#1st row is the begin index, 2nd row the end index
	index4lam <- rbind( cumsum(dimlam) - dimlam + 1, cumsum(dimlam) ) + n
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) ) + 0	
		
	index4xfoll <- as.vector(sapply(followers, function(i) index4x[1,i]:index4x[2,i]))
	index4xlead <- as.vector(sapply(leaders, function(i) index4x[1,i]:index4x[2,i]))
	index4lamfoll <- as.vector(sapply(followers, function(i) index4lam[1,i]:index4lam[2,i]))
	index4mufoll <- (1:length(init))[-(1:(n+m))]
	
	
	bestrespcount <- list(phicnt=0, jaccnt=0)
	objleadcount <- list(fn=0)

	#wrap objective gradient (x, i, j, arg)
	#1st arg=x, 2nd arg=id player, 3rd arg= deriv index, 4th arg=add. arg
	grobjfoll <- function(xfoll, play, d1, arg)
		transfoll(argtest1$grobj, xfoll, index4xfoll, index4xlead, followers, leaders, 
				  arg, arg$foll[play], arg$foll[d1], arg$add) 
				
	#wrap objective hessian (x, i, j, k, arg)
	#1st arg=x, 2nd arg=id player, 3rd arg= deriv index, 4th arg= deriv index, 5th arg=add. arg
	heobjfoll <- function(xfoll, play, d1, d2, arg)
		transfoll(argtest2$heobj, xfoll, index4xfoll, index4xlead, followers, leaders, 
				  arg, arg$foll[play], arg$foll[d1], arg$foll[d2], arg$add)
		
	#wrap constraint function (x, i, arg)
	#1st arg=x, 2nd arg=id player, 3rd arg= add. arg
	if(!is.null(argtest1$constr))
		constrfoll <- function(xfoll, play, arg)
			transfoll(argtest1$constr, xfoll, index4xfoll, index4xlead, followers, leaders, 
					  arg, arg$foll[play], arg$add)
	else
		constrfoll <- NULL
	#wrap constraint gradient (x, i, j, arg)
	#1st arg=x, 2nd arg=id player, 3rd arg= deriv index, 4th arg=add. arg
	if(!is.null(argtest1$grconstr))
		grconstrfoll <- function(xfoll, play, d1, arg)
			transfoll(argtest1$grconstr, xfoll, index4xfoll, index4xlead, followers, leaders, 
					  arg, arg$foll[play], arg$foll[d1], arg$add)
	else
		grconstrfoll <- NULL
	#wrap constraint hessian (x, i, j, k, arg)
	#1st arg=x, 2nd arg=id player, 3rd arg= deriv index, 4th arg= deriv index, 5th arg=add. arg
	if(!is.null(argtest2$heconstr))
		heconstrfoll <- function(xfoll, play, d1, d2, arg)
			transfoll(argtest2$heconstr, xfoll, index4xfoll, index4xlead, followers, leaders, 
					  arg, arg$foll[play], arg$foll[d1], arg$foll[d2], arg$add)
	else
		heconstrfoll <- NULL
	
	
	#wrap joint function (x, arg)
	#1st arg=x, 2nd arg=add. arg
	if(!is.null(argtest1$joint))
		jointfoll <- function(xfoll, arg)
			transfoll(argtest1$joint, xfoll, index4xfoll, index4xlead, followers, leaders, 
					  arg, arg$add)
	else
		jointfoll <- NULL
	#wrap joint gradient (x, j, arg)
	#1st arg=x, 2rd arg= deriv index, 3th arg=add. arg
	if(!is.null(argtest1$grjoint))
		grjointfoll <- function(xfoll, d1, arg)
			transfoll(argtest1$grjoint, xfoll, index4xfoll, index4xlead, followers, leaders, 
					  arg, arg$foll[d1], arg$add)
	else
		grjointfoll <- NULL
	#wrap joint hessian (x, j, k, arg)
	#1st arg=x, 3rd arg= deriv index, 4th arg= deriv index, 5th arg=add. arg
	if(!is.null(argtest2$hejoint))
		hejointfoll <- function(xfoll, d1, d2, arg)
			transfoll(argtest2$hejoint, xfoll, index4xfoll, index4xlead, followers, leaders, 
					  arg, arg$foll[d1], arg$foll[d2], arg$add)
	else
		hejointfoll <- NULL
	
	
	listfollfunc <- list(grobjfoll=grobjfoll, heobjfoll=heobjfoll,
		constrfoll=constrfoll, grconstrfoll=grconstrfoll, heconstrfoll=heconstrfoll,
		jointfoll=jointfoll, grjointfoll=grjointfoll, hejointfoll=hejointfoll)
	
#compute the objective of leaders for a strategy x of leaders
#and corresponding followers actions
	objleaders <- function(xlead, arg1, arg2, arg3, leaders, followers,
						   id4xfoll, id4lamfoll, id4mufoll, init, id4xlead,
						   follfun, method.follower, control.follower)
	{
#		cat("1-", id4xfoll, "index4xfoll\t", id4lamfoll, "index4lamfoll\t", id4mufoll, "index4mufoll", "\n")
#		
#		print(xlead)
#		print(attributes(arg1))
#		print(attributes(arg2))
#		print(attributes(arg3))	
#		print(leaders)
		n <- sum(arg1$dimx)
		nfoll <- sum(arg1$dimx[followers])
		x <- rep(NA, n)
		init2 <- c(rep(xlead, nbplay), rep(1e-3, length(init)-nbplay))
		
		foll <- bestresponse(xlead, arg1, arg2, leaders, followers,
							 id4xfoll, id4lamfoll, id4mufoll, init2, 
							 follfun, method.follower=method.follower, 
							 control.follower=control.follower, 
							 maxit.follower=maxit.follower, ...)
		bestrespcount$phicnt <<- bestrespcount$phicnt + foll$counts["phicnt"]
		bestrespcount$jaccnt <<- bestrespcount$jaccnt + foll$counts["jaccnt"]
		
		x[id4xlead] <- xlead
		x[id4xfoll] <- foll$par[1:nfoll]
		
		objleaders <- sapply(leaders, function(i) arg3$obj(x, i, arg3$argobj))
		
		objleadcount$fn <<- objleadcount$fn + 1
		

		if(length(objleaders) > 1)
			return( sqrt(sum(objleaders^2)) )
		else
			return( objleaders )
	}
	
	

	argfnlist <- list(arg1=argtest1, arg2=argtest2, arg3=argtest3, 
					  leaders=leaders, followers=followers,
					  id4xfoll=index4xfoll, id4lamfoll=index4lamfoll, 
					  id4mufoll=index4mufoll, init=init, follfun=listfollfunc,
					  method.follower=method.follower, control.follower=control.follower)

#	cat("\nblii\n")	
#	print(evalwitharglist(objleaders, init[index4xlead], argfnlist))
#	cat("blii\n\n")
#	stop("here")
	
#	bestresponse(xlead, arg1, arg2, leaders, followers, 
#				 id4xfoll, id4lamfoll, id4mufoll, init, follfun, 
#				 method.follower, control.follower, ...)
	
	#create the constraint function of leaders given best response of followers
#such that constr(x) <= 0 as in input we have constr(x) <= 0
	if(is.null(argtest1$joint) && is.null(argtest1$constr))
	{
		constrleaders <- NULL
	}else if(is.null(argtest1$joint) && !is.null(argtest1$constr)) 
	{
		if(!simpleconstr)
		{
			constrleaders <- function(xlead, arg1, arg2, leaders, followers,
								  id4xfoll, id4lamfoll, id4mufoll, init, 
								  follfun, method.follower, control.follower,
								  maxit.follower)
			{
			
			nfoll <- sum(arg1$dimx[followers])
			foll <- bestresponse(xlead, arg1, arg2, leaders, followers,
								 id4xfoll, id4lamfoll, id4mufoll, init, 
								 follfun, method.follower=method.follower, 
								 control.follower=control.follower,
								 maxit.follower=maxit.follower)
			bestrespcount$phicnt <<- bestrespcount$phicnt + foll$counts["phicnt"]
			bestrespcount$jaccnt <<- bestrespcount$jaccnt + foll$counts["jaccnt"]
		
			
			x <- rep(NA, n)
			x[index4xlead] <- xlead
			x[id4xfoll] <- foll$par[1:nfoll]	
			sapply(leaders, function(i) arg1$constr(x, i, arg1$argconstr))
			}	
			argconstrlist <- list(arg1=argtest1, arg2=argtest2, 
							  leaders=leaders, followers=followers,
							  id4xfoll=index4xfoll, id4lamfoll=index4lamfoll, 
							  id4mufoll=index4mufoll, init=init, follfun=listfollfunc,
							  method.follower=method.follower, 
							  control.follower=control.follower,
							  maxit.follower=maxit.follower)
		}else
		{
			constrleaders <- function(xlead, arg1, arg2, leaders, followers,
									  id4xfoll, id4lamfoll, id4mufoll, init, 
									  follfun, method.follower, control.follower,
									  maxit.follower)
			{
				
				nfoll <- sum(arg1$dimx[followers])
				foll <- rep(1, nfoll)
				
				x <- rep(NA, n)
				x[index4xlead] <- xlead
				x[id4xfoll] <- rep(1, nfoll)	
				sapply(leaders, function(i) arg1$constr(x, i, arg1$argconstr))
			}	
			argconstrlist <- list(arg1=argtest1, arg2=argtest2, 
								  leaders=leaders, followers=followers,
								  id4xfoll=index4xfoll, id4lamfoll=index4lamfoll, 
								  id4mufoll=index4mufoll, init=NULL, follfun=NULL,
								  method.follower=NULL, 
								  control.follower=NULL,
								  maxit.follower=NULL)			
		}
		
	}else if(!is.null(argtest1$joint) && is.null(argtest1$constr)) 
	{
		constrleaders <- function(xlead, arg1, arg2, leaders, followers,
								  id4xfoll, id4lamfoll, id4mufoll, init, 
								  follfun, method.follower, control.follower,
								  maxit.follower)
		{
			nfoll <- sum(arg1$dimx[followers])
			foll <- bestresponse(xlead, arg1, arg2, leaders, followers,
								 id4xfoll, id4lamfoll, id4mufoll, init, 
								 follfun, method.follower=method.follower, 
								 control.follower=control.follower,
								 maxit.follower=maxit.follower)
			bestrespcount$phicnt <<- bestrespcount$phicnt + foll$counts["phicnt"]
			bestrespcount$jaccnt <<- bestrespcount$jaccnt + foll$counts["jaccnt"]

			x <- rep(NA, n)
			x[index4xlead] <- xlead
			x[id4xfoll] <- foll$par[1:nfoll]	
			arg1$joint(x, arg1$argjoint)
		}
		argconstrlist <- list(arg1=argtest1, arg2=argtest2, 
							  leaders=leaders, followers=followers,
							  id4xfoll=index4xfoll, id4lamfoll=index4lamfoll, 
							  id4mufoll=index4mufoll, init=init, follfun=listfollfunc,
							  method.follower=method.follower, 
							  control.follower=control.follower,
							  maxit.follower=maxit.follower)		
	}else
	{
		constrleaders <- function(xlead, arg1, arg2, leaders, followers,
								  id4xfoll, id4lamfoll, id4mufoll, init, 
								  follfun, method.follower, control.follower,
								  maxit.follower)
		{
			nfoll <- sum(arg1$dimx[followers])
			foll <- bestresponse(xlead, arg1, arg2, leaders, followers,
								 id4xfoll, id4lamfoll, id4mufoll, init, 
								 follfun, method.follower=method.follower, 
								 control.follower=control.follower,
								 maxit.follower=maxit.follower)
			bestrespcount$phicnt <<- bestrespcount$phicnt + foll$counts["phicnt"]
			bestrespcount$jaccnt <<- bestrespcount$jaccnt + foll$counts["jaccnt"]

			x <- rep(NA, n)
			x[index4xlead] <- xlead
			x[id4xfoll] <- foll$par[1:nfoll]	
			y <- arg1$joint(x, arg1$argjoint)
			z <- sapply(leaders, function(i) arg1$constr(x, i, arg1$argconstr))
			c(y, z)
		}
		argconstrlist <- list(arg1=argtest1, arg2=argtest2, 
							  leaders=leaders, followers=followers,
							  id4xfoll=index4xfoll, id4lamfoll=index4lamfoll, 
							  id4mufoll=index4mufoll, init=init, follfun=listfollfunc,
							  method.follower=method.follower, 
							  control.follower=control.follower,
							  maxit.follower=maxit.follower)	
	}
	
	if(!silent)
		cat("start computation of SE\t")
	
	
	#computation of Stackelberg equilibria
	if(is.null(argtest1$joint) && is.null(argtest1$constr))
	{
		if(!silent)
			cat("no constraint function\n")
		
		
		reslead <- minpb(init[index4xlead], fn=objleaders, method=method.leader, 
			control=control.leader, argfn=argfnlist, ...)
		if(reslead$code == 100)
			return(reslead)
		else if(reslead$code != 0)
			warning("Non-optimal Stackelberg equilibrium.")

		
		parlead <- reslead$par[1:nlead]
		resval <- evalwitharglist(objleaders, parlead, argfnlist)
		resfoll <- bestresponse(parlead, argtest1, argtest2, leaders, followers,
								id4xfoll=index4xfoll, id4lamfoll=index4lamfoll, 
								id4mufoll=index4mufoll, init, listfollfunc, 
								method.follower=method.follower, 
								control.follower=control.follower,
								maxit.follower=maxit.follower)
		bestrespcount$phicnt <- bestrespcount$phicnt + resfoll$counts["phicnt"]
		bestrespcount$jaccnt <- bestrespcount$jaccnt + resfoll$counts["jaccnt"]
		
		parfol <- resfoll$par[1:nfoll]
			
		res <- list(par = c(parlead, parfol), value = parlead, 
			counts = list(leadfn= objleadcount$fn, follfn=bestrespcount$phicnt, 
				folljac=bestrespcount$jaccnt), 
			iter = reslead$iter, code = reslead$code, message = reslead$message)		
	
	}else 
	{
		if(!silent)
			cat("with constraint functions\n")

		if(is.null(constrleaders))
			stop("internal error in SE.nseq.")
				 
		
		reslead <- minpb(init[index4xlead], fn=objleaders, method=method.leader, 
						 hin=constrleaders, arghin=argconstrlist, 
						 control=control.leader, argfn=argfnlist, silent=silent, ...)
		if(reslead$code == 100)
			return(reslead)
		else if(reslead$code != 0)
			warning("Non-optimal Stackelberg equilibrium.")
		
		parlead <- reslead$par[1:nlead]
		resval <- evalwitharglist(objleaders, parlead, argfnlist)
		resfoll <- bestresponse(parlead, argtest1, argtest2, leaders, followers,
								id4xfoll=index4xfoll, id4lamfoll=index4lamfoll, 
								id4mufoll=index4mufoll, init, listfollfunc, 
								method.follower=method.follower, 
								control.follower=control.follower,
								maxit.follower=maxit.follower)
		bestrespcount$phicnt <- bestrespcount$phicnt + resfoll$counts["phicnt"]
		bestrespcount$jaccnt <- bestrespcount$jaccnt + resfoll$counts["jaccnt"]
		
		parfol <- resfoll$par[1:nfoll]
		
		res <- list(par = c(parlead, parfol), value = parlead, 
					counts = list(leadfn= objleadcount$fn, follfn=bestrespcount$phicnt, 
								  folljac=bestrespcount$jaccnt), 
					iter = reslead$iter, code = reslead$code, message = reslead$message)		
	}
		
	if(!silent)
		cat("end computation of SE\n")
	
	res
}

	
#compute the best response of followers for a given strategy x of leaders
bestresponse <- function(xlead, arg1, arg2, leaders, followers, 
	id4xfoll, id4lamfoll, id4mufoll, init, follfun, 
	method.follower, control.follower, maxit.follower=10, ...)
{
#		cat("2-", id4xfoll, "index4xfoll\t", id4lamfoll, "index4lamfoll\t", id4mufoll, "index4mufoll", "\n")
	
#	print(sapply(follfun, class))
	if(!is.list(follfun))
		stop("wrong type for follfun.")
	if(!all(sapply(follfun, is.function) | sapply(follfun, is.null)))
		stop("wrong argument follfun.")
	
	dimx <- arg1$dimx
	dimlam <- arg1$dimlam
	
	nfoll <- sum(dimx[followers])
	nbplay <- arg1$nplayer
	
	xfoll <- init[id4xfoll]
	lamfoll <- init[id4lamfoll]
	mufoll <- init[id4mufoll]
	
	if(is.null(follfun$constrfoll) && is.null(follfun$jointfoll))
		initfoll <- xfoll
	else if(!is.null(follfun$constrfoll) && is.null(follfun$jointfoll))
		initfoll <- c(xfoll, lamfoll)
	else if(is.null(follfun$constrfoll) && !is.null(follfun$jointfoll))
		initfoll <- c(xfoll, mufoll)
	else 
		initfoll <- c(xfoll, lamfoll, mufoll)
	
	
#		cat(xfoll, lamfoll, mufoll, "\n")
#		print(initfoll)
		
	
	arggrobjSE <- list(xlead=xlead, add=arg1$arggrobj, lead=leaders, foll=followers, nbplayer=nbplay)
	argheobjSE <- list(xlead=xlead, add=arg2$argheobj, lead=leaders, foll=followers, nbplayer=nbplay)
	argconstrSE <- list(xlead=xlead, add=arg1$argconstr, lead=leaders, foll=followers, nbplayer=nbplay)
	arggrconstrSE <- list(xlead=xlead, add=arg1$arggrconstr, lead=leaders, foll=followers, nbplayer=nbplay)
	argheconstrSE <- list(xlead=xlead, add=arg2$argheconstr, lead=leaders, foll=followers, nbplayer=nbplay)
	argjointSE <- list(xlead=xlead, add=arg1$argjoint, lead=leaders, foll=followers, nbplayer=nbplay)
	arggrjointSE <- list(xlead=xlead, add=arg1$arggrjoint, lead=leaders, foll=followers, nbplayer=nbplay)
	arghejointSE <- list(xlead=xlead, add=arg2$arghejoint, lead=leaders, foll=followers, nbplayer=nbplay)
	
	arg1SE <- list(dimx = dimx[followers], dimlam = dimlam[followers], 
				   grobj = follfun$grobjfoll, arggrobj = arggrobjSE, 
				   constr = follfun$constrfoll, argconstr = argconstrSE, 
				   grconstr = follfun$grconstrfoll, arggrconstr = arggrconstrSE, 
				   compl = arg1$compl, argcompl = arg1$argcompl, 
				   dimmu = arg1$dimmu, joint = follfun$jointfoll, 
				   argjoint = argjointSE, grjoint = follfun$grjointfoll, 
				   arggrjoint = arggrjointSE)
	
	arg2SE <- list(dimx = dimx[followers], dimlam = dimlam[followers], 
				   heobj = follfun$heobjfoll, argheobj = argheobjSE, 
				   constr = follfun$constrfoll, argconstr = argconstrSE, 
				   grconstr = follfun$grconstrfoll, arggrconstr = arggrconstrSE,
				   heconstr = follfun$heconstrfoll, argheconstr = argheconstrSE, 
				   gcompla = arg2$gcompla, gcomplb = arg2$gcomplb, argcompl = arg2$argcompl, 
				   dimmu = arg2$dimmu, joint = follfun$jointfoll, argjoint = argjointSE, 
				   grjoint = follfun$grjointfoll, arggrjoint = arggrjointSE, 
				   hejoint = follfun$hejointfoll, arghejoint = arghejointSE)	
#		cat("blii\n")
#		print(sapply(arg1SE, is.null))
#		print(sapply(arg1SE, length))
#		cat("blii2\n")
#		print(sapply(arg2SE, is.null))
#		print(sapply(arg2SE, length))
#		cat("dimx", arg1SE$dimx, "dimlam", arg1SE$dimlam, "dimmu", arg1SE$dimmu, "\n")
	
	checkerror <- length(initfoll) != sum(arg1SE$dimx)+sum(arg1SE$dimlam)+sum(arg1SE$dimmu)	 
	if(checkerror)
		stop("internal error in bestresponse.")
	
	#wrapped functions
	myfunSSR <- function(x, argfun, argjac)
	evalwitharglist(funSSR, x, argfun) 
	
	myjacSSR <- function(x, argfun, argjac)
	evalwitharglist(jacSSR, x, argjac)
	
	res <- list(code=99, value=Inf)
	iter <- 0
	
	while(res$code != 1 && iter < maxit.follower)
	{
		if(iter > 0)
			initfoll <- initfoll*(1+rnorm(length(initfoll), 0, .1))
		res2 <- nseq(initfoll, myfunSSR, myjacSSR, argfun=arg1SE, argjac=arg2SE, 
				method=method.follower, control=control.follower, ...)	
		iter <- iter + 1 
		if(res2$code == 100)
			stop(res2$message)

#		cat("iter", iter, res2$value,"\n")
#		print(cbind(init=initfoll, opt=res2$par))
		if(res2$value < res$value)
			res <- res2

	}
	if(res$code != 1)
		warning("Non-optimal point when computing best response")
	
	res
}	

#transform a function to a leader/follower setting
transfoll <- function(f, xfoll, id4xfoll, id4xlead, followers, leaders, arg, ...)
{
#		cat("n", n, "nlead+nfoll", nlead+nfoll, "bestrespcount$phicnt", bestrespcount$phicnt,"\n")
#print(index4xlead)		
#		print(attributes(arg))
#		print(length(arg))
	
	n <- sum(arg$dimx)
	nfoll <- sum(arg$dimx[followers])
	nlead <- sum(arg$dimx[leaders])
	
	x <- rep(NA, n)
	x[id4xlead] <- arg$xlead[1:nlead]
	x[id4xfoll] <- xfoll[1:nfoll]
	f(x, ...)
}
	



SE.objleaders <- function(x, leaders, init, dimx, dimlam, 
	obj, argobj, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	compl, gcompla, gcomplb, argcompl, 
	dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint, 
	method.follower="default", control.follower=list(), 
	maxit.follower=10, silent=TRUE, ...)
{
	if(method.follower == "default") method.follower <- "Newton"

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
	
	nbplay <- argtest1$nplayer
	
	argtest3 <- testarggapNIR(init[1:nbplay], dimx, obj, argobj)
	
	if(!is.numeric(leaders) || length(leaders) > nbplay-1)
		stop("wrong leaders argument.")
	if(any(!leaders %in% 1:nbplay))
		stop("wrong leaders argument.")
	followers <- (1:nbplay)[!(1:nbplay %in% leaders)]
	
	dimx <- argtest1$dimx
	n <- sum(dimx)
	nfoll <- sum(dimx[followers])
	nlead <- sum(dimx[leaders])
	dimlam <- argtest1$dimlam
	m <- sum(dimlam)
	#1st row is the begin index, 2nd row the end index
	index4lam <- rbind( cumsum(dimlam) - dimlam + 1, cumsum(dimlam) ) + n
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) ) + 0	
	
	index4xfoll <- as.vector(sapply(followers, function(i) index4x[1,i]:index4x[2,i]))
	index4xlead <- as.vector(sapply(leaders, function(i) index4x[1,i]:index4x[2,i]))
	index4lamfoll <- as.vector(sapply(followers, function(i) index4lam[1,i]:index4lam[2,i]))
	index4mufoll <- (1:length(init))[-(1:(n+m))]
	
	grobjfoll <- function(xfoll, play, d1, arg)
		transfoll(argtest1$grobj, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$foll[d1], arg$add) 
	
	heobjfoll <- function(xfoll, play, d1, d2, arg)
		transfoll(argtest2$heobj, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$foll[d1], arg$foll[d2], arg$add)
	

	if(!is.null(argtest1$constr))
		constrfoll <- function(xfoll, play, arg)
			transfoll(argtest1$constr, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$add)
	else
		constrfoll <- NULL

	if(!is.null(argtest1$grconstr))
		grconstrfoll <- function(xfoll, play, d1, arg)
			transfoll(argtest1$grconstr, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$foll[d1], arg$add)
	else
		grconstrfoll <- NULL

	if(!is.null(argtest2$heconstr))
		heconstrfoll <- function(xfoll, play, d1, d2, arg)
			transfoll(argtest2$heconstr, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$foll[d1], arg$foll[d2], arg$add)
	else
		heconstrfoll <- NULL
	
	if(!is.null(argtest1$joint))
		jointfoll <- function(xfoll, arg)
			transfoll(argtest1$joint, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$add)
	else
		jointfoll <- NULL

	if(!is.null(argtest1$grjoint))
		grjointfoll <- function(xfoll, d1, arg)
			transfoll(argtest1$grjoint, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[d1], arg$add)
	else
		grjointfoll <- NULL
	
	if(!is.null(argtest2$hejoint))
		hejointfoll <- function(xfoll, d1, d2, arg)
			transfoll(argtest2$hejoint, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[d1], arg$foll[d2], arg$add)
	else
		hejointfoll <- NULL
	
	listfollfunc <- list(grobjfoll=grobjfoll, heobjfoll=heobjfoll,
						 constrfoll=constrfoll, grconstrfoll=grconstrfoll, heconstrfoll=heconstrfoll,
						 jointfoll=jointfoll, grjointfoll=grjointfoll, hejointfoll=hejointfoll)
	
#compute the objective of leaders for a strategy x of leaders
#and corresponding followers actions
	objleaders <- function(xlead)
	{
		
#n <- sum(arg1$dimx)
#		nfoll <- sum(arg1$dimx[followers])
		z <- rep(NA, n)
		foll <- bestresponse(xlead, argtest1, argtest2, leaders, followers,
							 index4xfoll, index4lamfoll, index4mufoll, init, 
							 listfollfunc, method.follower=method.follower, 
							 control.follower=control.follower,
							 maxit.follower=maxit.follower, ...)
		
		z[index4xlead] <- xlead
		z[index4xfoll] <- foll$par[1:nfoll]
		
		objleaders <- sapply(leaders, function(i) argtest3$obj(z, i, argtest3$argobj))
		
		if(length(objleaders) > 1)
			return( sqrt(sum(objleaders^2)) )
		else
			return( objleaders )
	}
	
	sapply(x, objleaders)
}	



SE.bestresponse <- function(x, leaders, init, dimx, dimlam, 
	obj, argobj, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	compl, gcompla, gcomplb, argcompl, 
	dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint, 
	method.follower="default", control.follower=list(), maxit.follower, 
	silent=TRUE, ...)
{
	if(method.follower == "default") method.follower <- "Newton"
	if(!is.matrix(init))
		init <- matrix(init, length(x), length(init), byrow=TRUE)
	
	
	argtest1 <- testargfunSSR(init[1,], dimx, dimlam, grobj, arggrobj, constr, argconstr,  grconstr, arggrconstr, 
							  compl, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint)
	
#basic tests for funSSR
	test.try <- try( funSSR(init[1,], dimx, dimlam, grobj, arggrobj, constr, argconstr,  
							grconstr, arggrconstr, compl, argcompl, dimmu, joint, argjoint,
							grjoint, arggrjoint), silent=silent )
	
	if(class(test.try) == "try-error")
	return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate Phi(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
	return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Phi(init) has infinite or NaN values.", fvec=NA) )
	
	argtest2 <- testargjacSSR(init[1,], dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
							  heconstr, argheconstr, gcompla, gcomplb, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint,
							  hejoint, arghejoint)	
	
	
#basic tests for jacSSR
	test.try <- try( jacSSR(init[1,], dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
							heconstr, argheconstr, gcompla, gcomplb, argcompl, dimmu, joint, argjoint,
							grjoint, arggrjoint, hejoint, arghejoint), silent=silent )
	if(class(test.try) == "try-error")
	return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evaluate Jac Phi(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
	return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Jac Phi(init) has infinite or NaN values.", fvec=NA) )
	
	nbplay <- argtest1$nplayer
	
	argtest3 <- testarggapNIR(init[1,1:nbplay], dimx, obj, argobj)
	
	if(!is.numeric(leaders) || length(leaders) > nbplay-1)
	stop("wrong leaders argument.")
	if(any(!leaders %in% 1:nbplay))
	stop("wrong leaders argument.")
	followers <- (1:nbplay)[!(1:nbplay %in% leaders)]
	
	dimx <- argtest1$dimx
	n <- sum(dimx)
	nfoll <- sum(dimx[followers])
	nlead <- sum(dimx[leaders])
	dimlam <- argtest1$dimlam
	m <- sum(dimlam)
#1st row is the begin index, 2nd row the end index
	index4lam <- rbind( cumsum(dimlam) - dimlam + 1, cumsum(dimlam) ) + n
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) ) + 0	
	
	index4xfoll <- as.vector(sapply(followers, function(i) index4x[1,i]:index4x[2,i]))
	index4xlead <- as.vector(sapply(leaders, function(i) index4x[1,i]:index4x[2,i]))
	index4lamfoll <- as.vector(sapply(followers, function(i) index4lam[1,i]:index4lam[2,i]))
	index4mufoll <- (1:length(init))[-(1:(n+m))]
	
	grobjfoll <- function(xfoll, play, d1, arg)
	transfoll(argtest1$grobj, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$foll[d1], arg$add) 
	
	heobjfoll <- function(xfoll, play, d1, d2, arg)
	transfoll(argtest2$heobj, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$foll[d1], arg$foll[d2], arg$add)
	
	
	if(!is.null(argtest1$constr))
		constrfoll <- function(xfoll, play, arg)
			transfoll(argtest1$constr, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$add)
	else
		constrfoll <- NULL
	
	if(!is.null(argtest1$grconstr))
		grconstrfoll <- function(xfoll, play, d1, arg)
			transfoll(argtest1$grconstr, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$foll[d1], arg$add)
	else
		grconstrfoll <- NULL
	
	if(!is.null(argtest2$heconstr))
		heconstrfoll <- function(xfoll, play, d1, d2, arg)
			transfoll(argtest2$heconstr, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[play], arg$foll[d1], arg$foll[d2], arg$add)
	else
		heconstrfoll <- NULL
	
	if(!is.null(argtest1$joint))
		jointfoll <- function(xfoll, arg)
			transfoll(argtest1$joint, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$add)
	else
		jointfoll <- NULL
	
	if(!is.null(argtest1$grjoint))
		grjointfoll <- function(xfoll, d1, arg)
			transfoll(argtest1$grjoint, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[d1], arg$add)
	else
		grjointfoll <- NULL
	
	if(!is.null(argtest2$hejoint))
		hejointfoll <- function(xfoll, d1, d2, arg)
			transfoll(argtest2$hejoint, xfoll, index4xfoll, index4xlead, followers, leaders, 
			  arg, arg$foll[d1], arg$foll[d2], arg$add)
	else
		hejointfoll <- NULL
	
	listfollfunc <- list(grobjfoll=grobjfoll, heobjfoll=heobjfoll,
						 constrfoll=constrfoll, grconstrfoll=grconstrfoll, heconstrfoll=heconstrfoll,
						 jointfoll=jointfoll, grjointfoll=grjointfoll, hejointfoll=hejointfoll)
	
#compute the objective of leaders for a strategy x of leaders
#and corresponding followers actions
	bestresp <- function(i)
	{
		foll <- bestresponse(x[i], argtest1, argtest2, leaders, followers,
							 index4xfoll, index4lamfoll, index4mufoll, init[i,], 
							 listfollfunc, method.follower=method.follower, 
							 control.follower=control.follower,
							 maxit.follower=maxit.follower, ...)
		c(par=foll$par[1:nfoll], code=foll$code, value=foll$value, lagrmult=foll$par[-(1:nfoll)])
	}
	
	sapply(1:length(x), bestresp)
}	

