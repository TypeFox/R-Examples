# This free software is distributed under de GNU GPL version 2 license agreement
# This source code may be copied and modified since the author is cited.
# This software is free and it is provided with no warranty whatsoever
# This software can be cited as follow:
# cite the dissertation the source code is attached:
#
# Junger, W. L. (2004) Modelo Poisson-Gama Semi-Paramétrico: Uma Abordagem de Penalização por Rugosidade. MSc Thesis. Rio de Janeiro, PUC-Rio, Departamento de Engenharia Elétrica
#
# R code for Poisson-Gamma Additive Models (pgam)
# it is part of Washington Junger's MSc. dissertations
# 
# this class of model will handle only Poisson-Gamma family for a while
# in the future it is intended to handle Negative Binomial-Beta family also.
#
# started in 11/06/2003 (dd/mm/yyyy)
# last change: (date first)
# 18/10/2003  up and running (debut) - version 0.1.0
# 19/10/2003  log(1e-7) changed to log(0.1) when y=0 (very troubling in many zeros series) - version 0.1.1
# 21/10/2003  fixed some bugs, deviance to control convergence and partial deviance residuals in smoothing algorithm - version 0.1.2
# 22/10/2003  deviance to control convergence (not good) and deviance partial residuals in smoothing algorithm not working - version 0.1.2
# 24/10/2003  fixed some bugs, trying partial deviance residuals once more (works fine)
# 15/12/2003  fixed residuals degree of freedom: resdf <- n-kp-sum(sdf)-fnz-1
# 17/12/2003  some objects renamed
# 02/01/2004  backfitting(...) function returns the matrix of partial residuals now - version 0.1.3
# 03/01/2004  resource consuming envelope plot implemented now, some bugs fixed - version 0.1.4
# 12/01/2004  function call is returned now - version 0.1.5
# 05/02/2004  minor changes in envelope routine - version 0.1.6
# 12/02/2004  several functions re-written, major changes, finally a library - version 0.2.0
# 18/02/2004  C code linked, functions re-written, faster than never, documentation ok, functions removed, others created for perfomance improvement - version 0.3.0 (release version)
# 19/02/2004  fixed documentation problems - version 0.3.1
# 20/02/2004  one second faster, predict.pgam(), pgam.fit() and pgam.likelihood() are compiled C code now, new examples added to documentation - version 0.3.2
# 21/02/2004  X11() call removed
# 17/03/2004  some bugs fixed, some docs corrected, no more references to dispersion parameter, se of last seasonal factor, periodogram fixed, ... - version 0.3.3
# 15/04/2004  some bugs fixed, compliance with R version 1.9.0, default (preferred) optimization method moved to L-BFGS-B, approximate dispersion parameter back, light docs on envelope, ...
# 16/04/2004  removal of partial residuals options, full non-linear predictor, deals with NA's now, ... - version 0.4.0 beta
# 28/08/2004  minor bugs fixed - version 0.4.0
# 13/02/2005  fixed compiler flag, R version >= 1.9.0, broken missing data support removed and minor bugs fixed - version 0.4.1
# 05/03/2005  several changes in order to share functions with ocdm package and co-author included (Ponce) - v. 0.4.2
# 02/04/2005  fixed model estimated degrees of freedom -  v. 0.4.3
# 12/09/2005  now dealing with missing values, computation of null.deviance, some fixes, estimated degrees of freedom returned, help improvement, approximate significance test of smoothing terms - v. 0.4.4
# 16/06/2009  bug fixes to comply with R 2.9.0 - v. 0.4.7
#
#
# To do list:
# -----------
# * implement the analytical method to estimate information matrix
# * correct estimation of spectrum. done! but it needs some improvement!!!
# * extend to use multiples seasonal factors
# * compute gain for smoothed terms
# * implement an na.replace method
# * se for smooth terms is still missing
# * allow extraction of terms via predict.pgam
# * allow local regression smoothing
#
#
# functions begin here -----------------------

formparser <- function(formula,env=parent.frame())
# reads the model formula and splits it in two new formulae: one of parametric
# terms and another for smoothed terms. 
{
#.pgam.dataset <- get(".pgam.dataset",env=env)
#attach(.pgam.dataset)	# necessary because of f()

formula <- as.formula(formula)
formterms <- terms.formula(formula,specials=c("g","f"))	# s'ed terms is supposed to be smoothed (g) and factorised (f)
modterms <- attr(formterms,"term.labels") # assigns labels of model terms
nterms <- length(modterms)	# number of terms in model
fformula <- NULL
pformula <- NULL	# parametric formula
sformula <- NULL  # smoothing formula
oterm <- NULL  # offset term
offset <- NULL
sdf <- NULL  # smoothing df vector
findex <- NULL
pindex <- NULL
sindex <- NULL
oindex <- NULL
fnames <- NULL
fdata <- NULL
fperiod <- NULL
factorvar <- NULL

response <- attr(formterms,"response") 
if (!(response == 0))
	{
	response <- as.formula(paste("~",as.character(attr(formterms,"variables")[2]),sep=""))
	r <- 1
	}
else
	{
	response <- NULL
	r <- 0
	}
if (!is.null(attr(formterms,"specials")$f))
	findex <- attr(formterms,"specials")$f	# indeces of terms to be factorised in modterms array
if (!is.null(attr(formterms,"specials")$g))
	sindex <- attr(formterms,"specials")$g	# indeces of terms to be smoothed in modterms array
if (!is.null(attr(formterms,"offset")))
	oindex <- attr(formterms,"offset") # offset location

if (nterms)
	for (k in (r+1):(nterms+r))
		if (!(k %in% sindex) && !(k %in% oindex) && !(k %in% findex))
    		pindex <- c(pindex,k)  # gets parametric terms indeces

fterms <- length(findex)	# number of terms in factorised formula
pterms <- length(pindex)	# number of terms in parametric formula
sterms <- length(sindex)	# number of terms in smooth formula

if (fterms)
	{
    for (k in 1:fterms)  # seasonal factors formula - one occurrence only !!!! must change!!!
    	{
    	fed <- eval(parse(text=as.character(attr(formterms,"variables")[findex[k]+1]))) # f extraction
  		factorvar <- as.matrix(eval(parse(text=fed$factor)))
		colnames(factorvar) <- fed$factor
		n <- length(factorvar)
		period <- max(factorvar)
		factordata <- as.data.frame(matrix(NA,n,period))
		factornames <- NULL
		for (k in 1:period)
			{
 		   	factordata[k] <- 1*(factorvar == k)
 	 		factornames <- c(factornames,paste(fed$factor,k,sep="."))
  		  	}
		colnames(factordata) <- factornames
		for (d in 1:(period-1))
        	fformula <- paste(fformula,factornames[d],sep=ifelse(is.null(fformula),"~","+"))
    	fnames <- c(fnames,factornames)
    	fdata <- factordata  #c(fdata,factordata)
        fperiod <- c(fperiod,period)
		}
	fformula <- as.formula(fformula)
    fdata <- as.matrix(as.data.frame(fdata))
    }
    
if (pterms)
	{
	for (k in 1:pterms)  # parametric formula
		pformula <- paste(pformula,as.character(attr(formterms,"variables")[pindex[k]+1]),sep=ifelse(is.null(pformula),"~","+"))
	pformula <- as.formula(pformula)
    }

if (sterms)
	{
    for (k in 1:sterms)  # non-parametric formula
    	{
    	sed <- eval(parse(text=as.character(attr(formterms,"variables")[sindex[k]+1]))) # g extraction
		sformula <- paste(sformula,sed$var,sep=ifelse(is.null(sformula),"~","+"))
    	sdf <- c(sdf,sed$df)
		}
	sformula <- as.formula(sformula)
    }

if (!is.null(oindex))  # dealing with the offset term
	{
	oterm <- as.character(attr(formterms,"variables")[[oindex+1]])[2]
	# oterm <- names(attr(formterms,"factors")[attr(formterms,"offset")]) 
	# oterm <- substr(oterm <- substr(oterm,8,nchar(oterm)),1,nchar(oterm)-1)
	offset<- eval(parse(text=oterm))
	}
	
if (attr(formterms,"intercept") == 0)
    intercept <- FALSE
else
	intercept <- TRUE

#detach(.pgam.dataset)  # i don't like it!!!
#rm(.pgam.dataset) 
retval <- list(response=response,pformula=pformula,pterms=pterms,sformula=sformula,sdf=sdf,sterms=sterms,fterms=fterms,fformula=fformula,factor=factorvar,fnames=fnames,fdata=fdata,fperiod=fperiod,oterm=oterm,offset=offset,intercept=intercept,fullformula=formula)
class(retval) <- "split.formula"
return(retval)
}


pgam.filter <- function(w,y,eta)
# runs appropriate recursions for omega estimate
{
n <- length(y)

oo <- .C("pgam_filter",
		as.double(w),
		as.double(y),
		as.integer(n),
		as.double(eta),
		att1=double(n),
		btt1=double(n),
		PACKAGE="pgam")

return(list(att1=oo$att1,btt1=oo$btt1))
}


pgam.likelihood <- function(par,y,x,offset,fperiod,env=parent.frame())
# likelihood function to be maximized
{
# splitting parameter vector into omega and beta
fnz <- fnz(y)
n <- length(y)
psi <- pgam.par2psi(par,fperiod)

#print(psi$beta)  # for debugging

if (!is.null(psi$beta))
	eta <- x%*%psi$beta+offset  # parametric piece of predictor
else
	eta <- offset

filtered <- pgam.filter(psi$w,y,eta)
att1 <- filtered$att1
btt1 <- filtered$btt1

oo <- .C("pgam_loglik",
		as.double(y),
		as.integer(n),
		as.integer(fnz),
		as.double(att1),
		as.double(btt1),
		lvalue=double(1),
		NAOK=TRUE,
		PACKAGE="pgam")

loglik <- list(value=oo$lvalue,eta=eta,att1=att1,btt1=btt1)
assign("loglik",loglik,env=env)

return(oo$lvalue)
}


pgam.fit <- function(w,y,eta,partial.resid)
# estimates of y_hat_t|t-1
{
filtered <- pgam.filter(w,y,eta)
att1 <- filtered$att1
btt1 <- filtered$btt1

ynz <- y+0.1*(y==0)  # replaces zero counts by a small value (dangerous!!!)
n <- length(y)
fnz <- fnz(y)

oo <- .C("pgam_predict",
		as.double(y),
		as.integer(n),
		as.double(att1),
		as.double(btt1),
		yhat=double(n),
		vyhat=double(n),
		deviance=double(n),
		pearson=double(n),
		PACKAGE="pgam")
  
# partial residuals
if (partial.resid == "response")
   	presid <- (log(ynz)-log(oo$yhat))[(fnz+1):n]
#else if (partial.resid == "deviance")
#   	presid <- log(sign(y-oo$yhat)*sqrt(oo$deviance))[(fnz+1):n]
#else if (partial.resid == "pearson")
#	presid <- log((y-oo$yhat)/oo$vyhat)[(fnz+1):n]

return(list(yhat=oo$yhat,resid=presid))
}


pgam.psi2par <- function(w,beta,fperiod)
# puts hyperparameters into optmization form
{
alpha <- log(w/(1-w)) # transformation to ensure constraints on w
fp <- length(fperiod)
newbeta <- NULL

if (!is.null(beta))
	{
	fp <- length(fperiod)
    bk <- 1
    
    if (!is.null(fperiod))
    	{
        for (k in 1:fp)
    		{
	        newbetak <- beta[bk:(fperiod[k]-1)]
        	bk <- length(newbetak)+2  # skip the last seasonal factor ex.: saturday, december etc
        	newbeta <- c(newbeta,newbetak)
        	}
    	if (length(beta) > sum(fperiod))
    		newbeta <- c(newbeta,beta[bk:length(beta)])
		}
	else
    	newbeta <- beta
    }
par <- c(alpha,newbeta) # parameter vector to be optimized
return(par)
}


pgam.par2psi <- function(par,fperiod)
# puts parameters optmization form back into beta
{
alpha <- par[1]
w <- exp(alpha)/(1+exp(alpha))
beta <- NULL

if (length(par) > 1)
	{
	newbeta <- par[2:length(par)]
	fp <- length(fperiod)
	bk <- 1

    if (!is.null(fperiod))
    	{
		for (k in 1:fp)
    		{
        	newbetak <- newbeta[bk:(bk+fperiod[k]-2)]  # get pieces of beta from par
        	bk <- length(newbetak)+1
        	beta <- c(beta,newbetak,-sum(newbetak))
        	}
	    if (length(newbeta) > sum(fperiod))  # something odd is going on about here!
		    beta <- c(beta,newbeta[bk:length(newbeta)])
		}
	else
    	beta <- newbeta
    }
retval <- list(w=w,beta=beta)
return(retval)
}


pgam.hes2se <- function(hes,fperiod,se.estimation="numerical")
# puts parameters hessian matrix in the form of beta s.e.
{
sigma <- try(solve(-hes))

alpha <- sigma[1,1]
se.w <- sqrt(alpha*(exp(alpha)/(1+exp(alpha))^2)^2)  # correcting sigma^2_{omega} using delta rule James, B. (1996), Probabilidade:..., p.253

beta <- NULL
se.beta <- NULL

if (length(diag(sigma)) > 1)
	{
	p <- length(diag(sigma))
	
	sigmabeta <- sigma[2:p,2:p]
	if (p == 2)
		varbeta <- sigmabeta
	else
		varbeta <- diag(sigmabeta)
	
	if (!is.null(fperiod))
    	{
		fp <- length(fperiod)
		bk <- 1
   
		for (k in 1:fp)
    		{
        	varbetak <- varbeta[bk:(bk+fperiod[k]-2)]  # get pieces of beta from par
        	covarbetak <- sigmabeta[bk:(bk+fperiod[k]-2),bk:(bk+fperiod[k]-2)]
			bk <- length(varbetak)+1
        	diag(covarbetak) <- 0.0
			beta <- c(beta,varbetak,sum(varbetak)-sum(covarbetak))
        	}
	    if (length(varbeta) > sum(fperiod))  # something odd is going on about here!
		    beta <- c(beta,varbeta[bk:length(varbeta)])
		}
	else
		{
    	beta <- varbeta
		}
	se.beta <- sqrt(beta)  # extracts standard error
    }
	
retval <- list(w=se.w,beta=se.beta)
return(retval)
}


pgam <- function(formula,dataset,omega=0.8,beta=0.1,offset=1,digits=getOption("digits"),na.action="na.exclude",maxit=1e2,eps=1e-6,lfn.scale=1,control=list(),optim.method="L-BFGS-B",bkf.eps=1e-3,bkf.maxit=1e2,se.estimation="numerical",verbose=TRUE)
# estimates Poisson-Gamma Additive Models
{
st <- proc.time()
called <- match.call()
pgam.env <- new.env()  # new environment defined for pgam
#assign(".pgam.dataset",dataset,env=pgam.env)	# it is necessary because of f(factor)

if (is.null(formula))
	stop("Model formula is not specified.")
if (is.null(dataset))
	stop("Model dataset is not specified.")

# some fixed default options
partial.resid <- "response"
smoother <- "spline"
	
if (!verbose)
	{
	control$trace <- 0
	control$REPORT <- maxit
	}
control$fnscale <- lfn.scale*(-1)

parsed <- formparser(formula,env=pgam.env)  # parse the formula and get components
# building datasets and setting the model structure (no explanatory variables, seasonal factors, linear, additive, offset, drift)
y <- framebuilder(parsed$response,dataset)
yname <- names(y)
y <- as.matrix(y)
n.obs <- dim(y)[1]
etap <- NULL
px <- NULL
kx <- 0
pnames <- NULL
pform <- NULL
if (parsed$pterms)
	{
    pform <- parsed$pformula
	px <- framebuilder(pform,dataset)
    pnames <- names(px)
    px <- as.matrix(px)
    kx <- dim(px)[2]  # number of non-seasonal factor explanatory variables
    }

etas <- NULL
bkf <- NULL
sx <- NULL
sdf <- NULL
ks <- 0
snames <- NULL
sform <- NULL
sduplicated <- NULL
var.smox <- NULL
if (parsed$sterms)
	{
    sform <- parsed$sformula
	sx <- framebuilder(sform,dataset)
    snames <- names(sx)
    sx <- as.matrix(sx)
	sdf <- parsed$sdf
	ks <- dim(sx)[2] # number of non-parametric explanatory variables
	# getting duplicated values positions
	sduplicated <- apply(sx,2,sdup <- function(x){duplicated(sort(signif(x,6)))})	
	}

fx <- NULL
fp <- 0
kf <- 0
fperiod <- NULL
fnames <- NULL
fform <- NULL
if (parsed$fterms)
	{
    fform <- parsed$fformula
	fx <- parsed$fdata
    fnames <- parsed$fnames
    fperiod <- parsed$fperiod
    fp <- length(fperiod)  # number of seasonal factors terms
    kf <- sum(fperiod) # total count of seasonal factors
    }

offcoef <- offset
if (!is.null(parsed$oterm))
    {
	offset <- as.matrix(parsed$offset*offcoef)
	colnames(offset) <- parsed$oterm
	}
else
	offset <- double(n.obs)

if (parsed$intercept)
	{
    # to be worked out
    }

terms <- list(response=parsed$response,pform=pform,fform=fform,sform=sform,snames=snames,df=sdf)

px <- cbind(fx,px)
kp <- kf+kx # number of explanatory variables

# print(dim(sx));cat("\n");print(dim(px));cat("\n")
# getting rid of NA's
tempDataset <- cbind(y,offset,px,sx)
#	stop("Missing data is not handled for now. It will be fixed soon.")
tempDataset <- eval(parse(text=paste(na.action,"(tempDataset)",sep="")))
if (!is.null(attr(tempDataset,"na.action")))
	na.info <- list("na.action"=attr(tempDataset,"na.action"),"na.class"=attr(attr(tempDataset,"na.action"),"class"))
else
	na.info <- list("na.action"=attr(tempDataset,"na.action"),"na.class"=substr(na.action,4,nchar(na.action)))
#print(dim(tempDataset));cat("\n")
y <- as.matrix(tempDataset[,1])
n <- length(y)  # get number of valid obs
fnz <- fnz(y)
offset <- as.matrix(tempDataset[,2])
if (!is.null(px))
	px <- as.matrix(tempDataset[,3:(kp+2)])
if (!is.null(sx))
	sx <- as.matrix(tempDataset[(fnz+1):n,(kp+3):(kp+ks+2)])
#print(dim(sx));cat("\n");print(dim(px));cat("\n")

undef <- rep(0,fnz)

if (sum((y-as.integer(y))>0))
	{
    y <- as.integer(y) # non-integer values are truncated
    warning("Some values must have been truncated in order to ensure compliance with Poisson specification.")
    }
else
    y <- as.integer(y) # to avoid type conflict with integer functions
	
if (sum(y < 0) > 0)	# checking for negative values. if detected, program halts
	stop("Negative observations detected. The series does not comply with Poisson specification.")

if (length(beta) == 1)
	beta <- rep(beta,kp) # expansion of initial value of beta vector (assuming all the same)
if (kp == 0)
	beta <- NULL
#if (is.null(px) && !is.null(sx))
#   	stop("This application still not fit full non-linear predictor models.")

assign("loglik",NULL,env=pgam.env)

# mle estimation
par <- pgam.psi2par(omega,beta,fperiod)

optimized <- optim(par,pgam.likelihood,y=y,x=px,offset=offset,fperiod=fperiod,env=pgam.env,method=optim.method,hessian=FALSE,control=control)
newoffset <- offset  # case of full parametric model

# smoothing
k <- 0  # if k=0 at the end o function, this is a full parametric model
norm <-0  # same as k
if (!is.null(sx))
	{
    norm <- 1e35  # large number intialization of norm
    # oldetas <- 0
    loglik0 <- 0
    # deviance0 <- 0
    k <- 0
    for (k in 1:maxit)
     	{
        # getting useful information
		psi <- pgam.par2psi(optimized$par,fperiod)
        if (!is.null(beta))
			etap <- px%*%psi$beta  # parametric piece of predictor
		else
			etap <- double(n)
		presid <- pgam.fit(psi$w,y,etap+offset,partial.resid)$resid
		loglik1 <- (-1)*optimized$value
        # backfitting smoothing
        bkf <- backfitting(presid,sx,sdf,smoother=smoother,eps=bkf.eps,maxit=bkf.maxit,info=FALSE)
        etas <- c(undef,bkf$sumfx)
      	newoffset <- offset+etas  # preserves original offset in full eta
        norm <- abs((loglik1-loglik0)/loglik0)
		if (verbose)
			cat(paste("iter:",k," | ","norm:",round(norm,digits),"\n"))

        if (!(norm < eps))
            {
	    	if (k > maxit)
            	{
            	if (verbose)
					cat("\nNo convergence after last iteration of the estimation algorithm.\n")
                break
            	}
            # oldetas <- etas
            loglik0 <- loglik1
            # parametric fitting
            optimized <- optim(optimized$par,pgam.likelihood,y=y,x=px,offset=newoffset,fperiod=fperiod,env=pgam.env,method=optim.method,hessian=FALSE,control=control)
            }
         else
		 	{
			if (verbose)
				cat("\nSemiparametric model estimation algorithm has converged.\n")
            break
			}
        }
	}

# debugging
#mypsi<-pgam.par2psi(optimized$par,fperiod)
#print(mypsi)

# last running in order to get evaluated functions and hessian matrix
if (verbose)
	cat("\nFinal run: Getting estimated parameters, functions and numerical hessian matrix...\n")
numerical.se <- ifelse(se.estimation=="numerical",TRUE,FALSE)
optimized <- optim(optimized$par,pgam.likelihood,y=y,x=px,offset=newoffset,fperiod=fperiod,env=pgam.env,method=optim.method,hessian=numerical.se,control=control)
psi <- pgam.par2psi(optimized$par,fperiod)
if (!is.null(sx))
	{
    if (!is.null(beta))
		etap <- px%*%psi$beta  # parametric piece of predictor
	else
		etap <- double(n)
    presid <- pgam.fit(psi$w,y,etap+offset,partial.resid)$resid
    # backfitting smoothing
    bkf <- backfitting(presid,sx,sdf,smoother=smoother,eps=bkf.eps,maxit=bkf.maxit,info=TRUE)
    # restoring original size of elements in bkf and sx
	bkf$sumfx <- c(undef,bkf$sumfx)
    bkf$smox <- rbind(matrix(NA,fnz,ks),bkf$smox)
	dimnames(bkf$smox) <- list(NULL,snames)
    bkf$pres <- rbind(matrix(NA,fnz,ks),bkf$pres)
	dimnames(bkf$pres) <- list(NULL,snames)
	sx <- rbind(matrix(NA,fnz,ks),sx)
    # computing variance of smoothing functions
	sigma.smox <- apply((bkf$pres-bkf$smox)^2,2,sum,na.rm=TRUE)/(n-fnz-bkf$edf)
	var.smox <- matrix(NA,n,length(sigma.smox))	
	for (j in 1:length(sigma.smox))
		{
		var.smox.short <- sigma.smox[j]*bkf$lev[,j]^2
		hold <- 0
		for (i in 1:n)
			{
			if (sduplicated[i,j])
				hold <- hold+1
			var.smox[i,j] <- var.smox.short[i-hold]
			}
		}
	bkf$var.smox <- var.smox
	}

# getting useful information
omega <- psi$w
names(omega) <- c("(Discount)")
beta <- psi$beta
names(beta) <- c(fnames,pnames)
se.omega <- NA
se.beta <- NA
if (numerical.se)
	{
    se <- pgam.hes2se(optimized$hessian,fperiod)
    se.omega <- se$w
    names(se.omega) <- c("(Discount)")
    se.beta <- se$beta
    names(se.beta) <- c(fnames,pnames)
    }

iterations <- optimized$counts
if (verbose)
	{cat("Counts (fn | gr): ");cat(iterations);cat("\n")}

if (!is.null(optimized$convergence))
	if (optimized$convergence == 0)
		{
        convergence <- "converged"
        }
else
	convergence <- "not converged"

loglik.value <- (-1)*optimized$value
loglik <- get("loglik",env=pgam.env)
att1 <- loglik$att1
btt1 <- loglik$btt1

# corrected estimated residuals degrees of freedom (including ^w)
if (!is.null(bkf))
	edf <- length(beta)+fnz+sum(bkf$edf+1)  # correction suggested by Hastie and Tibshirani (1990)
else
	edf <- length(beta)+fnz+1
resdf <- n-edf

# including seasonal factors in returned dataset
if (parsed$fterms)
	{factor <- parsed$factor; factor[na.info$na.action] <- NA; factor <- na.omit(factor)}
else
	factor <- NULL
dataset <- cbind(tempDataset,factor)  # deparse(substitute(dataset))

et <- proc.time()
if (verbose)
	cat(paste("\nEstimation process took", elapsedtime(st,et)),"(hh:mm:ss)\n")

retval <- list(call=called,formula=formula,dataset=dataset,omega=omega,se.omega=se.omega,beta=beta,se.beta=se.beta,att1=att1,btt1=btt1,loglik=loglik.value,convergence=convergence,optim.method=optim.method,y=y,px=px,sx=sx,terms=terms,offset=offset,offcoef=offcoef,etap=etap,etas=etas,partial.resid=partial.resid,bkf=bkf,alg.k=k,alg.norm=norm,edf=edf,resdf=resdf,n.obs=n.obs,n=n,tau=fnz,na.info=na.info)
class(retval) <- "pgam"
return(retval)
}


backfitting <- function(y,x,df,smoother="spline",w=rep(1,length(y)),eps=1e-3,maxit=1e2,info=TRUE)
# ajust the non-parametric piece of predictor
{
# getting useful information
n <- dim(x)[1]	# number of observations
p <- dim(x)[2]	# number of variables to be smoothed

# initialization
smox <- matrix(0,n,p)  # create storage space for smoothed variables
pres <- matrix(0,n,p)  # create storage space for partial residuals variables
lev <- matrix(NA,n,p)  # create storage space for smoother leverage
edf <- double(p)  # create storage space for edf
norm <- 1e35	# large value
smooth <- double(n)
m <- 0

#backfitting
while ((norm >= eps) && (m <= maxit))
	{
	m <- m+1
    oldsmooth <- smooth
    smooth <- double(n)
    for (j in 1:p)
		{
		sumfx <- double(n)
        for (k in 1:p)
            smox[,k] <- (k!=j)*bkfsmooth(y,x[,k],df[k],smoother=smoother,w=w)$fitted
        for (k in 1:p)
            sumfx <- sumfx+(k!=j)*smox[,k]
        pres[,j] <- y-sumfx
        smoothed <- bkfsmooth(pres[,j],x[,j],df[j],smoother=smoother,w=w)
		smox[,j] <- smoothed$fitted
		lev[1:length(smoothed$lev),j] <- smoothed$lev
		edf[j] <- smoothed$df
        smooth <- smooth+smox[,j]
        }
    
	norm <- lpnorm(smooth,oldsmooth,p=2)/lpnorm(oldsmooth,p=2)
	}

sumfx <- apply(smox,1,sum)  # sum of fx
if (!info)
	{
	# short returning list - intended to be fast!
	retval <- list(sumfx=sumfx)
	}
else
	{
	# complete returning list - slower!
	retval <- list(sumfx=sumfx,smox=smox,pres=pres,lev=lev,edf=edf)
	}
class(retval) <- "backfitting"
return(retval)
}


bkfsmooth <- function(y,x,df,smoother="spline",w=rep(1,length(y)))
# smooths y against x wiht ds degrees of smoothness
{
if (smoother=="spline")
	{
	smoothed <- smooth.spline(x=x,y=y,w=w,df=df)
	fitted <- predict(smoothed,x)$y
	lev <- smoothed$lev
	df <- smoothed$df
	}
#else if (smoother=="loess")
#	{
#	tempds <- as.data.frame(cbind(y,x))	# temporary dataset
#	names(tempds) <- c("y","x")
#	fitted <- loess(as.formula("y~x"),tempds,weights=w,span=1/df)$fitted
#   retval <- list(fitted=fitted)
#	}
else
	stop(paste("Smoother",smoother,"not implemented yet."))

retval <- list(fitted=fitted,lev=lev,df=df)    
return(retval)
}


residuals.pgam <- function(object,type="deviance",...)
# method for residuals extraction of a pgam model object
{
predicted <- predict(object,...)
na.action <- object$na.info$na.action

# adopting GLM notation for residuals extraction
y <- napredict(na.action,object$y)
mu <- predicted$yhat
v.mu <- predicted$vyhat
d <- predicted$deviance
h <- predicted$hat
att1 <- napredict(na.action,object$att1)
btt1 <- napredict(na.action,object$btt1)
phi <- predicted$scale

raw <- y-mu  # getting raw residuals

if (type == "response")
	resid <- raw
else if (type == "pearson")
	resid <- raw/sqrt(v.mu)
else if (type == "deviance")
	resid <- sign(raw)*sqrt(d)
else if (type == "std_deviance")
	resid <- sign(raw)*sqrt(d/(1-h))
else if (type == "std_scl_deviance")
	resid <- sign(raw)*sqrt(d/(phi*(1-h)))
#else if (type == "adj_deviance")
#	{
#    skew <- (1/6)*(2-att1)/sqrt(btt1*(1-att1))  # skewness coeficient of negative binomial distribution
#	resid <- sign(raw)*sqrt(d)*skew
#   }
else
	stop(paste("Residuals of type",type,"is not implmented yet."))

retval <- naresid(na.action,resid)  # put NA back
return(retval)
}


predict.pgam <- function(object,forecast=FALSE,k=1,x=NULL,...)
# method for prediction for new values
{
# getting goods
n <- object$n
y <- object$y
etap <- object$etap
etas <- object$etas
w <- object$omega
na.action <- object$na.info$na.action

# estimates of y_hat_t|t-1
if (is.null(etap) && !(is.null(etas)))
	{
	model <- "nonparametric"
	etap <- double(n)
	}
else if (!is.null(etap) && is.null(etas))
	{
	model <- "fullparametric"
	etas <- double(n)
	}
else if (is.null(etap) && is.null(etas))
	{
	model <- "levelonly"
	etap <- double(n)
	etas <- double(n)
	}
else
	model <- "semiparametric"

eta <- etap+etas

filtered <- pgam.filter(w,y,eta)
att1 <- filtered$att1
btt1 <- filtered$btt1

n <- length(y)
hat <- double(n)
fnz <- fnz(y)
level <- NULL
yhats <- NULL

oo <- .C("pgam_predict",
		as.double(y),
		as.integer(n),
		as.double(att1),
		as.double(btt1),
		yhat=double(n),
		vyhat=double(n),
		deviance=double(n),
		pearson=double(n),
		PACKAGE="pgam")

for (t in 1:(n-1))
	{
	# pseudo hat matrix
   	hat[t] <- w*exp(eta[t+1])/sum(w^(0:(t-1))*exp(eta[t:1]),na.rm=TRUE)
	}

# an attempt of extraction of components
if (model == "semiparametric")
	{
    # semiparametric model
    plevel <- oo$yhat*exp(-etap)
	level <- plevel*exp(-etas)
	yhats <- plevel/level
    }
else if (model == "fullparametric")
	{
    # full parametric model
    level <- oo$yhat*exp(-etap)
    yhats <- NULL
    }
else if (model == "levelonly")
	level <- oo$yhat

# residuals are to be extracted elsewhere

# scale parameter based on generalized Pearson statitics
scale <- sum(oo$pearson,na.rm=TRUE)/object$resdf

# cleaning the house
oo$yhat[1:fnz] <- NA
oo$vyhat[1:fnz] <- NA
oo$deviance[1:fnz] <- NA
oo$pearson[1:fnz] <- NA

# put NA back
oo$yhat <- napredict(na.action,oo$yhat)
oo$vyhat <- napredict(na.action,oo$vyhat)
oo$deviance <- napredict(na.action,oo$deviance)
oo$pearson <- napredict(na.action,oo$pearson)
level <- napredict(na.action,level)
yhats <- napredict(na.action,yhats)
    
# forecasting
if (forecast)
	{
	stop("Sorry! Forecasting is broken for now.\n")
	if (model == "semiparametric")
		stop("Forecasting of semiparametric models is not implemented yet.\n")
	fc <- double(k)
	if (model == "levelonly")
		for (l in 1:k)
			fc[l] <- sum(w^(0:(n-1))*y[n:1],na.rm=TRUE)/sum(w^(0:(n-1)),na.rm=TRUE)
	else
		{
		if (!is.null(x))
			{
			etaf <- x%*%beta
			for (t in 1:k)
				fc[l] <- w*exp(etaf[n+l])*sum(w^(0:(n-1))*y[n:1],na.rm=TRUE)/sum(w^(0:(n-1))*exp(eta[n:1]),na.rm=TRUE)
			}
		else
			{
			cat("\nCovariates were not supplied. Forecasting is not possible.\n")
			fc <- NULL
			}
		}	
	forecast <- fc
	}

retval <- list(yhat=oo$yhat,vyhat=oo$vyhat,deviance=oo$deviance,pearson=oo$pearson,scale=scale,hat=hat,level=level,yhats=yhats,forecast=forecast)
return(retval)
}


fitted.pgam <- function(object,...)
# method for fitted extraction only
{
predicted <- predict(object,...)

retval <- predicted$yhat
return(retval)
}


coef.pgam <- function(object,...)
# method for parametric coefficients extraction. although omega is not a coef, it is output in the first slot.
{
retval <- c(object$omega,object$beta)
return(retval)
}


logLik.pgam <- function(object,...)
# logLik method for extraction of loglik from a pgam object
{
retval <- object$loglik
return(retval)
}


deviance.pgam <- function(object,...)
# deviance method for extraction of deviance from a pgam object
{
predicted <- predict(object,...)

retval <- sum(predicted$deviance,na.rm=TRUE)
return(retval)
}


AIC.pgam <- function(object,k=2,...)
# method for AIC estimation of a pgam model object
{
predicted <- predict(object,...)

# gathering necessary quantities
nprime <- object$n-object$tau
deviance <- sum(predicted$deviance,na.rm=TRUE)
edf <- object$edf
phi <- predicted$scale

retval <- (1/nprime)*(deviance+k*edf*phi)
# retval <- (1/nprime)*(deviance+k*edf)
return(retval)
}


summary.pgam <- function(object,smo.test=FALSE,...)
# method for output summary of pgam model object
{
predicted <- predict(object,...)

# general stuff
call <- object$call
formula <- object$formula
convergence <- object$convergence
optim.method <- object$optim.method
n.obs <- object$n.obs
n <- object$n
tau <- object$tau
alg.k <- object$alg.k
alg.norm <- object$alg.norm
fterms <- terms.formula(formula,specials=c("g","f"))
# get null.deviance
if (!any(attr(fterms,"offset")))
	nullform <- paste(as.character(formula)[2]," ~ NULL",sep="")
else
	nullform <- paste(as.character(formula)[2]," ~ NULL + offset(",as.character(attr(fterms,"variables")[[attr(fterms,"offset")+1]])[2],")",sep="")
null.model <- pgam(formula=nullform,omega=object$omega,offset=object$offcoef,dataset=as.data.frame(object$dataset),na.action=paste("na.",object$na.info$na.class,sep=""),optim.method=object$optim.method,se.estimation="none",verbose=FALSE)
null.deviance <- deviance.pgam(null.model)

# goodness-of-fit statistics
deviance <- sum(predicted$deviance,na.rm=TRUE)
pearson <- sum(predicted$pearson,na.rm=TRUE)
scale <- predicted$scale
edf <- object$edf
res.edf <- object$resdf
# maximum likelihood parameters and hypothesis testing
loglik <- object$loglik
coeff <- c(object$omega,object$beta)
se.coeff <- c(object$se.omega,object$se.beta)
t.coeff <- coeff/se.coeff
pt.coeff <- 2*(1-pt(abs(t.coeff),df=res.edf))
if (!is.null(object$bkf) && (smo.test))
	{
	# nonparametric stuff
	vars.smo <- dimnames(object$bkf$smox)[[2]]
	nvars.smo <- length(dimnames(object$bkf$smox)[[2]])
	edf.smo <- object$bkf$edf
	test.dev <- double(nvars.smo)
	# buiding parametric part formula
	baseform <- paste(as.character(formula)[2]," ~ f(", as.character(attr(fterms,"variables")[[attr(fterms,"specials")$f+1]])[2],") + ", as.character(object$terms$pform)[2],sep="")
	for (j in 1:nvars.smo)
		{
		devform <- paste(baseform," + ",vars.smo[j],sep="")
		vars.smo.j <- vars.smo; vars.smo.j[j] <- NA; vars.smo.j <- na.omit(vars.smo.j)
		edf.smo.j <- object$terms$df; edf.smo.j[j] <- NA; edf.smo.j <- na.omit(edf.smo.j)
		devform <- paste(devform," + g(",vars.smo.j,",",edf.smo.j,")",sep="")
		dev.mod <- 
		pgam(formula=devform,omega=object$omega,beta=c(object$beta,mean(object$beta)),offset=object$offcoef,dataset=as.data.frame(object$dataset),na.action=paste("na.",object$na.info$na.class,sep=""),optim.method=object$optim.method,se.estimation="none",verbose=FALSE)
		test.dev[j] <- deviance.pgam(dev.mod)
		}
	# approximate F-tests must be inserted at this point (soon!)
	chi.smo <- abs(deviance-test.dev)
	pchi.smo <- 1-pchisq(chi.smo,edf.smo-1)
	}
else if (!is.null(object$bkf))
	{
	# nonparametric stuff
	vars.smo <- dimnames(object$bkf$smox)[[2]]
	edf.smo <- object$bkf$edf
	chi.smo <- NULL
	pchi.smo <- NULL
	}
else
	{
	vars.smo <- NULL
	edf.smo <- NULL
	chi.smo <- NULL
	pchi.smo <- NULL
	}

retval <- list(call=call,formula=formula,convergence=convergence,optim.method=optim.method,n.obs=n.obs,n=n,tau=tau,alg.k=alg.k,alg.norm=alg.norm,deviance=deviance,pearson=pearson,edf=edf,res.edf=res.edf,scale=scale,loglik=loglik,null.deviance=null.deviance,coeff=coeff,se.coeff=se.coeff,t.coeff=t.coeff,pt.coeff=pt.coeff,vars.smo=vars.smo,edf.smo=edf.smo,chi.smo=chi.smo,pchi.smo=pchi.smo,smo.test=smo.test)
class(retval) <- "summary.pgam"
return(retval)
}


print.summary.pgam <- function(x,digits=getOption("digits"),...)
# method for summary printing
{
cat("Function call:\n")
print(x$call)
cat("\nModel formula:\n")
print(x$formula)
cat("\nParametric coefficients:\n")
width <- max(nchar(names(x$coeff)))+2
cat(rep(" ",width),"     Estimate     std. err.     t ratio      Pr(>|t|)\n",sep="")
for (i in 1:length(x$coeff))
    cat(formatC(names(x$coeff)[i],width=width)," ",formatC(x$coeff[i],width=digits+6,digits=digits)," ",formatC(x$se.coeff[i],width=digits+6,digits=digits)," ",formatC(x$t.coeff[i],width=digits+6,digits=digits),"    ",format.pval(x$pt.coeff[i]),"\n",sep="")
# nonparametric partition of the model
if (!is.null(x$vars.smo) && (x$smo.test))
	{
	cat("\nApproximate significance of nonparametric terms:\n")
    width <- max(nchar(x$vars.smo))+2
    cat(rep(" ",width),"       edf            chi.sq         p-value\n",sep="")
    for (i in 1:length(x$vars.smo))
		cat(formatC(x$vars.smo[i],width=width)," ",formatC(x$edf.smo[i]-1,width=digits+6,digits=digits),"   ",formatC(x$chi.smo[i],width=digits+6,digits=digits),"     ",format.pval(x$pchi.smo[i]),"\n",sep="")
	}
else if (!is.null(x$vars.smo))
	{
	cat("\nNonparametric terms:\n")
    width <- max(nchar(x$vars.smo))+2
    cat(rep(" ",width),"       edf\n",sep="")
    for (i in 1:length(x$vars.smo))
		cat(formatC(x$vars.smo[i],width=width)," ",formatC(x$edf.smo[i]-1,width=digits+6,digits=digits),"   "
		,"\n",sep="")
	}

cat("\nLog-likelihood value is ",round(x$loglik,digits)," after ",x$optim.method," has ",x$convergence,".\n",sep="")
if (x$alg.k > 0)
	cat("\nEstimation process of semiparametric model stopped after ",x$alg.k," iterations at ",round(x$alg.norm,digits),".\n",sep="")
else
	cat("\nFull parametric model estimated.\n")
cat("\nNull deviance is ",round(x$null.deviance,digits)," on ",x$n-x$tau-1," degrees of freedom.\n",sep="")
cat("\nResidual deviance is ",round(x$deviance,digits)," on ",x$res.edf," degrees of freedom.\n",sep="")
cat("\nApproximate dispersion parameter equals ",round(x$scale,digits)," based on generalized Pearson statistics ",round(x$pearson,digits),".\n",sep="")
# cat("\nGeneralized Pearson statistics is ",round(x$pearson,digits),".\n",sep="")
cat("\nDiffuse initialization wasted the ",x$tau," first observation(s). \nIn addition, ",x$n.obs-x$n," observation(s) lost due to missingness.\n",sep="")
}


print.pgam <- function(x,digits=getOption("digits"),...)
{
cat("Function call:\n")
print(x$call)
cat("\nModel formula:\n")
print(x$formula)
cat("\nLog-likelihood value is ",round(x$loglik,digits)," after ",x$optim.method," has ",x$convergence,".\n",sep="")
if (x$alg.k > 0)
	cat("\nEstimation process of semiparametric model stopped after ",x$alg.k," iterations at ",round(x$alg.norm,digits),".\n",sep="")
else
	cat("\nFull parametric model estimated.\n")
cat("\nDiffuse initialization wasted the ",x$tau," first observation(s). \nIn addition, ",x$n.obs-x$n," observation(s) lost due to missingness.\n",sep="")
invisible(x)
}


plot.pgam <- function(x,rug=TRUE,se=TRUE,at.once=FALSE,scaled=FALSE,...)
# method for smooth terms plotting
{
predicted <- predict(x,...)

#gathering some useful information
na.action <- x$na.info$na.action
y.level <- predicted$level
x.level <- seq(1:length(y.level))
if (!is.null(x$bkf$edf))
	edf <- round(x$bkf$edf,0)  # rounding for better visualisation
else
	edf <- x$bkf$edf
# se <- napredict(na.action,sqrt(x$bkf$var.smox))

#warn.opt <- getOption("warn"); options(warn=-1)
# plotting level
plot(x.level,y.level,type="l",xlab="time",ylab="local level",...)
if (rug)
	rug(x.level)

# plotting smoothed covariates in predictor scale
if (!is.null(edf))
	{
	s <- length(edf)
	vars.smo <- dimnames(x$bkf$smox)[[2]]
	x.g <- napredict(na.action,x$sx)
	y.g <- napredict(na.action,x$bkf$smox)
	
#	lb.se <- y.g-2*se
#	ub.se <- y.g+2*se
	#setting labels
	y.g.lab <- paste("g(",vars.smo,",",edf,")",sep="")
	for (i in 1:s)
		{
		if (at.once)
			getOption("device")()
		else
			if (interactive())
				{
				prompt <- readline("Press ENTER for next page or X to exit and keep this page... ")
				if((prompt == "x") || (prompt == "X"))
					break
				}
		
		# setting scale to smoothed covariates plots
#		y.g.lim <- c(min(lb.se[,i],na.rm=TRUE),max(ub.se[,i],na.rm=TRUE))
		if (scaled)
			y.g.lim <- c(min(y.g,na.rm=TRUE),max(y.g,na.rm=TRUE))
		else
			y.g.lim <- c(min(y.g[,i],na.rm=TRUE),max(y.g[,i],na.rm=TRUE))
		# must be in appropriate order
		x <- as.data.frame(x.g[order(x.g[,i]),])
		y <- as.data.frame(y.g[order(x.g[,i]),])
#		lb <- as.data.frame(lb.se[order(x.g[,i]),])
#		ub <- as.data.frame(ub.se[order(x.g[,i]),])
		
		plot(x[,i],y[,i],type="l",ylim=y.g.lim,xlab=vars.smo[i],ylab=y.g.lab[i],...)
#		lines(x[,i],lb[,i],type="l",col="red",ylim=y.g.lim,xlab=vars.smo[i],ylab=y.g.lab[i],...)
#		lines(x[,i],ub[,i],type="l",col="red",ylim=y.g.lim,xlab=vars.smo[i],ylab=y.g.lab[i],...)
		rug(x.g[,i])
		}
	}
else
	cat("\nFull parametric model. Nothing left to do.\n")
#options(warn=warn.opt)
}


f <- function(factorvar)
# builds factor data matrix
{
factorname <- deparse(substitute(factorvar))
retval <- list(factor=factorname)
class(retval) <- "factor.info"
return(retval)
}


g <- function(var,df=NULL)
# extracts spline information from the formula term g(var,df,fx)
{
varname <- deparse(substitute(var))
retval <- list(var=varname,df=df)
class(retval) <- "spline.info"
return(retval)
}


fnz <- function(y)
# returns first non-zero observation index (base 1)
{
for (t in 1:length(y))
	if (y[t] > 0)
    	return(t)
}


framebuilder <- function(formula,dataset)
# builds a data frame from the dataset given the formula
{
if (!is.null(formula))
	frame <- model.frame(formula,dataset,na.action=na.pass)
else
	frame <- NULL
retval <- frame
class(retval) <- "data.frame"
return(retval)
}


elapsedtime <- function(st,et)
# Computes the elapsed time between st (start time) and et (end time)
{
time <- et[3]-st[3]	# gets time from the third position of the vectors
h <- trunc(time/3600)
if (h<10)
	hs <- paste("0",h,sep="",collapse=" ")
else
	hs <- h
time <- time-h*3600
min <- trunc(time/60)
if (min<10)
	mins <- paste("0",min,sep="",collapse=" ")
else
	mins <- min
time <- time-min*60
sec <- trunc(time)
if (sec<10)
	secs <- paste("0",sec,sep="",collapse=" ")
else
	secs <- sec

retval <- paste(hs,":",mins,":",secs,sep="",collapse=" ")
return(retval)
}


link <- function(x,link="log",inv=FALSE)
# apllies the link function
{
if (link=="log")
	{
	if (!inv)
		linked <- log(x)
	else
		linked <- exp(x)
	}
else
	stop(paste("Link funtion",link,"not implemented yet."))
return(linked)
}


lpnorm <- function(seq1,seq2=0,p=0)
# returns the Lp-norm. If p=0 then Infinity norm is returned
{
if (p == 0)
	norm <- max(abs(seq1-seq2),na.rm=TRUE)
if (p == 1)
	norm <- sum(abs(seq1-seq2),na.rm=TRUE)
if (p >= 2)
	norm <- sum((seq1-seq2)^p,na.rm=TRUE)
if (p > 2)
	warning("Lp-norm where p is greater than 2 is quite unusual.")

return(norm)
}


intensity <- function(w,x)
# spectral analysis of series x
{
n <- length(x)
t <- seq(1:n)
sp <- ((sum(x*cos(w*t)))^2+(sum(x*sin(w*t)))^2)/n
return(sp)
}


periodogram <- function(y,rows=trunc(length(na.omit(y))/2-1),plot=TRUE,...)
# creates and plots periodogram of series x
{
# initialization
if (rows > trunc(length(na.omit(y))/2-1))
	rows <- trunc(length(na.omit(y))/2-1)
x <- na.omit(y)
n <- length(x)
i <- seq(1:trunc(n/2-1))
omega <- (2*pi*i)/n

# Iomega <- sapply(i,function(i){intensity(omega[i], x=x)})
Iomega <- sapply(omega[i],intensity,x=x)
period <- (2*pi)/omega
period.max <- round(max(period),2)
period.min <- round(min(period),2)

if (plot)
	{
	# plots the periodogram
	plot(omega, Iomega, xlab="Angular frequency omega (rad)\n[Period on top axis]",ylab="I(omega)",...)
	axis(3,at=c(min(omega),0.5,1.0,1.5,2.0,2.5,3.0,max(omega)),        labels=c(period.max,12.57,6.28,4.19,3.14,2.51,2.09,period.min))
	title(main=paste("Raw Periodogram of Series",deparse(substitute(y))))
	lines(omega,Iomega,type="h")
	}

# returns periodogram
periodogram <- cbind.data.frame(period,omega,Iomega)
periodogram <- periodogram[order(periodogram$Iomega,decreasing=TRUE),]

retval <- periodogram[1:rows,]
return(retval)
}


envelope <- function(object,...)
# generic function for simulated envelope generation
	UseMethod("envelope")
	

envelope.pgam <- function(object,type="deviance",size=.95,rep=19,optim.method=NULL,epsilon=1e-3,maxit=1e2,plot=TRUE,title="Simulated Envelope of Residuals",verbose=FALSE,...)
# simulates and plots an envelope of residuals based on A. C. Atkinson book
{
if (class(object) == "pgam")
	{
	st <- proc.time()
	if (rep < 19)
        stop("Number of replications must equal or be greater than 19.")
	
	if (is.null(optim.method))
		optim.method <- object$optim.method
	n <- object$n
	tau <- object$tau
	dataset <- as.data.frame(eval(parse(text=object$dataset)))
	fitted <- fitted(object)[(tau+1):n]
	formula <- as.formula(paste("SIMRESP",object$formula[1],object$formula[3]))
	resid <- resid(object,type)[(tau+1):n]
	e <- matrix(NA,(n-tau),rep)
	
	attach(dataset)
	for (i in 1:rep)
        {
        if (verbose)
			cat(paste("\nReplication:",i,"\n"))
		else
			cat(paste("[",i,"]",sep=""))
		SIMRESP <- rpois(object$n.obs,fitted)
        runningmodel <- pgam(formula,dataset=cbind.data.frame(SIMRESP,dataset),omega=object$omega,beta=object$beta,offset=object$offset,maxit=1e2,eps=1e-4,optim.method=optim.method,se.estimation="none",verbose=verbose,lfn.scale=1e3)
        # for now these will be the defaults ---> must be changed!
		runningresid <- resid(runningmodel,type)[(tau+1):n]
		e[,i] <- sort(runningresid,na.last=TRUE,method="shell")
		}

	e1 <- numeric(n-tau)
	e2 <- numeric(n-tau)

	if (rep ==19)
		for (i in 1:(n-tau))
        	{
        	eo <- sort(e[i,])
        	e1[i] <- min(eo)
        	e2[i] <- max(eo)
        	}
	else
		for (i in 1:(n-tau))
        	{
        	eo <- sort(e[i,])
        	e1[i] <- eo[ceiling(((1-size)/2)*rep)]
        	e2[i] <- eo[ceiling(((1+size)/2)*rep)]
        	}
	residmean <- apply(e,1,mean,na.rm=TRUE)
	band <- range(resid,e1,e2,na.rm=TRUE)

	et <- proc.time()
	if (verbose)
		cat(paste("\nOverall envelope simulation process took", elapsedtime(st,et)),"(hh:mm:ss)\n")
	else
		cat("\n")

	retval <- list(lb=e1,ub=e2,mean=residmean,residuals=resid)
	class(retval) <- "envelope"
	}
else if (class(object) == "envelope")
	{
	e1 <- object$lb
	e2 <- object$ub
	residmean <- object$mean
	resid <- object$residuals
	band <- range(resid,e1,e2,na.rm=TRUE)
	}

# plotting the envelope
if (plot)
	{
	qqnorm(e1,axes=FALSE,main="",xlab="",ylab="",type="l",ylim=band,lty=1,lwd=1,col="red",bg="white")
	par(new=TRUE)
	qqnorm(e2,axes=FALSE,main="",xlab="",ylab="",type="l",ylim=band,lty=1,lwd=1,col="red",bg="transparent")
	par(new=TRUE)
	qqnorm(residmean,axes=FALSE,main="",xlab="",ylab="",type="l",ylim=band,lty=2,col="blue",bg="transparent")
	par(new=TRUE)
	qqnorm(resid,main=title, ylab="Residual Components",xlab="Standard Normal Quantiles",ylim=band,bg="transparent",...)
	}
if (class(object) == "pgam")
	return(retval)
}


tbl2tex <- function(tbl,label="tbl:label(must_be_changed!)",caption="Table generated with tbl2tex.",centered=TRUE,alignment="center",digits=getOption("digits"),hline=TRUE,vline=TRUE,file="",topleftcell="   ")
# outputs a data matrix in LaTeX format
{
tbl <- as.matrix(tbl)
r <- dim(tbl)[[1]]
c <- dim(tbl)[[2]]
rnames <- dimnames(tbl)[[1]]
cnames <- dimnames(tbl)[[2]]
align <- rep(substr(alignment,1,1),c+1)

# vertical lines
if (vline)
	pos <- paste("|",paste(align,sep="",collapse="|"),"|",sep="")
else
	pos <- paste(align,collapse=" ")

# headers
cat("% File generated in R environment by tbl2tex()\n% Creation date: ",date(),"\n% Suggestions to Washington Junger <wjunger@ims.uerj.br>\n",sep="",file=file,append=FALSE)
cat("\\begin{table}\n",ifelse(centered,"\\centering\n",""),"\\caption{",caption,"}\n\\label{",label,"}\n",sep="",file=file,append=TRUE)
cat("    \\begin{tabular}{",pos,"} ",ifelse(hline,"\\hline",""),"\n",sep="",file=file,append=TRUE)
cat("    ",topleftcell," & ",paste(cnames,sep="",collapse=" & ")," \\\\ ",ifelse(hline,"\\hline",""),"\n",sep="",file=file,append=TRUE)
# data
for(i in 1:r)
	cat("    ",rnames[i]," & ",paste(round(tbl[i,],digits),sep="",collapse=" & ")," \\\\ ",ifelse(hline,"\\hline",""),"\n",sep="",file=file,append=TRUE)
# footer
cat("    \\end{tabular}\n\\end{table}\n",sep="",file=file,append=TRUE)
cat("% End of tbl2tex() generated file.\n",sep="",file=file,append=TRUE)
}


# functions end here -----------------------

