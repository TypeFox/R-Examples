##########################################################################
## SAVE main Method
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 15, November 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################

`SAVE` <-
function (response.name=NULL, controllable.names=NULL,
			calibration.names=NULL,field.data=NULL,
            model.data=NULL,mean.formula=~1,bestguess=NULL,
            kriging.controls = SAVE.controls(),verbose=FALSE){

	model <- new ("SAVE")
	dprct <-.deprecate.parameters(call=sys.call(sys.parent(1)))
  	model@call <- as.call(dprct)
	
	#establish the wd:
	wd <- tempdir()
	model@wd <- paste(wd,'/',sep='')	
				

	if (is.null(model.data) || length(model.data)==0)
		stop ("Model data cannot be NULL\n")
	
	if (is.null(response.name) || length(response.name)==0)
		stop ("Response name cannot be NULL\n")
		
	if (length(response.name)>1)
		stop ("Response has to be univariate\n")

	if (sum(response.name==controllable.names)>0 || sum(response.name==calibration.names)>0)
		stop ("Response variable cannot be a calibration nor a controllable input\n")
				
	if ((length(unique(controllable.names)) != length(controllable.names)) || 
		(length(unique(calibration.names)) != length(calibration.names)))
		stop ("Calibration or Controllable parameters contain duplicated inputs\n")
	
	model@responsename <- as.character(response.name)
	   
	#####
	#Field data:
  	if (is.null(field.data) || length(field.data)==0){
		# By now we are not allowing the field.data to be null
		# in the future we will remove this stop
		stop("Field data cannot be NULL\n")
  	  	if (!is.null(bestguess) & length(bestguess)!=0)
  			cat("The parameter bestguess unused and set to NULL since field data is not provided\n")
  		cat("The result wont have Stage II parameters.\n")
		bestguess <- as.numeric(NULL)
		# To change the call and remove bestguess from the call since field data is not provided 
		aux <- which(names(model@call)=="bestguess")
		model@call <- model@call[-aux]
	}else {
		#Select those columns corresponding to the response variable:
        	yf <- as.vector(t(field.data[,model@responsename]))
        	model@yf <- yf
	}	

    if (!is.null(controllable.names) && 
		(length(unique(field.data[,controllable.names]))==1)) 
			stop("The current field data indicates that the controllable parameters are fixed at a single value in the field experiment.\n 
					 In this case, controllable.names parameter should be set to NULL.\n")

	if (is.null(calibration.names) || length(calibration.names)==0){
		stop ("Calibration parameters cannot be NULL in the current implementation of the package.")
		# We are currently working on allowing the package to work with no Calibration parameters
		# the alternative using bestguess is not fully tested yet, so that is why we are are stopping here.
		# April, 29th, 2014
		if (!is.null(bestguess) & length(bestguess)!=0)
			cat("The parameter bestguess unused and set to NULL since calibration parameters have not been provided\n")
		bestguess <- as.numeric(NULL)
		# To change the call and remove bestguess from the call since calibration parameters have not been provided 
		aux <- which(names(model@call)=="bestguess")
		model@call <- model@call[-aux]
	}else {
		if (is.null(bestguess) || length(bestguess)==0)
			stop("bestguess cannot be NULL when there are calibration parameters.")
	}	
	model@calibrationnames <- as.character(calibration.names)

    model@constant.controllables <- F # The Flag for constan controllable inputs is set to its default value
                                      # i.e. the controllables are not constant by default.
	if (is.null(controllable.names)){model@constant.controllables <- T}
    ####
	if (!model@constant.controllables){
		df <- as.matrix(field.data[,controllable.names])		
		colnames(df) <- controllable.names
		model@df <- df
		model@controllablenames <- as.character(controllable.names)
	}

	model@meanformula <- mean.formula
				
	####
    # Estimate stage I parameters using DiceKriging

    #We differentiate between the two cases: varying controllable inputs and constant controllable inputs
    #in the second case model.data to be given to the optimizer must only contain calibration inputs.
    if (!model@constant.controllables){ #Case 1
		#Model data:
        dm <- as.matrix(model.data[,c(model@controllablenames,model@calibrationnames)])
        colnames(dm) <- c(model@controllablenames,model@calibrationnames)
	    x.unique <- as.data.frame(unique(model@df))
	    names(x.unique) <- model@controllablenames
	    model@xf <- model.matrix(model@meanformula,x.unique)
    } else{ #Case 2: Constant controllable inputs
		#Model data:
		model@xf <- matrix(1, ncol=1, nrow=1)
        dm <- as.matrix(model.data[,model@calibrationnames])
        colnames(dm) <- model@calibrationnames
    }

    #####
    #####
    #Model data:
    model@dm <- dm
    #write.table(as.data.frame(model@dm), file=paste(model@wd,"/model_inputs.dat",sep=""),
    #            col.names=F, row.names=F)

    ym <- as.vector(t(model.data[,model@responsename]))
    model@ym <- ym
    #write.table(as.vector(model@ym), file=paste(model@wd,"/model_data.dat",sep=""),
    #            col.names=F, row.names=F)

    #Detect if mean.formula contains terms that involve calibration inputs and remove them with a warning
	#Also, keep the explicit expression for formula
	mean.formula<- drop.response(as.formula(mean.formula), data=model.data[,c(controllable.names,calibration.names)])
	if (!is.null(calibration.names) & length(calibration.names)>0){
		termstoremove<- integer(0)
		for (i in 1:length(calibration.names)){
			termstoremove<- c(termstoremove, grep(calibration.names[i],attr(terms(mean.formula),"term.labels")))
		}
		if(length(termstoremove)>0){
			warning("Terms involving calibration inputs have been suppressed in the mean formula\n")
			mean.formula<- drop.terms(termobj=terms(mean.formula), dropx = termstoremove)
			mean.formula<- reformulate(attr(terms(mean.formula), "term.labels"))
		}
	}

    #write.table(model@xm,file=paste(model@wd,"/mcmc.field.design.M.matrix.dat",sep=""),
    #            col.names=F, row.names=F)
    #write.table(model@xf,file=paste(model@wd,"/mcmc.field.design.F.matrix.dat",sep=""),
    #            col.names=F, row.names=F)

    ####
    # Estimate stage I parameters using DiceKriging
    #m <- km(model@meanformula,design=model@dm,response=model@ym,covtype="powexp")
    # Now we have incorporated the possibility of controlling the parameters
    # lower, upper, optim.method and parinit on km()
    m <- do.call("km",c(list("formula"=model@meanformula,"design"=model@dm, "response"=model@ym,"covtype"="powexp"),kriging.controls))
    alphaM <- (m@covariance)@shape.val
    betaM <- ((m@covariance)@range.val)^(-alphaM)
    lambdaM <- 1/((m@covariance)@sd2)
    thtaM <- c(lambdaM,betaM,alphaM)
    
    #We distinguish between the cases: varying and constant controllable inputs:
    if (!model@constant.controllables){ #Case 1
        names (thtaM) <- c("lambdaM",paste("betaM",c(model@controllablenames,model@calibrationnames),sep='.'),paste("alphaM",c(model@controllablenames,model@calibrationnames),sep='.'))
    }
    else{ #Case 2: Constant controllable inputs
        names (thtaM) <- c("lambdaM",paste("betaM",model@calibrationnames,sep='.'),paste("alphaM",model@calibrationnames,sep='.'))	
	}
    thtaL <- m@trend.coef

    ####
    # Estimate stage II parameters

    # predict code at field design augmented with guess of u

    # produce design
		#takes the values on the list u.guess and puts on the right order (the one given by calibration.names)
	if (!is.null(bestguess) & length(bestguess)!=0){
			model@bestguess<- vapply(bestguess, FUN=function(x){x[1]}, FUN.VALUE=c(0))[model@calibrationnames]
		}
		else model@bestguess<- as.numeric(NULL)
				
		if (is.null(model@calibrationnames) || length(model@calibrationnames)==0){
			xnew<- model@df
		}
		else{
		    aux <- matrix(rep(model@bestguess,length(model@yf)),ncol=length(model@bestguess),byrow=T)
 			aux <- as.data.frame(aux)
			names(aux) <- model@calibrationnames
                    #Here again we distinguish between varying or constant controllable inputs
                    if (!model@constant.controllables){
                        xnew <- data.frame(model@df,aux)
                        names(xnew) <- c(model@controllablenames,model@calibrationnames)
                    }
                    else{ # Constant Controllable inputs
                        xnew <- data.frame(aux)
                        names(xnew) <- model@calibrationnames
                    }
        }
    # predict code
    ymnew <- predict(m,newdata=xnew,type="UK",se.compute=FALSE,
                         light.return=TRUE)

    # resulting bias
    biashat <- model@yf-ymnew$mean
    # print(df)

	####################	
	#load some auxiliary functions
	duplicates <- function(dat)
	{
		s <- do.call("order", as.data.frame(dat))
		if(dim(as.matrix(dat))[2]==1){
			non.dup <- !duplicated(as.matrix(dat[s]))
		}
		else{
			non.dup <- !duplicated(as.matrix(dat[s,]))
		}
		orig.ind <- s[non.dup]
		first.occ <- orig.ind[cumsum(non.dup)]
		#first.occ[non.dup] <- NA
		first.occ[order(s)]
	}
	###############end of loading auxiliary funcs.	

	if(model@constant.controllables){ # if we only have one distinct input value we use optimize.c
		#####
		#Values of the parameters:
		sum2 <- sum(biashat^2) # sum of squares
		tot2 <- sum(biashat)^2 # total square
		n <- length(biashat) # sample size
		write(c(n,sum2,tot2),file=paste(model@wd,"biasaux.tmp",sep="/"),ncolumns=3)
		maxiterations = 1000;
		eps = 1.E-7;
		err = 1.;
		psi0 = 0.5;
		
		#call optimize.c
		output <- .C('optimize',as.integer(verbose),
			as.integer(maxiterations),as.double(eps),
			as.double(err),as.double(psi0),
			as.character(model@wd))
	
		maux <- scan(file=paste(model@wd,"thetaF_mle.dat",sep="/"))
		lambdab <- 1/(maux[1])
		alphab <- maux[2]
		betab <- maux[3]
		lambdaF <- 1/(maux[4])
		thtaF <- c(lambdab,lambdaF)
		names (thtaF) <- c("lambdaB","lambdaF")		
	} else { # we call again DiceKriging
		# estimate parameters
    	# trend <- 0.0
    	maux <- km(formula=~1,design=model@df,response=biashat,covtype="gauss",
                   # coef.trend=trend,
                   nugget.estim=TRUE)
		betab <- ((maux@covariance)@range.val)^(-2)/2
		lambdab <- 1/((maux@covariance)@sd2)
		lambdaF <- 1/((maux@covariance)@nugget)
		alphab <- rep(2.0,length(betab))
		thtaF <- c(lambdab,betab,alphab,lambdaF)
		names (thtaF) <- c("lambdaB",paste("betaB",model@controllablenames,sep='.'),paste("alphab",model@controllablenames,sep='.'),"lambdaF")	    
		
	   }
    	
	model@mle <-list(thetaL=thtaL,thetaM=thtaM,thetaF=thtaF)
				
	
	## End of gaspfit
    #####
	#####
	#Model data:
	x.m<- as.data.frame(model.data[,c(model@controllablenames,model@calibrationnames)])
	names(x.m) <- c(model@controllablenames,model@calibrationnames)
	#####
	#####
	#Mean response:
    model@xm <- model.matrix(model@meanformula, x.m)

    if (!model@constant.controllables){dthetaM <- length(model@calibrationnames) + length(model@controllablenames)}
    else{dthetaM <- length(model@calibrationnames)}
    
    dthetaM <- dthetaM*2+1
    dthetaL <- ncol(model@xm)
    
    if (!model@constant.controllables){dthetaF <- 2*length(model@controllablenames)+2}
    else{dthetaF <- 2}
    
    if(length(model@mle$thetaM)!=dthetaM)
    	{#model@mle<- NULL;
        stop("Dimension of MLE for the covariance of the computer model is incorrect\n")}
    if(length(model@mle$thetaL)!=dthetaL)
    	{#model@mle<- NULL;
        stop("Dimension of MLE for the mean of the computer model is incorrect\n")}
    if(length(model@mle$thetaF)!=dthetaF)
    	{#model@mle<- NULL;
        stop("Dimension of MLE for the second stage parameters is incorrect\n")}


#validObject(model, complete=TRUE)
				
	unlink(paste0(model@wd,'/*'))

	return(model)	
}

####################
#load some auxiliary functions
normal<- function(var.name, mean, sd, lower, upper){
	if ((upper >= lower) && (mean >=lower) && (mean <= upper) && (sd >= 0))
		c(var.name, 1, lower, upper, mean, sd^2)
	else stop ("Error: the specified parameters are not appropriate for a normal.\n")
}

uniform<- function(var.name, lower, upper){
	if (lower <= upper)
			c(var.name, 0, lower, upper, 0, 0)
	else stop ("Error: the specified lower parameter is greater than the upper.\n")
}

.expand.call <- function(call=sys.call(sys.parent(2)), expand.dots = TRUE)
{
	ans <- as.list(call)
	frmls <- formals(deparse(ans[[1]]))
	# To remove '...'
	if (names(frmls[length(frmls)])=='...') frmls <- frmls[-length(frmls)]
	add <- which(!(names(frmls) %in% names(ans)))
	return(as.call(c(ans, frmls[add])))
}

.deprecate.parameters <- function(call=sys.call(sys.parent(1)), expand.dots = TRUE)
{
	ans <- as.list(.expand.call(call=call))
	frmls <- formals(deparse(ans[[1]]))
	aans <- ans [-1]
	if (names(frmls[length(frmls)])=='...') frmls <- frmls[-length(frmls)]
	#print(names(ans))
	#print(names(frmls))
        deprct <- which(!(names(aans) %in% names(frmls)))
        if (length(deprct) != 0) {
            if (as.character(names(aans)[deprct])!=""){ 
				cat('Warning!: Deprecated parameters:\n')
				print(as.character(names(aans)[deprct]))
			}
	    ans <- ans[-(deprct+1)]
        }
        #return (deprct)
        #return(call)
        #print(ans)
        return (ans)
}
