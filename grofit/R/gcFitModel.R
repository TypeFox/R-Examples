gcFitModel <-
function(time, data, gcID ="undefined", control=grofit.control())
{

# /// check input parameters
if (is(control)!="grofit.control") stop("control must be of class grofit.control!")
if (!control$fit.opt%in%c("m","b")) stop("Fit option is not set for a model fit. See grofit.control()")

# /// conversion to handle even data.frame inputs
time <- as.vector(as.numeric(as.matrix(time)))
data    <- as.vector(as.numeric(as.matrix(data)))

# /// check length of input data
if (length(time)!=length(data)) stop("gcFitModel: length of input vectors differ!")

# /// determine which values are not valid
bad.values <- (is.na(time))|(time<0)|(is.na(data))|(data<0)|(!is.numeric(time))|(!is.numeric(data))

# /// remove bad values or stop program
if (TRUE%in%bad.values)
{
   if (control$neg.nan.act==FALSE)
   {
      time    <- time[!bad.values]
      data    <- data[!bad.values]
   }
   else{
   stop("Bad values in gcFitModel")
   }
}

# /// check if there are enough data points
if (length(data)<5){
	warning("gcFitModel: There is not enough valid data. Must have at least 5 unique values!")
        gcFitModel   <- list(raw.time = time, raw.data = data, gcID = gcID, fit.time = NA, fit.data = NA, parameters = list(A=NA, mu=NA, lambda=NA, integral=NA), model = NA, nls = NA, reliable=NULL, fitFlag=FALSE, control = control)
	class(gcFitModel) <- "gcFitModel"
        return(gcFitModel)
}
else{
	# /// apply transformation
	if (control$log.x.gc==TRUE){ time <- log(1+time)}
	if (control$log.y.gc==TRUE){ data <- log(1+data)}
	
	# fitting parametric growth curves
	y.model     <- NULL
	bestAIC     <- NULL
	best        <- NULL
	used        <- NULL
	
	# starting values for parametric fitting from spline fit
	nonpara     <- gcFitSpline(time, data, control)
	mu.low      <- nonpara$parametersLowess$mu
	lambda.low  <- nonpara$parametersLowess$lambda
	A.low       <- nonpara$parametersLowess$A

	# /// determine length of model names
	l               <- 10
	for (i in 1:length(control$model.type)){
		l[i] <- nchar((control$model.type)[i])
	}
	lmax <- max(l)
	
	# /// loop over all parametric models requested
	for (i in 1:length(control$model.type)){
		if (control$suppress.messages==FALSE){
			cat(paste("--> Try to fit model", (control$model.type)[i]))
		}
		initmodel    <- paste("init",(control$model.type)[i],sep="")
		formulamodel <- as.formula(paste("data ~ ", (control$model.type)[i], "(time, A, mu, lambda, addpar)", sep=""))
		if ( (exists((control$model.type)[i])) && (exists(initmodel))){
			init.model  <- do.call(initmodel, list(y=data, time=time, A=A.low, mu=mu.low, lambda=lambda.low))
			try(y.model <- nls(formulamodel, start=init.model), silent = TRUE)
			if(!(TRUE%in%is.null(y.model))){
				AIC       <- AIC(y.model)
			}

		if (control$suppress.messages==FALSE){
			if (class(y.model)=="nls"){
				if (y.model$convInfo$isConv==TRUE){
					message(paste(rep(".",lmax+3-l[i])), " OK")	
				}
				else{
					warning(paste(rep(".",lmax+3-l[i]), " nls() failed to converge with stopCode ", as.character(y.model$convInfo$stopCode)))
				}
			}
			else{
                                        message(paste(rep(".",lmax+3-l[i]))," ERROR in nls(). For further information see help(gcFitModel)")
			}
		}
			if (FALSE%in%is.null(AIC)){
				if (is.null(best) || AIC<bestAIC){
				bestAIC <- AIC
				best    <- y.model
				used    <- (control$model.type)[i]
			}
			}
		} # of if ( (exists((control$model.type)[i])) ...
		else{
		cat((control$model.type)[i])
		cat("\n")
		cat(initmodel)
		cat("\n")
		stop("The model definition above does not exist! Spelling error?")
		}
y.model <- NULL
	}
	
	if (control$suppress.messages==FALSE) cat("\n") 
	
	# /// extract parameters from data fit
	if(is.null(best)==FALSE){
		Abest      <- summary(best)$parameters["A",1:2]
		mubest     <- summary(best)$parameters["mu",1:2]
		lambdabest <- summary(best)$parameters["lambda",1:2]
		fitFlag    <- TRUE
		if (length(time)==length(as.numeric(fitted.values(best)))){
			Integralbest <- low.integrate(time, as.numeric(fitted.values(best)) )
		}
		else{
			Integralbest <- NA
		}
	}
	else{
		warning("gcFitModel: Unable to fit this curve parametrically!")
		Abest        <- c(NA,NA)
		mubest       <- c(NA,NA)
		lambdabest   <- c(NA,NA)
		Integralbest <- NA
		fitFlag      <- FALSE
	}
} # end of else if(length(y)<5)

gcFitModel <- list(raw.time = time, raw.data = data, gcID = gcID, fit.time = time, fit.data = as.numeric(fitted.values(best)), parameters = list(A=Abest, mu=mubest, lambda=lambdabest, integral=Integralbest), model = used, nls = best, reliable=NULL, fitFlag=fitFlag, control = control)

class(gcFitModel) <- "gcFitModel"

gcFitModel
}

