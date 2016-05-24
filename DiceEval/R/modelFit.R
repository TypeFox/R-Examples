modelFit <- function (X,Y, type, ...)
{
  data <- data.frame(X,Y)
  X    <- data.frame(X)
  f    <- dim(X)[2]
  
  if (any(sapply(data,"is.factor")))
	stop("columns of class \"factor\" in arg 'data' is not allowed")

  argList<-list(...)
  if (type == "Linear"){
	if(is.null(argList$formula)==TRUE){
	  warning("[type==\"Linear\"] argument \'formula\' not found, set at \'Y~.\'")
	  fmla <- as.formula(paste(names(data[f+1]),"~.",sep=""))
	} else 	fmla <- argList$formula
	model <- lm(fmla,data)
	out   <- list(data=list(X=X,Y=data[,f+1,drop = TRUE]),type=type,formula=fmla,model=model)
 } else if(type == "Additive"){
	if (is.null(argList$formula)==TRUE){
	  warning("[type==\"Additive\"] argument \'formula\' not found, set at \'Y~~ s(X1)+...+s(Xp)\'")
	  fmla <- formulaAm(X,Y) 
	} else 	fmla <- argList$formula
	model <- gam::gam(formula=fmla, data=data)
	out <- list(data=list(X=X,Y=data[,f+1,drop = TRUE]),type=type,formula=fmla,model=model)
  } else if (type == "MARS"){
	if(is.null(argList$degree)==TRUE){
	  warning("[type==\"MARS\"] argument \'degree\' not found, set at \'2\'")
	  degree <- 2
	} else 	degree <- argList$degree;
	model <- mda::mars(data[,1:f],data[,f+1],degree=degree)
	out <- list(data=list(X=X,Y=data[,f+1,drop = TRUE]),type=type,degree=degree,model=model)
  } else if (type == "PolyMARS"){
	if(is.null(argList$gcv)==TRUE){
	  warning("[type==\"PolyMARS\"] argument \'gcv\' not found, set at \'4\'")
	  a <- 4
	} else 	a <- argList$gcv;	
	model <- polspline::polymars(data[,f+1],data[,1:f],gcv=a)
	out <- list(data=list(X=X,Y=data[,f+1,drop = TRUE]),type=type,gcv=a,model=model)
  } else if (type == "StepLinear"){
	if(is.null(argList$formula)==TRUE){
	  warning("[type==\"StepLinear\"] argument \'formula\' not found, set at \'Y~.\'")
	  fmla <- as.formula(paste(names(data[f+1]),"~.",sep=""))	  
	} else 	fmla <- argList$formula
	if(is.null(argList$penalty)==TRUE){
	  warning("[type==\"StepLinear\"] argument \'penalty\' not found, set at \'2\' (AIC criteria)")
	  p <- 2
	} else p <- argList$penalty
	init <- lm(fmla,data=data)
	if (length(init$coefficients)>dim(data)[1]) {
		stop("[type==\"StepLinear\"] There are too many terms into the full length model")
	}
	model <- step(init,trace=FALSE,k=p)
	out <- list(data=list(X=X,Y=data[,f+1,drop = TRUE]), type=type, formula=fmla, penalty=p, model=model)
  } else if (type == "Kriging"){
	model <- km(design=data[,1:f],response=data[,f+1],...)
	out <- list(data=list(X=X,Y=data[,f+1,drop = TRUE]),type=type,model=model)
  } else stop("The argument 'type' must be 'Linear', 'Additive', 'MARS', 'PolyMARS', 'StepLinear' or 'Kriging'")
	return(out)
  } 