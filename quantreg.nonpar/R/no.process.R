no.process <-
function(data=data, taus, formula, basis=NULL, var, load, rearrange=F, rearrange.vars="quantile", 
                       average=T, nderivs=1, method="fn")
{

  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)] 
  mf$drop.unused.levels <- TRUE 
  mf[[1]] <- as.name("model.frame") 
  mf <- eval.parent(mf) 
  mt<-attr(mf,"terms")
  y<-model.response(mf,"numeric")
  x<-model.matrix(mt,mf)
  if(all(x[,1]==rep(1,dim(x)[1]))){x <- x[,-1]} # Remove intercept from x matrix
  
  
  z<-as.vector(data[,var]) # Save vector with variable of interest
  f<-formula(y~x) # Save variables locally and create formula
  length.taus=length(taus)
  basis.type<-"noBasis"
  
  if (is.basis(basis)){
    if (basis$type=="bspline" | basis$type=="fourier"){
      basis.type <- "fda"
    } else {
      stop("basis must be a bspline basis, fourier basis, polynomial basis, or factor variable")
    }
  } else if (is.factor(basis)) {
    basis.type <- "factor"
  } else if (attr(basis,"class")[1]=="poly") {
    basis.type <- "polynomial"
  } else {
    stop("basis must be a bspline basis, fourier basis, polynomial basis, or factor variable")
  }
  
  if (rearrange==T & (nderivs!=0 | average!=F)) {
    stop( "Rearrangement may only be used when nderivs=0 and average=F") }
  
  if (rearrange==T & rearrange.vars!="quantile" & rearrange.vars!="var" & rearrange.vars!="both") {
    stop( "rearrange.vars must be quantile, var, or both")}
  
  if (average==T & nderivs==0) {
    stop( "This method is not available for average=T and nderivs=0")}
  
  if (is.null(basis)) {
    stop( "Please input a basis")
  }
  if(basis.type=="factor" & (average!=F | nderivs!=0)){
    stop( "Fully saturated indicator series approximation is only available when average=F and nderivs=0")
  }
  if((basis.type=="polynomial") & nderivs>2) {
    stop("Estimation using polynomials is not available for nderivs>2")
  }
  
  if (basis.type=="fda") {
    basvecs<-eval.basis(z,basis)[,-1] # Create vectors of nonparametric bases evaluated at observed variable values
  } else {
    basvecs<-basis
  }
  
  f<-update.formula(f, ~ . + basvecs) # Include series regressors in formula; supress intercept (included in nonparametric basis)
  if (is.null(n)) { # Set number of bootstrap repetitions to number of observations, if not user provided
    n <- dim(x)[1]
  }
  nObs <- dim(x)[1]
  
  # Generate loading vector if not input by user
  if (is.null(load)) {
    
    modelMatrix <- model.matrix(f, data)
    nVars <- dim(modelMatrix)[2]-dim(basvecs)[2]-1 # number of model variables, excluding the intercept and the var of interest
    
    
    if (average==T) { # generate loading vector
      
      # Get derivative of formula, not including the nonparametric components
      parDeriv <- formulaDeriv(formula,var,data=data,nDerivs=nderivs)   
      if (basis.type=="fda") {
        load <- cbind(matrix(apply(parDeriv, 2, mean), nrow=1), matrix(apply(as.matrix(eval.basis(z,basis,Lfdobj=nderivs)[,-1]), 2, mean), nrow=1)) 
      } else {
        load <- cbind(matrix(apply(parDeriv, 2, mean), nrow=1), matrix(apply(as.matrix(poly.wrap(x=z,degree=max(attr(basis,"degree")),nderivs=nderivs,coefs=attributes(basis)$coefs)), 2, mean), nrow=1)) 
      }
    } else {
      
      if (nderivs==0) {
        load <- as.matrix(aggregate(modelMatrix, by=list(z), FUN=load.sum)[,-1])
      } else {
        # Get derivative of formula, not including the nonparametric components
        parDeriv <- formulaDeriv(formula,var,data=data,nDerivs=nderivs)
        if (basis.type=="fda"){
          load <- as.matrix(aggregate(cbind(parDeriv, eval.basis(z,basis,Lfdobj=nderivs)[,-1]), by=list(z), FUN=mean)[,-1])
        } else {
          load <- as.matrix(aggregate(cbind(parDeriv, poly.wrap(x=z,degree=max(attr(basis,"degree")),nderivs=nderivs,coefs=attributes(basis)$coefs)), by=list(z), FUN=mean)[,-1])
        }
      }  
    }
    
  }
  
  length.taus=length(taus)
  loadLength=dim(load)[1]
  z.unique <- as.matrix(aggregate(z, by=list(z), FUN=mean)[,-1]) # generate vector with unique values of x
  no.proc  <- list(qfits=NULL, point.est=NULL, var.unique=NULL, load=NULL) # set up output vector
  no.proc$qfits <- vector("list",length.taus) # quantile regression fits
  
  err<-F
  for (i in 1:length.taus) {
    no.proc$qfits[[i]]<-tryCatch(expr=rq(f,tau=taus[i],data,method=method),error=function(e) {
      err<-T
      fit<-rq(f,tau=taus[i],data)
      return(fit)
    })
  }
  if (err==T) {
    message("Warning: rq method chosen produced an error; using method=br when necessary")
  }
  
  # Generate output, including t statistics
  fittedValues <- matrix(0, nrow=loadLength, ncol=length.taus) # predicted values from quantile regressions
  
  for (i in 1:length.taus){ # predict fitted values, evaluated at loading vector
    fittedValues[,i] <- load %*% no.proc$qfits[[i]]$coef
  }
  if (rearrange==T) { # perform rearrangement if requested
    if (length(taus)==1){
      fittedValues <- matrix(rearrangement(list(z.unique),as.vector(fittedValues)),nrow=loadLength,ncol=1)
    } else {
      fittedValues <- rearrangement(list(z.unique,taus),fittedValues)
    }
  }
  
  no.proc$var.unique <- z.unique
  
  no.proc$point.est <- matrix(NA, nrow=loadLength, ncol=length.taus) 
  no.proc$point.est <- fittedValues
  no.proc$load <- load
  
  return(no.proc);
}
