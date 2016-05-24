predict.manyglm=function(object,newdata = NULL,type = c("link", "response", 
    "terms"), se.fit = FALSE, dispersion=NULL, terms=NULL, 
    na.action=na.pass, ...)
{

        if (object$family=="gaussian") {
	    if (type=="link") {
	        stop("Possible type of predict.manylm is 'response' or 'term'.")
	    }
	    else {
	    return(predict.manylm(object, newdata=newdata, se.fit=se.fit, type=type, 
		              terms=terms, na.action = na.pass, ... ))
	    }	      
        }
	# set up matrices
	nVar <- NCOL(object$fitted.values)
	nObs <- NROW(object$fitted.values)
        if(is.null(newdata)==F)
            nObs = NROW(newdata)
 	ses <- fts <- matrix(NA,nObs,nVar)
	if(is.null(newdata))
	   dimnames(fts)[[1]] = dimnames(object$fitted.values)[[1]]
        else
           dimnames(fts)[[1]] = rownames(newdata)
  dimnames(fts)[[2]] = dimnames(object$fitted.values)[[2]]
        
	type <- match.arg(type)
	na.act <- object$na.action
	object$na.action <- NULL
 	
	fm <- formula(object)

	if( "K" %in% names(object$call) ) #for predictions when K is non-zero:
	{
	  K = eval(object$call$K)
	  if(K>1 & object$family!="binomial")
	  {
	    warning("Argument K should only be specified when family is binomial, K reset to 1")
	    K=1
	  }      
	}
	else
	  K = 1
    
	# use predict.glm to compute each column one at a time
	for(iVar in 1:nVar)
	{
	       fam = switch(object$family,
	             "binomial(link=logit)"=binomial(),
	             "binomial(link=cloglog)"=binomial("cloglog"),
	             "poisson"=poisson(),
	             "gaussian"=gaussian(),
	             "negative.binomial"=negative.binomial(theta=object$theta[iVar])
        	)
	        if(is.null(object$data))
	            dat.i = model.frame(object)
	        else
              dat.i = data.frame(object$y[,iVar], object$data)
          if(K>1)
          {
            dat.i$K=K
            form <- as.formula(paste("object$y[ ,", iVar, "]/K ~ ", fm[3]))
            object.i = glm(form, family=fam, data=dat.i, weights=K, start=as.vector(object$coef[,iVar])) 
          }
          else
          {
            form <- as.formula(paste("object$y[ ,", iVar, "] ~ ", fm[3]))
            object.i = glm(form, family=fam, data=dat.i, start=as.vector(object$coef[,iVar]))
            #DW, 18/11/14:  starting value included in predict to avoid errors with sparse data
          }
	       
          ft.i <- predict.glm(object.i, newdata=newdata, se.fit=se.fit, 
                                 type=type, terms = terms, na.action = na.action)
          if(se.fit==T)
		      {
	            fts[,iVar] = ft.i$fit
         	    ses[,iVar] = ft.i$se
      		}
		      else
		        fts[,iVar] = ft.i
  }
	if(se.fit)
		out = list(fit=fts,se.fit=ses)
	else
		out=fts
	return(out)
}
