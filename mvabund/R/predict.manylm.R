predict.manylm=function(object,newdata = NULL, se.fit=FALSE, 
type = c("response", "terms"), terms=NULL, na.action=na.pass, ...)
{
	# set up matrices
	nVar <- NCOL(object$fitted.values)
	nObs <- NROW(object$fitted.values)
        if(is.null(newdata)==F)
            nObs = NROW(newdata)
 	ses <- fts <- matrix(NA,nObs,nVar)
  dimnames(fts)=vector(2,mode="list") #to avoid error making a null thing a list.
	dimnames(fts)[[2]] = dimnames(object$fitted.values)[[2]]
	if(is.null(newdata))
	   dimnames(fts)[[1]] = dimnames(object$fitted.values)[[1]]
        else
           dimnames(fts)[[1]] = rownames(newdata)
	
	type <- match.arg(type)
	na.act <- object$na.action
	object$na.action <- NULL
 	
	fm <- formula(object)
        
	# use predict.glm to compute each column one at a time
	for(iVar in 1:nVar)
	{
                form <- as.formula(paste("object$y[ ,", iVar, "] ~ ", fm[3]))
		if (is.null(object$data))
		    dat.i = model.frame(object)
                else 
		    dat.i  = data.frame(object$y[, iVar], object$data) 
                object.i = lm(form, data=dat.i)
                ft.i <- predict.lm(object.i, newdata=newdata, se.fit=se.fit, 
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
