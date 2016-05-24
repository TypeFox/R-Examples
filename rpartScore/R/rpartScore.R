rpartScore <-
function(formula, data, weights, subset, na.action = na.rpart, 
    split="abs",prune="mc", model = FALSE, x = FALSE, y = TRUE, control,...) {
	call <- match.call()
    	if (is.data.frame(model)) {
        	m <- model
        	model <- FALSE
    	}
    	else {
        	m <- match.call(expand.dots = FALSE)
        	m$model <- m$split<-m$prune <- m$control <- NULL
        	m$x <- m$y <- m$parms <- m$... <- m$xval<-NULL
        	m$cost <- NULL
        	m$na.action <- na.action
        	m[[1L]] <- as.name("model.frame")
	
		m <- eval(m, parent.frame())
		wt<-model.extract(m, "weights")
		sum.wt <- ifelse(length(wt) == 0L,nrow(m),sum(model.extract(m, "weights")))
		}
	extraArgs <- list(...)
    	if (length(extraArgs)) {
        	controlargs <- names(formals(rpart.control))
        	indx <- match(names(extraArgs), controlargs, nomatch = 0)
        	if (any(indx == 0)) 
            	stop("Argument ", names(extraArgs)[indx == 0], "not matched")
    	}
    	controls <- rpart.control(...)
    	if (!missing(control)) 
        	controls[names(control)] <- control

	split.int <- pmatch(split, c("abs", "quad"))
     	if (is.na(split.int)) 
            stop("Invalid splitting criterion")

	prune.int <- pmatch(prune, c("mc", "mr"))
     	if (is.na(prune.int)) 
            stop("Invalid pruning criterion")

	
	splitnew<-function(y, wt, x, parms, continuous, minb) splitAbs(y, wt, x, parms, continuous,minb=controls$minbucket)
	methodnew<-list(eval=evalMedian,split=splitnew,init=initScore)
	
	if (prune=="mr") methodnew$eval<-evalMode
	if (split=="quad") {
		splitQuad.n<-function(y, wt, x, parms, continuous, sumwt) splitQuad(y, wt, x, parms, continuous,sumwt=sum.wt)
		methodnew$split<-splitQuad.n
	}

	c<-rpart(model=m,method=methodnew,control=controls)
	c[[3]]<-match.call()
	c[[length(c)+1]]<-m
	names(c)[length(c)]<-"model"
	
	cross<-c$control$xval
	
	if (((length(cross)==1) && (cross>1))|((length(cross)>1)& (var(cross)!=0))){

	Y<-model.extract(m, "response")
	wt <- model.extract(m, "weights")
	if (length(wt) == 0L) wt <- rep(1, nrow(m))

	nobs<-length(Y)
	if (length(cross) == 1L) {
        	xgroups <- sample(rep(1:cross, length = nobs), nobs, replace = FALSE)
    	}
    	else if (length(cross) == nobs) {
        	xgroups <- cross
        	cross <- length(unique(xgroups))
    	}
    	else {
        	if (!is.null(attr(m,"na.action"))) {
            	temp <- as.integer(attr(m,"na.action"))
            	cross <- cross[-temp]
            	if (length(cross) == nobs) {
                	xgroups <- cross
                	cross <- length(unique(xgroups))
            	}
            else stop("Wrong length for xval")
        	}
        	else stop("Wrong length for xval")
    	}
	
	xpred.cross<-xpred.rpart(c,xval=xgroups)

	if (prune=="mc") {
		xpred.abs.diff.wt<-abs(xpred.cross-Y)*wt
		xerror<-apply(xpred.abs.diff.wt,2,sum)/c$frame$dev[1]
		xerror.within<-matrix(unlist(by(xpred.abs.diff.wt,xgroups,colSums)),cross,nrow(c$cptable),byrow=TRUE)
		xstd<-sqrt(apply(xerror.within,2,var)*(cross-1))/c$frame$dev[1]
	}
	else {
		xpred.err.wt<-(xpred.cross!=Y)*wt
		xerror<-apply(xpred.err.wt,2,sum)/c$frame$dev[1]
		xerror.within<-matrix(unlist(by(xpred.err.wt,xgroups,colSums)),cross,nrow(c$cptable),byrow=TRUE)
		xstd<-sqrt(apply(xerror.within,2,var)*(cross-1))/c$frame$dev[1]


	}

	c$cptable<-cbind(c$cptable,xerror,xstd)
	}
	if (model!=TRUE) {
		c$model<-NULL
	}
	class(c)<-"rpart"
	c
}

