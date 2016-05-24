esStat<-function(logw=NULL, w.ctrl=NULL, gbm1=NULL, i=1, data,
                  sampw, rule.summary, na.action="level", vars,
                  treat.var, collapse.by.var=FALSE, verbose=FALSE,
                  estimand, multinom){
   
	if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")
	if(estimand=="ATT") {
		if(is.null(gbm1) && is.null(w.ctrl) && is.null(logw))
			stop("No weights given. logw, gbm1, and w.ctrl cannot all be NULL.")
		if(!is.null(rule.summary)) rule.summary <- match.fun(rule.summary)
			w1 <- rep(1/sum(data[,treat.var]==1), nrow(data))
		
		if(is.null(gbm1)){
			if (!is.null(logw)) w.ctrl <- exp(logw)
      		w1[data[,treat.var]==0] <- w.ctrl
      		}
      	else {
      		w <- exp(predict(gbm1,newdata=data[data[,treat.var]==0,],n.trees=i))
      		w1[data[,treat.var]==0] <- w
      		}
      	w1 <- w1*sampw

		# compute effect sizes
		es <- lapply(data[,vars], ps.summary.new2, t=data[,treat.var],
                w=w1, sampw = sampw, get.means=TRUE,
                get.ks=FALSE, na.action=na.action,
                estimand=estimand, multinom = multinom)

		if(collapse.by.var){
			es <- sapply(es,function(x,rule.summary){rule.summary(x$std.eff.sz,na.rm=TRUE)},
                   rule.summary=rule.summary)
            }
		else {
			es <- unlist(sapply(es,function(x){x$std.eff.sz}))
		}
		
		if(!is.null(rule.summary)){
			if(verbose) print(rule.summary(es,na.rm=TRUE))
			return(rule.summary(abs(es),na.rm=TRUE))
			} 
		else{
			return(es)
			}
	}

	if(estimand=="ATE") {
		if (is.null(gbm1) && is.null(w.ctrl) && is.null(logw)) 
        	stop("No weights given. logw, gbm1, and w.ctrl cannot all be NULL.")
		if (!is.null(rule.summary)) 
        	rule.summary <- match.fun(rule.summary)
    	w1 <- rep(1/nrow(data), nrow(data))
    	if (is.null(gbm1)) {
        	if (!is.null(logw)) w.ctrl <- exp(logw)
        w1[data[, treat.var] == 0] <- w.ctrl
    }
    else {
       w <- exp(predict(gbm1, newdata = data, n.trees = i))/(1+exp(predict(gbm1, newdata = data, n.trees = i)))
        w1[data[, treat.var] == 0] <- 1/(1-w[data[, treat.var] == 0])
		w1[data[, treat.var] == 1] <- 1/(w[data[, treat.var] == 1])
    }
    w1 <- w1 * sampw
    es <- lapply(data[, vars], ps.summary.new2, t = data[, treat.var], 
        w = w1, sampw = sampw, get.means = TRUE, get.ks = FALSE, na.action = na.action,
        estimand=estimand, multinom = multinom)
    if (collapse.by.var) {
        es <- sapply(es, function(x, rule.summary) {
            rule.summary(x$std.eff.sz, na.rm = TRUE)
        }, rule.summary = rule.summary)
    }
    else {
        es <- unlist(sapply(es, function(x) {
            x$std.eff.sz
        }))
    }
    if (!is.null(rule.summary)) {
        if (verbose) 
            print(rule.summary(es, na.rm = TRUE))
        return(rule.summary(abs(es), na.rm = TRUE))
    }
    else {
        return(es)
    }
}
}

