

# computes KS statistics from weights given in a gbm object, log(w), or control weights
ksStat<-function(logw=NULL,
                  w.ctrl=NULL,
                  gbm1=NULL,
                  i=1,
                  data,
#                  sampw=rep(1,nrow(data)),
                  sampw,                  
#                  rule.summary=mean,
                  rule.summary,                  
                  na.action="level",
                  vars,
                  treat.var,
                  collapse.by.var=FALSE,
                  verbose=FALSE, estimand, multinom)
{

if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")

if(estimand=="ATT") {
	if(is.null(gbm1) && is.null(w.ctrl) && is.null(logw))
		stop("No weights given. logw, gbm1, and w.ctrl cannot all be NULL.")
	if(!is.null(rule.summary)) rule.summary <- match.fun(rule.summary)
	
	w1 <- rep(1/sum(data[,treat.var]==1), nrow(data))
	
	if(is.null(gbm1)){
      if (!is.null(logw)) w.ctrl<-exp(logw)
      	w1[data[,treat.var]==0] <- w.ctrl
      } 
      else {
      	w <- exp(predict(gbm1,newdata=data[data[,treat.var]==0,],
                           n.trees=i))
        w1[data[,treat.var]==0]<- w
   }

   # compute KS statistics
   
   w1 <- w1*sampw
   ks <- lapply(data[,vars], ps.summary.new2,
                t=data[,treat.var],
                w=w1,
                get.means=FALSE,
                sampw = sampw,
                get.ks=TRUE,
                na.action=na.action,
                estimand=estimand, multinom=multinom)
   ks <- lapply(ks, function(x){x$ks})

   if(collapse.by.var) ks <- sapply(ks,max)
   else ks <- unlist(ks)

   if(!is.null(rule.summary)){
      if(verbose) print(rule.summary(ks))
      return(rule.summary(ks))
   } 
   else{
      return(ks)
   }
}

if(estimand=="ATE") {
    if (is.null(gbm1) && is.null(w.ctrl) && is.null(logw)) 
        stop("No weights given. logw, gbm1, and w.ctrl cannot all be NULL.")
    if (!is.null(rule.summary)) 
        rule.summary <- match.fun(rule.summary)
    w1 <- rep(1/nrow(data), nrow(data))
    if (is.null(gbm1)) {
        if (!is.null(logw)) {
            w.ctrl <- exp(logw)
        }
        w1[data[, treat.var] == 0] <- w.ctrl
    }
    else {
        w <- exp(predict(gbm1, newdata = data, n.trees = i))/(1+exp(predict(gbm1, newdata = data, n.trees = i)))
        w1[data[, treat.var] == 0] <- 1/(1-w[data[, treat.var] == 0])
	w1[data[, treat.var] == 1] <- 1/(w[data[, treat.var] == 1])
    }
    w1 <- w1 * sampw
    ks <- lapply(data[, vars], ps.summary.new2, t = data[, treat.var], sampw = sampw, 
        w = w1, get.means = FALSE, get.ks = TRUE, na.action = na.action, estimand=estimand, multinom = multinom)
    ks <- lapply(ks, function(x) {
        x$ks
    })
    if (collapse.by.var) 
        ks <- sapply(ks, max)
    else ks <- unlist(ks)
    if (!is.null(rule.summary)) {
        if (verbose) 
            print(rule.summary(ks))
        return(rule.summary(ks))
    }
    else {
        return(ks)
    }

}


}