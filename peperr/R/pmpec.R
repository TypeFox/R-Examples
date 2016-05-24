
pmpec <- function(object, response=NULL, x=NULL, times, model.args=NULL, type=c("PErr","NoInf"), 
                  external.time=NULL, external.status=NULL, data=NULL) {
    require(survival)
    type <- match.arg(type)

    if(is.null(response)&&is.null(x)) {
       if(!is.null(data)){
          time <- data$time
          status <- data$status
          x <- as.matrix(data[,!(names(data) %in% c("time", "status")), drop=FALSE])
       } else {
          stop("Please pass either arguments 'response' and 'x' or 'data'")
       }
    } else {
    time <- response[, "time"]
    status <- response[, "status"]
    }

    fit.time <- time
    fit.status <- status
    
    if (!is.null(external.time) && !is.null(external.status)) {
        fit.time <- external.time
        fit.status <- external.status
    }

    censdata <- data.frame(time=fit.time,status=ifelse(fit.status == 0,1,0))
    cens.fit <- survival::survfit(Surv(time,status) ~ 1,data=censdata)

    eval.cens.prob <- summary(cens.fit,times=times)$surv
    
    if (is.null(external.time) || is.null(external.status)) {
        sortsurv <- summary(cens.fit,times=sort(time))$surv
        invcensprob <- rep(NA,length(status))
        invcensprob[order(time)] <- c(1,sortsurv[1:length(sortsurv)-1])
    } else {
        sort.time <- sort(fit.time)
        sortsurv <- c(1,summary(cens.fit,times=sort.time)$surv)
        sort.time <- c(-999,sort.time)
        sortsurv <- c(1,sortsurv)
        invcensprob <- unlist(lapply(time,function(actual.time) sortsurv[which(sort.time >= actual.time)[1]-1]))
    }

    mod.time <- time
    if (any(status != 0 | status != 1)) {
        mod.time[status > 1] <- max(time,times) + 1
    }

#     status.mat <- matrix(times,length(time),length(times),byrow=TRUE) < mod.time
# 
#     probmat <- do.call("predictProb", c(list(object=object, response=response, x=x ,times=times),
#          model.args))
# 
#     weightmat <- t(t(status.mat) / eval.cens.prob) + (1 - status.mat) * matrix(status,length(status),length(times)) * matrix(1/invcensprob,length(status),length(times))

    status.mat <- matrix(times, length(time), length(times), byrow = TRUE) < mod.time
    weight.status.mat <- matrix(times, length(time), length(times), byrow = TRUE) < time
    probmat <- do.call("predictProb", c(list(object = object, response = response, x = x, times = times), model.args))
    weightmat <- t(t(weight.status.mat)/eval.cens.prob) + (1 - weight.status.mat) * matrix(status != 0, length(status), length(times)) * matrix(1/invcensprob,  length(status), length(times))

    if (type == "PErr") {
        return(apply(weightmat*(status.mat - probmat)^2,2,mean))
    } else {
        
        #   vectorized version (that actually seems to be slower!)
        # return(apply(matrix(as.vector(t(weightmat)[,rep(1:nrow(data),rep(nrow(data),nrow(data)))])*(
        #                     as.vector(t(status.mat)[,rep(1:nrow(data),rep(nrow(data),nrow(data)))]) - as.vector(t(probmat)))^2
        #                     ,ncol=length(times),byrow=TRUE),2,mean))
        
        res <- .C("noinf",
                  as.double(as.vector(weightmat)),
                  as.integer(as.vector(status.mat)),
                  as.double(as.vector(probmat)),
                  as.integer(nrow(x)),
                  as.integer(length(times)),
                  err=double(length(times))
                  )
        
        return(res$err)
                
        # t.status.mat  <- t(status.mat)
        # t.weightmat <- t(weightmat)        
        # errmat <- matrix(unlist(lapply(1:nrow(probmat),function(i) apply(t.weightmat*(t.status.mat - probmat[i,])^2,1,mean))),nrow(probmat),length(times),byrow=TRUE)
        # return(apply(errmat,2,mean))
    }
}

