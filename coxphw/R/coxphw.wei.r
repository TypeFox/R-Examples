qeq<-function(a,b,tol=1e-6){
 abs(a-b)<tol
 }
 

coxphw.wei <- function
(
 formula,
 data,
 obj,
 censcorr,
 normalize,
 breslow,
 prentice,
 taroneware,
 scale.weights,
 trunc.weights,
 id,
 cov.name,
 NFP
 )
### calculate weights for coxphw
### 2008-05
### if fp( ) terms are present, then all effects will be weighted (provided any effect is specified in breslow, prentice of tarone; 
###         it does not matter which effects are specified there)
{
        NTDE <- obj$NTDE
        k <- ncol(obj$mm1)              # number covariates
        n <- nrow(data)
        
        weights <- matrix(1, n, k+NTDE) # weight matrix
        w.raw <- weights[,1]
        w <- w.raw
        formulaAll <- as.formula(paste(as.character(formula)[2], "1", sep="~"))
        ## fit <- survfit(formulaAll, data)
        my.survfit<-getFromNamespace("survfit","survival")
        fit <- my.survfit(Surv(obj$resp[,1],obj$resp[,2],obj$resp[,3])~1, data)
        event <- event1 <- obj$resp[, 3]
        event1[-1][diff(obj$resp[, 2]) == 0] <- 0   # first of simultan. events
        ties <- (diff(c(-999, obj$resp[, 2])) == 0) ## change GH 080411
        ## ties <- (diff(c(-999, obj$resp[, 2])) == 0) & event

   if(censcorr) {
          time.obs.km <- (-(1:max(id)))
          cens.obs.km <- ((1:max(id)))
          for(id.i in 1:max(id)){
            time.obs.km[id.i]<-max(obj$resp[id==id.i,2])
            cens.obs.km[id.i]<-1-max(obj$resp[id==id.i,3])
           }
           obskm<-my.survfit(Surv(time.obs.km,cens.obs.km)~1)
           obskm$surv<-cbind(1,obskm$surv)
           obskm$time<-cbind(0,obskm$time)
           if(is.na(breslow)&is.na(prentice)&is.na(taroneware)) {
             w.obskm<-w
             for(i in 1:n) {
              w.obskm[i] <- 1/min(obskm$surv[obskm$time<obj$resp[i,2]])
              w[i]<-w.raw[i]*w.obskm[i]
             }
            }
            if(normalize)    w <- w * sum(event1) / sum((w * event1))
            weights<-matrix(w,length(w),k+NTDE)
          }
        if(!is.na(taroneware)) {
                if(taroneware != TRUE) {
                    formulaTW <- as.formula(paste(as.character(formula)[2],
                                              as.character(taroneware)[2], sep="~"))
                    objTW <- decomposeSurv(formulaTW, data, sort=TRUE)
                    cov.nameTW <- objTW$covnames[objTW$ind]
                } else cov.nameTW <- cov.name
                ## n at risk with counting process style
                w <- (n:1)   #fit$n.risk
                for (i in 1:n){
                        if(event1[i]==1){
                                w[i] <- sum((obj$resp[,1] < obj$resp[i,2]) & (obj$resp[,2] >= obj$resp[i,2]))
                        }
                        else w[i]<-0
                }
                w <- w ^0.5
                
                w[!event] <- 0
                w[ties] <- w[c(ties[-1], FALSE)]
                w.raw<-w
                if(censcorr) {
                        w.obskm<-w
                        for(i in 1:n) {
                                w.obskm[i] <- 1/min(obskm$surv[obskm$time<obj$resp[i,2]])
                                w[i]<-w.raw[i]*w.obskm[i]
                        }
                }

                if(normalize)    w <- w * sum(event1) / sum((w * event1))
                if(NFP==0) weights[, match(cov.nameTW, cov.name)] <- w else weights<-matrix(w,length(w),k+NTDE)
#                weights[, match(cov.nameTW, cov.name)] <- w
        }
        if(!is.na(breslow)) {
                if(breslow!=TRUE) {
                 formulaB <- as.formula(paste(as.character(formula)[2],
                                             as.character(breslow)[2], sep="~"))
                 objB <- decomposeSurv(formulaB, data, sort=TRUE)
                 cov.nameB <- objB$covnames[objB$ind]
                } else cov.nameB <- cov.name
                ## n at risk with counting process style
                ## count n rows with time1<time2[i] and time2>=time2[i]
                w <- (n:1)   #fit$n.risk
                for (i in 1:n){
                        if(event[i]==1){
                                w[i] <- sum((obj$resp[,1] < obj$resp[i,2]) & (obj$resp[,2] >= obj$resp[i,2]))
                        }
                        else w[i]<-0
                }
                ##  w <- w * sum(w * event1) / sum((w * event1) ^ 2)
                w[!event] <- 0
                w[ties] <- w[c(ties[-1], FALSE)]
                w.raw<-w
                if(censcorr) {
                        w.obskm<-w
                        for(i in 1:n) {
                                w.obskm[i] <- 1/min(obskm$surv[obskm$time<obj$resp[i,2]])
                                w[i]<-w.raw[i]*w.obskm[i]
                        }
                }
                if(normalize)    w <- w * sum(event1) / sum((w * event1))
                if(NFP==0) weights[, match(cov.nameB, cov.name)] <- w else weights<-matrix(w,length(w),k+NTDE)
        }
        if(!is.na(prentice)) {
                if(prentice != TRUE) {
                    formulaP <- as.formula(paste(as.character(formula)[2],
                                             as.character(prentice)[2], sep="~"))
                    objP <- decomposeSurv(formulaP, data, sort=TRUE)
                    cov.nameP <- objP$covnames[objP$ind]
                } else cov.nameP <- cov.name
                ## w <- c(1, fit$surv[-length(fit$surv)])[cumsum(1 - ties)] #S_t
                w.sorted <- c(1, fit$surv[-length(fit$surv)]) #S_t
                w<-obj$resp[,3]
                for(i in 1:n) {
                        w[i]<-w[i]*w.sorted[qeq(fit$time,obj$resp[i,2])]
                }
                w[!event] <- 0
                w[ties] <- w[c(ties[-1], FALSE)]
                w.raw<-w
                if(censcorr) {
                        w.obskm<-w
                        for(i in 1:n) {
                                w.obskm[i] <- 1/min(obskm$surv[obskm$time<obj$resp[i,2]])
                                w[i]<-w.raw[i]*w.obskm[i]
                        }
                }
                if(normalize)    w <- w * sum(event1) / sum((w * event1))
                if(NFP==0) weights[, match(cov.nameP, cov.name)] <- w else weights<-matrix(w,length(w),k+NTDE)
#               weights[, match(cov.nameP, cov.name)] <- w
        }
        
        if(!censcorr)
          w.obskm <- rep(1,n)
        
        ## calculate weight matrix
        ### weight truncation (30-04-2009)
        quant <- function(x)  { as.vector(quantile(x, probs=trunc.weights)) }
        quant.w<-quant(w)        
        w[w>quant.w]<-quant.w
        w.matrix <- cbind(time=obj$resp[,2], w.raw, w.obskm, w)[obj$resp[,3] != 0,]
        orderw <- rank(w.matrix[,1])
        w.matrix <- w.matrix[orderw, , drop=FALSE]

        ## are all weight columns the same?
        same <- all(weights == weights[, 1])
       truncat <- apply(weights, 2, quant)
  
       for (i.wt in 1:ncol(weights))
       {
         weights[weights[,i.wt]>truncat[i.wt], i.wt] <- truncat[i.wt]
       }
       
        ## number of weighted variables
        NGV <- sum(weights[1, ] != 1 | weights[n, ] != 1)

        weights <- weights * scale.weights
        
        ## return list
        list(
             weights=weights,
             w.matrix=w.matrix,
             const=same,                # are all weights the same?
             NGV=NGV                    # number of weighted variables
             )
}
