

## This is the top level of the survival function which is called. The first if clause is for Brier and it repeats the lower level survival function the number of times 
##corresponding to the length of the brier.vec. It then sums the risk over these times. The second if clause is for IPCW and just runs the lower risk function once.
survival.overall.risk <- function(bas.fx,y,wt,opts){    
    if(opts$loss.fx=="Brier"){ 
        ##Defines the function to be called for each cutpoint in Brier.vec 
        get.each.risk <- function(index, bas.fx, y, wt, opts){    
            ##creates a binary y vector given a specific values of T.star
            T.star=opts$brier.vec[index]
            y.Brier=as.numeric(y[,1]>T.star)
                  
            ##Finds the corresponding weights for this y vector
            wt.Brier=wt[,index]
            ##Gets the risk for this particular y, wt, and bas.fx combination
            get.risk=survival.risk(bas.fx,y.Brier,wt.Brier,opts)
            return(get.risk)
        }
        ##risks is a list with each element corresponding to the risk for each cutpoint in brier.vec   
        risks=lapply(1:length(opts$brier.vec),get.each.risk, bas.fx=bas.fx, y=y, wt=wt, opts=opts)
        ## get the overall risk by summing risks for each individual brier cutpoint. added in multiplication by IBS.wt
        get.overall.risk=sum(unlist(risks) * opts$IBS.wt)
        
    } else if(opts$loss.fx=="IPCW" ){
        ## For IPCW, just called the survival risk function once.
        get.overall.risk=survival.risk(bas.fx,y[,1],wt,opts)
    }
    return(get.overall.risk)
}


## Note that y is just the survival time (without the censoring component).  For IPCW, this will be a continuous variable and for brier it will be binary
## derived from the continuous variable and the cutpoint. 
survival.risk<-function(bas.fx,y,wt,opts){
    
    ## probably don't need this "if" statement
    if (opts$loss.fx=="IPCW"|opts$loss.fx=="Brier"){
        
        ## calculate the average within each basis function vector by first multiplying the weights by the y-values
        numerats<-y*wt
     	denomins<-wt
        
        ## Do some formatting
        if(!is.matrix(bas.fx)){
            if(is.vector(bas.fx)){
                bas.fx<-as.matrix(bas.fx,ncol=1)
            }else bas.fx<-as.matrix(bas.fx,ncol=ncol(bas.fx))
        }
        if (!is.vector(numerats)){
            numerats=as.vector(numerats)
        }
        if (!is.vector(denomins)){
            denomins=as.vector(denomins)
        }
        
        ## Then sum the numerats within each basis function and divide by the total weight in that basis function. betas will be predicted values in each basis function.
        
        betas<- apply(bas.fx*numerats,2,sum)/(apply(bas.fx*denomins,2,sum))

### Added by KL on 10/25/11 to deal with case when all weights in node are 0
        if(length(which(is.na(betas)))>0){
            toReplace=which(is.na(betas))
            for(k in toReplace){
                betas[k]=max(y[which(bas.fx[,k]==1)])		
            }     
        }
        
        ## Get predicted value for each observation by multiplying bas.fx by betas.  
        pred.val <- bas.fx %*% betas
        
        ##calculate loss by comparing y to predicted value
        get.loss <- sum(wt * (y - pred.val) ^ 2)/sum(wt)   
        return(get.loss)
    }
}



 ##### get.at.surv.times #### from AM
"get.at.surv.times"<-
    function(surv.cens, times)
{
                                        #
                                        # surv.cens is an object created by survfit
                                        # times is a vector of times at which you want 
                                        # an estimate of the survival function
                                        #
    nt <- length(times)
    outs <- rep(0, nt)
    survv <- surv.cens$surv
    ns <- length(survv)
    timev <- surv.cens$time
    for(i in 1:nt) {
        if(times[i] < timev[1]) {
            outs[i] <- 1
        }
        else if(times[i] >= timev[ns]) {
            outs[i] <- survv[ns]
        }
        else {
            outs[i] <- survv[timev == max(timev[timev <= times[i]])
                             ][1]
        }
    }
    no <- length(outs[outs == 0])
    outs[outs == 0] <- rep(survv[ns - 1], no)
    return(outs)
}


#### get.init.am ### from AM
"get.init.am"<-
    function(surv.cens, coeff.cox, w, delta, ttilde,deltaMod=delta){
        w <- ifelse(is.na(w),0,w)
        if(!is.matrix(w)){
            w <- matrix(w,nrow=length(w),ncol=1)
        }
      	nn <- length(ttilde)
        coeff.cox <- matrix(coeff.cox,nrow=length(coeff.cox),ncol=1)
        coeff.cox[is.na(coeff.cox)] <- 0
      	linpred <- w %*%coeff.cox
        
      	sum.surv.cond <- surv.cens^(exp(linpred))
        sum.surv.cond <- ifelse(sum.surv.cond<.2,.2,sum.surv.cond)
        if(is.na(min(sum.surv.cond))){
            if(delta[is.na(sum.surv.cond)]==0){
                sum.surv.cond[is.na(sum.surv.cond)] <- 1
            } 
        }
        ##### For IPCW, deltaMod is just delta####
	B <- (deltaMod)/sum.surv.cond
	return(B)
    }


################## Assign Survival Weights ###################################################################

## Note x is only passed in if we are doing IPCW using the cox method to calculate the weights.
assign.surv.wts <- function (x,y,opts,T.star=T.star, cox.vec=cox.vec){
    ##note: default for IPCW in opts will be cox
    cens=y[,2]
    censC=1-cens
    Time=y[,1]

  #case 1 with loss function as IPCW 
    if (opts$loss.fx=="IPCW") {
        ## The method for calculating the weights for IPCW can be cox or KM, first option is cox
        if(opts$wt.method=="Cox"){
            ## if condition means that at least some patients are censored
            if(sum(censC)>0){
                mm <- model.matrix(~.,data.frame(x[,cox.vec])) 
                keep.for.cox <- apply(rbind(apply(data.frame(x[,cox.vec]),1,is.na)),2,sum)
                Time.cox<- Time[keep.for.cox==0]
                censC.cox <- censC[keep.for.cox==0]
                cens.cox <- cens[keep.for.cox==0]
                cov.tr <- as.matrix(mm[,-1])
                surv.cox <- coxph(Surv(Time.cox, censC.cox) ~ cov.tr)
                                        # Get coefficients from cox model
                coeff.cox <- surv.cox$coefficients
                                        # Get baseline survival
                                        # AMM changed 10/15/11 - DUe to change Therneau made in survfit for baseline for mean instead of 0
                                        #surv.cens <- survfit(surv.cox, newdata = data.frame(cov.tr = 0),type = "kaplan-meier")
                surv.cens <- survfit(surv.cox, type = "kaplan-meier")
                                        #  Vector of baseline survival
                basesurv <- get.at.surv.times(surv.cens, Time.cox)
                                        # UG(D) for each subject
                ic0 <- get.init.am(basesurv, coeff.cox, cov.tr, cens.cox, Time.cox)
                
                ## what to do with missing data in \bar{G} now give them a one or 0
                cox.wt1 <- rep(NA,length(keep.for.cox))
                cox.wt1[keep.for.cox==0] <- ic0
                cox.wt1[keep.for.cox==1] <- cens[keep.for.cox==1]
            }
            else cox.wt1<-cens
            return(as.vector(cox.wt1))
            
        }else if(opts$wt.method=="KM"){
                                        #case 2 where loss function is IPCW and wt.method is KM
                                        #create kaplan meier model
            KMfit=survfit(Surv(Time,censC)~1)
                                        #calculate weights by finding the "survival" probability corresponding to the given Brier time
            G=NULL
            for (i in 1:length(Time)){
                                        # FIXME: floating point equality is failing, causing exception
                                        #G[i]<-KMfit$surv[which(KMfit$time==Time[i])]
                G[i]<-KMfit$surv[which.min(abs(KMfit$time - Time[i]))]  # XXX must be better solution
                if(G[i]<.2){G[i]=.2}
            }
                                        #the actual weight is 1/G
            KM.wt1<-1/G
                                        #set those people censored to have weight 0
            KM.wt1[which(cens==0)] <- 0
            return (KM.wt1)     
        }
    }else if(opts$loss.fx=="Brier"){
                                        #now with loss function as brier, we use kaplan meier or cox to calculate weights

                                        # set T.star to be the next smallest Time value
        T.star.closest=Time[(Time-T.star)>=0][which.min((Time-T.star)[(Time-T.star)>=0])]
                                        #For people with event time after T.star, set their brier time equal to T* in order to find the weight according to Graf equation
        Time_Brier=ifelse(Time>T.star,T.star.closest,Time)
        if (opts$wt.method=="Cox"){
                       Time_Brier=ifelse(Time>T.star,T.star,Time)
                       ## if condition means that at least some patients are censored
                       if(sum(censC)>0){
                           mm <- model.matrix(~.,data.frame(x[,cox.vec]))
                           keep.for.cox <- apply(rbind(apply(data.frame(x[,cox.vec]),1,is.na)),2,sum)
                           Time.cox<- Time_Brier[keep.for.cox==0]
                           censC.cox <- censC[keep.for.cox==0]
                           cens.cox <- cens[keep.for.cox==0]
                           cov.tr <- as.matrix(mm[,-1])
                           surv.cox <- coxph(Surv(Time.cox, censC.cox) ~ cov.tr)
                                        # Get coefficients from cox model
                           coeff.cox <- surv.cox$coefficients
                                        # Get baseline survival
                                        # AMM changed 10/15/11 - DUe to change Therneau made in survfit for baseline for mean instead of 0
                                        #surv.cens <- survfit(surv.cox, newdata = data.frame(cov.tr = 0),type = "kaplan-meier")
                           surv.cens <- survfit(surv.cox, type = "kaplan-meier")
                                        #  Vector of baseline survival
                           basesurv <- get.at.surv.times(surv.cens, Time.cox)
                                        # UG(D) for each subject
                           ic0 <- get.init.am(basesurv, coeff.cox, cov.tr, cens.cox, Time_Brier,deltaMod=rep(1,length(cens.cox)))
                           ## what to do with missing data in \bar{G} now give them a one or 0
                           cox.wt1 <- rep(NA,length(keep.for.cox))
                           cox.wt1[keep.for.cox==0] <- ic0
                           cox.wt1[keep.for.cox==1] <- cens[keep.for.cox==1]
                           cox.wt1[which(Time<=T.star& cens==0)] <- 0
                       }else {
                           cox.wt1<-cens
                       }
                       return(as.vector(cox.wt1))
                   }else if(opts$wt.method=="KM"){
                                        #create kaplan meier model
                       KMfit=survfit(Surv(Time,censC)~1)
                                        #calculate weights by finding the "survival" probability corresponding to the given Brier time
                       G=NULL
                       for (i in 1:length(Time)){
                                        # FIXME: floating point equality is failing, causing exception
                                        #G[i]<-KMfit$surv[which(KMfit$time==Time_Brier[i])]
                           G[i]<-KMfit$surv[which.min(abs(KMfit$time - Time_Brier[i]))]  # XXX must be better solution
                           if(G[i]<.2){G[i]=.2}
                       }
                                        #the actual weight is 1/G
                       KM.wt1<-1/G
                                        #set those people censored before T.star to have weight 0
                       KM.wt1[which(Time<=T.star& cens==0)] <- 0
                       return (KM.wt1)  
                   }        
    }
}


### We will only go into this function for the brier loss function. We will reassign the coefficients in the model constructed on the
### training set by recalculating the predicted values based on the binary y and the weight for a given time point.
calculate.brier.risk <- function(model, x, y, wt, x.test, y.test, wt.test, opts){
	
    test.set.risk.DSA=array(0,length(model$coefficients))
                                        #brier.vec is a vector containing the times we are looking at for brier
    brier.vec=opts$brier.vec
    for (i in 1:length(brier.vec)){
        T.star=brier.vec[i]
                                        #create binary versions of y based on this particular T.star value
        y.Brier=as.numeric(y[,1]>T.star)
        y.test.Brier=as.numeric(y.test[,1]>T.star)
               
 	
                                        #find weights using this particular T.star value 
        wt.Brier=wt[,i]
        wt.test.Brier=wt.test[,i]
                                        #using these weights and y-values, modify the field ty$coefficients to represent the new predicted values for each node 
        fun <- function(k) get.coefs(k, dat=x, y=y.Brier, wt=wt.Brier)$coef
        get.coef.from.training <- lapply(model$IkPn, fun)
        model$coefficients <- get.coef.from.training
        
                                        #use this updated ty, to predict values based on the test set
        pred.test.DSA <- predict.dsa(model,x.test)
        
                                        #calculate the loss and make cv.risk the sum over all values in brier.vec
        tmp <- as.vector(wt.test.Brier )* (as.vector(y.test.Brier) - pred.test.DSA) ^ 2
        test.set.risk.DSA <- test.set.risk.DSA + opts$IBS.wt[i] * apply(tmp, 2, sum) / sum(wt.test.Brier)
    }
    return (test.set.risk.DSA)
}
