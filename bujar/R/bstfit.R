### Different methods for BJ estimation
###fit boosting or penalized regression within a BJ iteration 
bstfit <- function(tuning, x, ynew, nu, mselect, mstop2, twin, center, interaction, degree, learner, l2, nfold, n.cores, cv, tuningSwitch, k, trace, gamma, lambda, lamb, whichlambda, method, rng){
    dat1 <- as.data.frame(cbind(ynew,x))
    x <- as.matrix(x)
    if(learner=="linear.regression"){
        ctrl <- boost_control(nu=nu,mstop=mselect)
	#ctrl <- boost_control(nu=nu,mstop=mstop)
        dat1.glm <- glmboost(ynew~x,center=center, control=ctrl)
    }
    if(learner=="pspline"){
        ctrl <- boost_control(nu=nu,mstop=mselect)
        dat1.glm <- gamboost(ynew~., data=dat1,control=ctrl)
    }
    if(learner=="mars"){
        penalty <- ifelse(degree==2,3,2) #This is unneccessary since this is the default in the mars function. We follow the advice in [8] to set the cost for each basis function optimization
                                        # to be 3 in the MARS for two-way interaction models. (see Friedman 1991, and Yi Lin and Helen Zhang, AOS, 2006)
        if(tuningSwitch){
            dat1.glm <- mars(x=as.data.frame(x),y=ynew,degree=degree,penalty=penalty)
            nk <- dat1.glm$nk
        }
        else   dat1.glm <- mars(x=as.data.frame(x),y=ynew,degree=degree,penalty=penalty,nk=nk,prune=FALSE)   #should I add "prune=FALSE" here
        dat1.glm <- mars.to.earth(dat1.glm, trace=FALSE)
    }
    if(learner=="tree"){
        datfit <- as.data.frame(list(y=ynew,x=x))
        colnames(datfit) <- c("y",colnames(x));
        if(!is.null(rng))
	set.seed(rng)
	dat1.glm <- gbm(y~.,data=datfit,distribution="gaussian",interaction.depth=degree,shrinkage=nu,cv.folds=nfold,n.trees=mselect,verbose=FALSE, n.cores=n.cores)
	###############boosting trees################
        if(tuning && cv && tuningSwitch){
            mselect <- gbm.perf(dat1.glm,plot.it=FALSE,method="cv")
            if(!is.null(rng))
	    set.seed(rng)
            dat1.glm <- gbm(y~.,data=datfit,distribution="gaussian",interaction.depth=degree,shrinkage=nu,n.trees=mselect,verbose=FALSE)
        }
        if(trace) cat("k=",k," number of trees utilized",mselect,"\n")
    }
                                        #    if(learner=="acosso"){
                                        #      best.pow <- wt.pow
                                        #      dat1.glm <- acosso(X=as.data.frame(x),y=ynew,order=degree,wt.pow=best.pow)  #use 2 from Section 5, line 4 of Adaptive COSSO; 0 for COSSO
                                        #    }
    if(learner=="enet"){
        dat1.glm <- enet(x=as.matrix(x), y=ynew, lambda=lamb)
    }
    if(learner %in% c("enet2", "mnet","snet")){
        dat1.glm <- glmreg(x=as.matrix(x), y=ynew, penalty=l2, gamma=gamma, lambda=lambda)
    }
######################################################################################
### Choose tuning parameter for penalized regression methods, NOT for learner "enet", "mars","acosso" ####
###################################################################################### 
    if(learner %in% c("enet2", "mnet","snet")){
        if(tuning){
            cv.fit <- cv.glmreg(x=as.matrix(x), y=ynew, penalty=l2, gamma=gamma, nfolds=nfold, plot.it=FALSE)
            mselect <- cv.min <- cv.fit$lambda.which
### needs something with fixed lambda
        }
      else mselect <- cv.min <- whichlambda
    }
    if(learner%in% c("linear.regression","pspline") && tuning){
        if(!cv){
            if(learner=="pspline") df <- "trace"
            else df <- "actset"
            aic <- AIC(dat1.glm,method=method,df=df) #model selection
            if(trace)
                cat("aic=",aic, "mstop(aic)=",mstop(aic),"\n")
	    if(mstop(aic) > 1)
            mselect <- mstop(aic)
            dat1.glm <- dat1.glm[mselect]
        }
      else {
          if(learner %in% c("linear.regression","pspline") ){                      
              file <- NULL
              mselect <- cvmboost(obj=dat1.glm,nrow(dat1),nfold=nfold,figname=file);
              if(trace) cat("CV selects mstop", mselect,"\n")
              dat1.glm <- dat1.glm[mselect]
          }
      }
    }
####### twin boosting #######
    if(twin && learner=="linear.regression"){
        xselect.init = unique(dat1.glm$xselect()) - 1 ### it seems to treat intercept as first variable
        coefir=coef(dat1.glm, which = 1:length(variable.names(dat1.glm)))[-1]
        dat1.glm <- bst(x, ynew, learner=l2, ctrl = bst_control(twinboost=TRUE,
                                                                coefir=coefir, xselect.init = xselect.init, mstop=mstop2, center=TRUE, intercept=TRUE))
    }	    
    if(twin && learner %in% c("pspline", "tree")){
        if(learner=="pspline"){
            f.init <- predict(dat1.glm)
            xselect.init = unique(dat1.glm$xselect())
        }	   else{
            f.init <- dat1.glm$fit
            xselect <- summary(dat1.glm,order=FALSE,plotit=FALSE)[,2]
            xselect.init <- which(xselect > 0)
        }
        dat1.glm <- bst(x, ynew, learner=l2, ctrl = bst_control(twinboost=TRUE,
                                                                f.init=f.init, xselect.init = xselect.init, mstop=mstop2, center=TRUE, intercept=TRUE))
    }
    RET <- list(dat1.glm=dat1.glm, mselect=mselect)
}
