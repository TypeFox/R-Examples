##################################################################
train_predict_mix <- function(
         test,train,k,    
         theta0=0,alpha.shape=0.5,alpha.rate=5,no.alpha=30,
	 common.alpha=FALSE,no.alpha0=100,
         mc.iters=200,iters.labeltheta=10,
	 iters.theta=20,width.theta=0.1,
         correction=TRUE,no.theta.adj=30,approxim=TRUE,
         pred.start=100)
	 
{   theta.range <- c(theta0,1-theta0)
    p <- ncol(train) - 1
    no.test <- nrow(test)
    if(k < p)
    {   info.sel <- selfth.abscor(k,x=train[,-1],y=train[,1])
        fth.sel <- info.sel$fth+1
        gamma.sel <- info.sel$abscors[k]    
    }
    else 
    {  fth.sel <- 1:p+1
       gamma.sel <- 0
       correct <- FALSE
    }
    spl <- training(
             train[,c(1,fth.sel),drop=FALSE],gamma.sel,p-k,
             mc.iters,iters.theta,width.theta,iters.labeltheta,
             theta.range,alpha.shape,alpha.rate,no.alpha,
             correction,no.theta.adj,common.alpha,approxim,no.alpha0)
    
    prob.pred <- Predict.mix(test[,c(1,fth.sel),drop=FALSE],
                             spl,pred.start)
    
    #report result
    states.pred <- 1*(prob.pred>0.5) 
    wrong <- 1*(states.pred != test[,1])
    error.rate <- mean(wrong) 
    # calculate average minus log probs
    aml <- -mean(dbinom(test[,1],1,prob.pred,log=TRUE))
    # tabulate the predictive probs
    summary.pred <- present(prob.pred,test[,1]) 
    # calculate mean square error
    mse <- mean((prob.pred-test[,1])^2)   
    
    #change the names in spl
    names(spl)[names(spl)=="r"] <- "alpha"
    names(spl)[names(spl)=="r1"] <- "alpha0"
    names(spl)[names(spl)=="phi"] <- "theta"
    
    #output
    c( list(
    	    prediction=cbind(true=test[,1],pred=states.pred,
	                  prob.pred=prob.pred,wrong),
	    aml=aml,error.rate=error.rate,mse=mse,
            summary.pred=summary.pred,
	    
	    features.selected = c(feature=fth.sel,cutoff=gamma.sel)),
       spl)
}

#spl elements: 
#list(label=label,I1=I1,I2=I2,N1=N1,N2=N2,phi=phi,r=r,r1=r1,alpha_set= 
#    rspace,alpha0_set=rspace1)
##################################################################
#calculating the absolute correlation between x and y
abs_cor <- function(x,y)
{
    if(sum(abs(x-mean(x)))==0 || sum(abs(y-mean(y)))==0 )
	return(0)
    else
	return(abs(cor(y,x)))
}
#this function selects features from x by the absolate correlation with y
selfth.abscor <- function(k,x,y)
{  
   abscors <- apply(x,2,abs_cor,y)
   info.sort <- sort(abscors,decreasing=TRUE,index.return=TRUE)
   list(fth=info.sort$ix[1:k], abscors=info.sort$x[1:k])
}
