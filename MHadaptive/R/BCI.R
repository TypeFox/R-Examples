BCI <-
function(mcmc_object,interval=c(0.025,0.975))
{
    num_params<-length(mcmc_object$par_names)
    if(num_params > 1)
    {
        bci<-array(dim=c(num_params,length(interval)))
        post_mean<-numeric(num_params)
        for(i in 1:num_params)
        {
            bci[i,]<-quantile(mcmc_object$trace[,i],probs=interval)
            post_mean[i]<-mean(mcmc_object$trace[,i])
        }
        rownames(bci)<-mcmc_object$par_names
        colnames(bci)<-as.character(interval)
    }
    if(num_params == 1)
    {  
            bci<-quantile(mcmc_object$trace,probs=interval)
	    post_mean<-mean(mcmc_object$trace)
    }
    return(cbind(bci,post_mean))
}

