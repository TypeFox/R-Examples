plotMH <-
function(mcmc_object,correlogram=TRUE)
{
    par(pty='s',mar=c(3.1,4,3,1),mgp=c(2,1,0))
    num_params<-length(mcmc_object$par_names)
    par(mfrow=c(num_params,2))
    
    if(num_params > 1)
    {
        for(i in 1:num_params)
        {
            par<-mcmc_object$par_names[i]
            hist(mcmc_object$trace[,i],breaks=30,probability=TRUE,xlab=par,ylab='Posterior Density',main=paste('Posterior distribution of ',par) )
            plot(mcmc_object$trace[,i],type='l',main=paste('Trace of ',par),xlab='Iteration',ylab=par)
        }
      
        if(correlogram)
        {
            par(ask=TRUE)
            par(mfrow=c(num_params,num_params),mar=c(2.5,2,1,1),mgp=c(2,1,0))
            for(i in 1:num_params)
            {
                for(j in 1:num_params)
                {
                    if(i == j)
                    {
                        frame()
                        text(0.5,0.5,mcmc_object$par_names[i],cex=2)
                    }
                    if(i != j)
                        plot(mcmc_object$trace[,j],mcmc_object$trace[,i],xlab='',ylab='',main='')
                }
            }
        }
    }
    if(num_params == 1)
    {
            par<-mcmc_object$par_names
            hist(mcmc_object$trace,breaks=30,probability=TRUE,xlab=par,ylab='Posterior Density',main=paste('Posterior distribution of ',par) )
            plot(mcmc_object$trace,type='l',main=paste('Trace of ',par),xlab='Iteration',ylab=par)
    }
    par(ask=FALSE)
}

