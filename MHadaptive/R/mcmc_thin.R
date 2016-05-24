mcmc_thin <-
function(mcmc_object,thin=5)
{
    num_params<-length(mcmc_object$par_names)
    if(num_params > 1)
    {
        ind<-seq(1,length(mcmc_object$trace[,1]),thin)
        mcmc_object$trace<-mcmc_object$trace[ind,]
    }else
    {
        ind<-seq(1,length(mcmc_object$trace),thin)
        mcmc_object$trace<-mcmc_object$trace[ind]
    }
        return(mcmc_object)
}
