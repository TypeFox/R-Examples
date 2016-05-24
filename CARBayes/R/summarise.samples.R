summarise.samples <- function(samples, columns=NULL, quantiles=0.5, exceedences=NULL)
{
#### Compute exceedence probabilities and posterior quantiles
#### for  a set of given samples

#### Check for initial errors
     if(class(samples) != "mcmc") stop("The samples object is not a mcmc type.", call.=FALSE)    

     if(!is.null(exceedences))
     {
          if(!is.numeric(exceedences) ) stop("the exceedences argument is not numeric.", call.=FALSE)          
     }else
     {
     }
     
     if(!is.null(quantiles))
     {
          if(!is.numeric(quantiles)) stop("the quantiles argument is not numeric.", call.=FALSE)
          if(min(quantiles) < 0) stop("some of the quantiles are less than zero.", call.=FALSE)
          if(max(quantiles) > 1) stop("some of the quantiles are greater than one.", call.=FALSE)
     }else
     {
     }
     
     

     if(is.null(columns))
     {
     N.all <- ncol(samples)
     columns <- 1:N.all     
     }else
     {
     N.total <- ncol(samples)
     N.all <- length(columns)
          if(!is.numeric(columns)) stop("the columns argument is not numeric.", call.=FALSE)
          if(sum(floor(columns)==ceiling(columns)) < N.all) stop("the columns argument is not all integers.", call.=FALSE)
          if(min(columns)<1) stop("some of the columns chosen are less than 1.", call.=FALSE)
          if(max(columns)>N.total) stop("some of the columns chosen are greater than the size of the samples object.", call.=FALSE)
          }
     

#### Compute the summaries     
n.sample <- nrow(samples)
     if(missing(exceedences))
     {
     N.exceedence <- 0     
     }else
     {
     N.exceedence <- length(exceedences)
     }
          
     if(missing(quantiles))
     {
     N.quantiles <- 0     
     }else
     {
     N.quantiles <- length(quantiles)
     }

     if(N.exceedence>0 & N.quantiles>0) 
     {
     posterior.exceedence <- array(NA, c(N.all,N.exceedence))
     posterior.quantile <- array(NA, c(N.all,N.quantiles))
     colnames(posterior.exceedence) <- exceedences
     colnames(posterior.quantile) <- quantiles
             
          for(j in 1:N.all)
          {
          mcmc <- samples[ ,columns[j]]
          posterior.quantile[j, ] <- quantile(mcmc, quantiles)
               for(k in 1:N.exceedence)
               {
               posterior.exceedence[j, k] <- length(which(mcmc>exceedences[k]))/n.sample       
               }
          }
     }else if(N.exceedence==0 & N.quantiles>0)
     {
     posterior.exceedence <- NULL
     posterior.quantile <- array(NA, c(N.all,N.quantiles))
     colnames(posterior.quantile) <- quantiles
             
           for(j in 1:N.all)
           {
            mcmc <- samples[ ,columns[j]]
            posterior.quantile[j, ] <- quantile(mcmc, quantiles)
           }    
     }else if(N.exceedence>0 & N.quantiles==0)
     {
     posterior.exceedence <- array(NA, c(N.all,N.exceedence))
     posterior.quantile <- NULL
     colnames(posterior.exceedence) <- exceedences
             
           for(j in 1:N.all)
           {
            mcmc <- samples[ ,columns[j]]
                for(k in 1:N.exceedence)
                {
                posterior.exceedence[j, k] <- length(which(mcmc>exceedences[k]))/n.sample       
                }
           }
     }else
     {
     posterior.quantile <- NULL
     posterior.exceedence <- NULL
     }
       
results <- list(quantiles=posterior.quantile, exceedences=posterior.exceedence)  
return(results)
}
 
