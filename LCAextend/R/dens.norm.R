dens.norm <-
function(y.x,param,var.list=NULL)
{
    y <- y.x[1:length(param$mu[[1]])]   
    if(length(y.x)==length(param$mu[[1]])) x <- NULL
    else x <- y.x[(length(param$mu[[1]])+1):length(y.x)]
    res <- rep(1,times=length(param$mu))
    if(all(is.na(y))) stop("there is an affected subject with no measurements")
    miss <- which(is.na(y))
    for(i in 1:length(param$mu)) 
    {
        if(length(miss)==0) res[i] <- dmvnorm(y,param$mu[[i]],param$sigma[[i]])
        else
        {
            if(length(y[-miss])==1) res[i] <- dnorm(y,param$mu[[i]],param$sigma[[i]])
            else res[i] <- dmvnorm(y[-miss],param$mu[[i]][-miss],param$sigma[[i]][-miss,-miss])
        }
    }
    res
}

