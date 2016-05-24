optim.indep.norm <-
function(y,status,weight,param,x=NULL,var.list=NULL)
{
    minp <- 0.00000001
    mu <- sigma <- list()
    mu.mat.old <- mu.mat <- matrix(NA,nrow=ncol(weight),ncol=ncol(y))
    sigma.array.old <- sigma.array <- array(0,c(ncol(weight),ncol(y),ncol(y)))
    for(k in 1:ncol(weight))
    {
        mu.mat.old[k,] <- param$mu[[k]]
        sigma.array.old[k,,] <- param$sigma[[k]]
    }
    for(j in 1:ncol(y))
    {
        indivs <- which(is.finite(y[,j]))
        weight.aff <- matrix(weight[status==2,],nrow=sum(status==2),ncol=ncol(weight))
        weight.miss <- matrix(weight[status==0,],nrow=sum(status==0),ncol=ncol(weight))
        weight.aff.j <- matrix(weight.aff[indivs,],nrow=length(indivs),ncol=ncol(weight))
        mu.mat[,j] <- t(weight.aff.j)%*%y[indivs,j]
        mu.mat[,j] <- mu.mat[,j]+apply(weight.miss,2,sum)*mu.mat.old[,j]
        mu.mat[,j] <- mu.mat[,j]/(apply(weight.aff.j,2,sum)+apply(weight.miss,2,sum))
    }
    for(j in 1:ncol(y))
    {
        indivs <- which(is.finite(y[,j]))
        weight.aff <- matrix(weight[status==2,],nrow=sum(status==2),ncol=ncol(weight))
        weight.miss <- matrix(weight[status==0,],nrow=sum(status==0),ncol=ncol(weight))
        weight.aff.j <- matrix(weight.aff[indivs,],nrow=length(indivs),ncol=ncol(weight))
        sigma.array[,j,j] <- diag(t(weight.aff.j)%*%outer(y[indivs,j],mu.mat[,j],"-")^2)
        sigma.array[,j,j] <- sigma.array[,j,j]+apply(weight.miss,2,sum)*(sigma.array.old[,j,j]+(mu.mat.old[,j]-mu.mat[,j])^2)
        sigma.array[,j,j] <- sigma.array[,j,j]/(apply(weight.aff.j,2,sum)+apply(weight.miss,2,sum))
        sigma.array[,j,j] <- pmax(sigma.array[,j,j],minp)
    }
    for(k in 1:ncol(weight))
    {
        mu[[k]] <- mu.mat[k,]
        sigma[[k]] <- matrix(sigma.array[k,,],nrow=ncol(y),ncol=ncol(y))
    }
    param <- list("mu"=mu,"sigma"=sigma)
    param
}

