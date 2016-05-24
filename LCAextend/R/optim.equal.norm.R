optim.equal.norm <-
function(y,status,weight,param,x=NULL,var.list=NULL)
{
    mu <- sigma <- list()
    mu.mat.old <- mu.mat <- matrix(NA,nrow=ncol(weight),ncol=ncol(y))
    sigma.array.old <- sigma.array <- array(NA,c(ncol(weight),ncol(y),ncol(y)))
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
    for(j1 in 1:ncol(y)) for(j2 in 1:ncol(y))
    {
        indivs <- intersect(which(is.finite(y[,j1])),which(is.finite(y[,j2])))
        weight.aff <- matrix(weight[status==2,],nrow=sum(status==2),ncol=ncol(weight))
        weight.miss <- matrix(weight[status==0,],nrow=sum(status==0),ncol=ncol(weight))
        weight.aff.j1j2 <- matrix(weight.aff[indivs,],nrow=length(indivs),ncol=ncol(weight))
        sigma.array[,j1,j2] <- diag(t(weight.aff.j1j2)%*%(outer(y[indivs,j1],mu.mat[,j1],"-")*outer(y[indivs,j2],mu.mat[,j2],"-")))
        sigma.array[,j1,j2] <- sigma.array[,j1,j2]+apply(weight.miss,2,sum)*(sigma.array.old[,j1,j2]+(mu.mat.old[,j1]-mu.mat[,j1])*(mu.mat.old[,j2]-mu.mat[,j2]))
        sigma.array[,j1,j2] <- sigma.array[,j1,j2]/(apply(weight.aff.j1j2,2,sum)+apply(weight.miss,2,sum))
    }
    if(ncol(y)==1)
    {
        for(k in 1:ncol(weight))
        {
            mu[[k]] <- mu.mat[k,]
            sigma[[k]] <- matrix(sigma.array[k,,],nrow=ncol(y),ncol=ncol(y))
        }
    }
    else
    {
        diags <- apply(apply(sigma.array,1,diag),2,sum)
        no.diags <- (apply(sigma.array,1,sum)-diags)/(ncol(y)*(ncol(y)-1))
        diags <- diags/ncol(y)
        for(k in 1:ncol(weight))
        {
            mu[[k]] <- mu.mat[k,]
            sigma[[k]] <- matrix(no.diags[k],nrow=ncol(y),ncol=ncol(y))
            diag(sigma[[k]]) <- diags[k]
        }
    }
    param <- list("mu"=mu,"sigma"=sigma) 
    param
}

