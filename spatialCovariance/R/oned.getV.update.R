precompute.update <- function(info,cat.level=0,aniso=1)
  {
    preCompTime <- proc.time()
    
    ## aniso value has changed.
    ## we need to change the limits of integration

    ## divide out by the old aniso
    info$rowwidth <- info$rowwidth/info$aniso
    info$rowsep <- info$rowsep/info$aniso
    info$indices[,1] <- info$indices[,1]/info$aniso

    ## multiply by the new one
    info$aniso <- aniso
    info$rowwidth <- info$rowwidth*info$aniso
    info$rowsep <- info$rowsep*info$aniso
    info$indices[,1] <- info$indices[,1]*info$aniso

    ## compute new limits of integration
    ## the preLimits are always computed with aniso = 1
    indices <- info$indices.preLimits
    indices <- cbind(indices[[1]],indices[[2]],indices[[3]],indices[[4]],indices[[5]],indices[[6]])
    indices[,1] <- indices[,1]*info$aniso  ## rowsep information via ax

    c.pts <- cbind((indices[,1]-info$rowwidth)^2+(indices[,2]-info$colwidth)^2, indices[,1]^2+(indices[,2]-info$colwidth)^2, indices[,1]^2+indices[,2]^2, (indices[,1]-info$rowwidth)^2+indices[,2]^2, indices[,1]^2+(indices[,2]-info$colwidth)^2, (indices[,1]+info$rowwidth)^2+(indices[,2]-info$colwidth)^2, (indices[,1]+info$rowwidth)^2+indices[,2]^2, indices[,1]^2+indices[,2]^2,  indices[,1]^2+indices[,2]^2, (indices[,1]+info$rowwidth)^2+indices[,2]^2, (indices[,1]+info$rowwidth)^2+(indices[,2]+info$colwidth)^2, indices[,1]^2+(indices[,2]+info$colwidth)^2, (indices[,1]-info$rowwidth)^2+indices[,2]^2, indices[,1]^2+indices[,2]^2, indices[,1]^2+(indices[,2]+info$colwidth)^2, (indices[,1]-info$rowwidth)^2+(indices[,2]+info$colwidth)^2)
    c.pts <- sqrt(c.pts)

    ## expand out the indices so each one corresponds to one integral,
    ## then stick them back together again
    extract.Limits <- function(vals) as.numeric(levels(as.factor(vals)))
    c.pts <- apply(c.pts,1,extract.Limits)
    lengths <- unlist(lapply(c.pts,length))  ## number of limits for each integral
    indices <- cbind(indices,lengths)

    ## expand the current indices to make numerical integration better
    indicesNEW <- NULL
    ptm <- proc.time()
    for(count.i in 1:length(c.pts))
      {
        c.pts[[count.i]] <- matrix(sort(c(c.pts[[count.i]],c.pts[[count.i]]))[2:(lengths[count.i]*2-1)],lengths[count.i]-1,2,byrow=T)
        ## indicator to say which rows are sufficient for analytic results
        ## 1 means evaluate analytic result
        ## 0 means do not, instead give value of 0
        c.pts[[count.i]] <- cbind(c(1,rep(0,dim(c.pts[[count.i]])[1]-1)),c.pts[[count.i]])
        indicesNEW <- rbind(indicesNEW,cbind(matrix(rep(indices[count.i,],lengths[count.i]-1),lengths[count.i]-1,dim(indices)[2],byrow=T),c.pts[[count.i]]))
      }
    if(cat.level>=1) cat("expanding the limits of integration",(proc.time()-ptm)[3],"seconds\n")

    info <- list(rowwidth=info$rowwidth,colwidth=info$colwidth,rowsep=info$rowsep,colsep=info$colsep,nrows=info$nrows,ncols=info$ncols,lengths=lengths,rowReps=info$rowReps,locations=info$locations,indices=indicesNEW,aniso=aniso,indices.preLimits=info$indices.preLimits)

    if(cat.level) cat("precompute update takes",(proc.time()-preCompTime)[3],"seconds\n")

    info
  }
