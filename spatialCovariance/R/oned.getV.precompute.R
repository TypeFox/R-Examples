## April 24 2004

precompute <- function(nrows,ncols,rowwidth,colwidth,rowsep,colsep,cat.level=0)
  {
    if(nrows==1 && ncols==1) {
      ## this is for when there is one row and one column
      ## rarely used but there for completion
      
      info <- precompute(2,2,rowwidth,colwidth,rowsep,colsep,cat.level=0)
      info <- oneByoneHack(info)
    } else {
      preCompTime <- proc.time()
      
      ## April 15th 2004
      ## this code doesn't work when nrows=1 and ncols=1
      
      ## Jan 21st changes
      ## somethings I compute are integers, but R stores them as real numbers
      ## accounting for this in the initial steps at least cut back on memory allocation
      ##
      
      n <- nrows*ncols
      if(length(rowsep)==1) rowsep <- rep(rowsep,nrows-1)
      if(length(colsep)==1) colsep <- rep(colsep,ncols-1)
      nrow <- rep(seq(1:nrows), ncols) - 1
      ncol <- rep(seq(1:ncols), rep(nrows, ncols)) - 1
      t1 <- proc.time()
      num.ifail <- 0

      ## vector of x-coordinates:
      ax.vals <- 0
      for(gen.i in 1:(nrows-1))
        ax.vals <- c(ax.vals,ax.vals[gen.i]+rowwidth+rowsep[gen.i])
      if(cat.level>2) cat(ax.vals,"\n")

      bx.vals <- 0
      for(gen.i in 1:(ncols-1))
        bx.vals <- c(bx.vals,bx.vals[gen.i]+colwidth+colsep[gen.i])
      if(cat.level>2) cat(bx.vals,"\n")

      indices <- cbind(gl(n,n,n^2),gl(n,1,n^2))
      ## only use pairs where second is to the right of the first
      indices <- indices[indices[,2]>indices[,1],]
      indices <- rbind(as.integer(c(1,1)),indices)
      if(cat.level>1) {
        z <- sapply(ls(), function(x) object.size(get(x)))
        cat("Memory Allocation is approx:",sum(z)/1000000,"MB\n")
      }

      ## following what I used to do
      getV.i.bl <- as.integer(indices[,1]%%nrows)
      getV.i.bl[getV.i.bl==0] <- as.integer(nrows)
      getV.j.bl <- as.integer((indices[,1]-getV.i.bl)/nrows+1)
      getV.i.nt <- as.integer(indices[,2]%%nrows)
      getV.i.nt[getV.i.nt==0] <- as.integer(nrows)
      getV.j.nt <- as.integer((indices[,2]-getV.i.nt)/nrows+1)
      getV.i <- as.integer(abs(getV.i.nt-getV.i.bl)+1)
      getV.j <- as.integer(abs(getV.j.nt-getV.j.bl)+1)

      ax <- abs(ax.vals[getV.i.bl]-ax.vals[getV.i.nt])
      bx <- abs(bx.vals[getV.j.bl]-bx.vals[getV.j.nt])
      ax <- round(ax,12)
      bx <- round(bx,12)
      rm(getV.i.bl, getV.j.bl)
      rm(getV.i.nt, getV.j.nt)
      if(cat.level>1) {
        z <- sapply(ls(), function(x) object.size(get(x)))
        cat("Memory Allocation is approx:",sum(z)/1000000,"MB\n")
      }

      ## perhaps I can improve this step
      ## firstly getV.i, getV.j and indices are all integers
      ## but ax and bx are not, is it possible to store
      ## these in two groups to save memory
      ## if so, how do I use the apply command later
      ## later I will tack on another 5 columns,
      ## 3 of type integer and 2 of real
      ## now use lists
      indices <- list(ax=ax,bx=bx,getV.i=getV.i,getV.j=getV.j,indices=indices)
      rm(ax,bx,getV.i,getV.j)
      if(cat.level>1) {
        z <- sapply(ls(), function(x) object.size(get(x)))
        cat("Memory Allocation is approx:",sum(z)/1000000,"MB\n")
      }

      ## which of the rows need evaluation
      ## does not take row column symmetry into account
      evalFactor <- as.integer(as.factor(indices[[1]]):as.factor(indices[[2]]))

      ## now sort according to the numerical value of evalFactor
      indices$evalFactor <- evalFactor
      ord <- order(evalFactor)

      for(count.i in 1:length(indices))
        {
          if(is.matrix(indices[[count.i]])) {
            indices[[count.i]] <- indices[[count.i]][ord,]
          } else {
            indices[[count.i]] <- indices[[count.i]][ord]
          }
        }
      
      ## now make an old copy, and in the current one only keep the distinct rows
      ## all that is needed from the current version is the list of locations
      locations <- indices$indices
      if(cat.level>1) {
        z <- sapply(ls(), function(x) object.size(get(x)))
        cat("Memory Allocation is approx:",sum(z)/1000000,"MB\n")
      }

      keepRows <- (c(1,diff(indices$evalFactor))==1)
      keepRows <- (1:length(indices[[1]]))[keepRows]
      for(count.i in 1:length(indices))
        {
          if(is.matrix(indices[[count.i]])) {
            indices[[count.i]] <- indices[[count.i]][keepRows,]
          } else {
            indices[[count.i]] <- indices[[count.i]][keepRows]
          }
        }
      if(cat.level>1) {
        z <- sapply(ls(), function(x) object.size(get(x)))
        cat("Memory Allocation is approx:",sum(z)/1000000,"MB\n")
      }

      indices.preLimits <- indices

      ## now convert indices into a matrix again:
      indices <- cbind(indices[[1]],indices[[2]],indices[[3]],indices[[4]],indices[[5]],indices[[6]])

      ## limits of integration, may or may not be used
      ## how to incorporate aniso here...
      c.pts <- cbind((indices[,1]-rowwidth)^2+(indices[,2]-colwidth)^2, indices[,1]^2+(indices[,2]-colwidth)^2, indices[,1]^2+indices[,2]^2, (indices[,1]-rowwidth)^2+indices[,2]^2, indices[,1]^2+(indices[,2]-colwidth)^2, (indices[,1]+rowwidth)^2+(indices[,2]-colwidth)^2, (indices[,1]+rowwidth)^2+indices[,2]^2, indices[,1]^2+indices[,2]^2,  indices[,1]^2+indices[,2]^2, (indices[,1]+rowwidth)^2+indices[,2]^2, (indices[,1]+rowwidth)^2+(indices[,2]+colwidth)^2, indices[,1]^2+(indices[,2]+colwidth)^2, (indices[,1]-rowwidth)^2+indices[,2]^2, indices[,1]^2+indices[,2]^2, indices[,1]^2+(indices[,2]+colwidth)^2, (indices[,1]-rowwidth)^2+(indices[,2]+colwidth)^2)
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

      ## function used to collapse the partial integrals into the corrent number of results
      collapse.sum <- function(vect,positions)
        {
          res <- NULL
          pos <- 1
          for(count.i in 1:length(positions))
            {
              res <- c(res,sum(vect[pos:(pos+positions[count.i]-2)]))
              pos <- pos+positions[count.i]-1
            }
          res
        }

      rowReps <- as.integer(summary.factor(as.factor(evalFactor),maxsum=length(evalFactor)))
      ## this indicates how many times a particular entry will have to be recorded in the matrix
      ## it has nothing to do with the rows in the lattice, more to do with the number of repetitions in the rows of the matrix that contains the results, unfortunate variable name

      ## what is needed from here:
      ## geometric info associated with the results below
      ## indicesNEW - information on what needs to be computed
      ## collapse.sum - converts the computed into the correct results
      ##              - by adding the parts from the different integrals
      ##              - or in the case analytic results, the result and some zeros
      ## rowReps      - needed to expand out the results
      ## locations
      ## lengths - know how much to expand the results by

      ## relabel these:
      ## info
      ## collapse.sum
      ## locations
      ## lengths as part of info, to be used in collapse.sum

      info <- list(rowwidth=rowwidth,colwidth=colwidth,rowsep=rowsep,colsep=colsep,nrows=nrows,ncols=ncols,lengths=lengths,rowReps=rowReps,locations=locations,indices=indicesNEW,aniso=1,indices.preLimits=indices.preLimits)
      
      if(cat.level) cat("precompute takes",(proc.time()-preCompTime)[3],"seconds\n")
    }
    info
  }
