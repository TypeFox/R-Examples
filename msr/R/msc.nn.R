#type 1 indiacte to use persistence for merging, type 2 uses
#R^2 based merging criterion
msc.nn <- function (y, x, knn = ncol(x), pLevelP = 0.2, pLevel, 
                    nLevels, type = 1, smooth = FALSE, eps=0.01) 
{
    this.call <- match.call()
  
    if(is.null(nrow(x))){
      x <- as.matrix(x, ncol=1)
    }
    else{ 
      x <- as.matrix(x)
    }
    if(missing(pLevel)){
      if( pLevelP > 0 && pLevelP < 1 ){
        pLevel <- ( max(y) - min(y) )*pLevelP
      }
      else{
        pLevel <- 0
      }
    }

    #add minimal amount of noise to avoid duplicated values and failure of the MS computation
    if(anyDuplicated(y)){
      tmp <- max(y) - min(y)
      y <- y + runif(length(y))*tmp*0.000001
    }

    l <- length(y)
    nr <- nrow(x)
    nc <- ncol(x)
    if( nr != l ){
      stop("x and y don't have the same number of observations") 
    }

   
    #call C morse-smale function
    mscl <- c()
    pers <- NULL
    if(missing(nLevels) || is.null(nLevels)){
       nLevels <- 1
       res <- NULL
       if(type == 1){
         res <- .Call("nnmspartition", as.integer(nc), as.integer(nr), as.double(y),
                       as.double(t(x)), as.integer(knn), as.double(pLevel), 
					   as.logical(smooth), as.double(eps))
       }
       else{
         res <- .Call("nnmscR2", as.integer(nc), as.integer(nr), as.double(y),
                      as.double(t(x)), as.integer(knn), as.logical(smooth), as.double(eps))
       }
       pers <- as.vector(res[[4]]) 
       mscl[[1]] <- structure(list(call=call, pLevel=-1,
            partition = as.vector(res[[1]]), maxs=as.vector(res[[2]]), 
            mins=as.vector(res[[3]])), class = "msc.level")
        mscl[[1]]$partitionSize <- c()
        for(k in 1:length(mscl[[1]]$maxs)){
          mscl[[1]]$partitionSize[[k]] = sum(mscl[[1]]$partition == k) 
        } 
    }
    else{
      if(type == 1){
        res <- .Call("nnmsc2", as.integer(nc), as.integer(nr), as.double(y),
                      as.double(t(x)), as.integer(knn), as.integer(nLevels), 
					  as.logical(smooth), as.double(eps))
      }
      else{
        res <- .Call("nnmsc2R2", as.integer(nc), as.integer(nr), as.double(y),
                      as.double(t(x)), as.integer(knn), as.integer(nLevels), 
					  as.logical(smooth), as.double(eps))
      }

      pers <- as.vector(res[[1]])
      pl <- length(pers)

      if(nLevels > pl){
        nLevels <- pl
      }
      for(i in 0:(nLevels-1)){
        p <- pers[pl-i]
        msc <- structure(list(call=call, pLevel=p, 
                  partition = as.vector(res[[i*3+2]]), maxs = as.vector(res[[i*3+3]]), 
                  mins=as.vector(res[[i*3+4]])), class = "msc.level")
        msc$partitionSize <- c()
        for(k in 1:length(msc$maxs)){
          msc$partitionSize[[k]] = sum(msc$partition == k) 
        } 
        mscl[[i+1]] <- msc
      }

    }
    
    obj <- structure(list(call=call, y=y, x=x, knn = knn, level=mscl, 
            persistence=pers, nLevels=nLevels, predictLevel = nLevels), class="msc") 
    obj
}




#Extract a subset of the persistence hierarchy of ms
msc.sublevels <- function(ms, startLevel = ms$predictLevel, endLevel = startLevel){
  ms$level <- ms$level[startLevel:endLevel]
  ms$predictLevel <- endLevel-startLevel+1
  ms$nLevels <- endLevel - startLevel + 1
  ms
}




