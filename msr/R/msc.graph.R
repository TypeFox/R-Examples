msc.graph <- function (y, x, knn, knnd, nLevels, smooth = FALSE) 
{
    this.call <- match.call()
  
    if(is.null(nrow(x))){
      x <- as.matrix(x, ncol=1)
    }
    else{ 
      x <- as.matrix(x)
    }


    #add minimal amount of noise to avoid duplicated values and failure of the MS computation
    if(anyDuplicated(y)){
      tmp <- max(y) - min(y)
      y <- y + runif(length(y))*tmp*0.000001
    }

    l <- length(y)
    nr <- nrow(x)
    nc <- ncol(x)
    k <- ncol(knn)
    
    if( nr != l ){
      stop("x and y don't have the same number of observations") 
    }

   
    #call C morse-smale function
    mscl <- c()
    pers <- NULL
    res <- .Call("graphmsc", as.integer(nc), as.integer(nr), as.integer(k), 
                             as.double(y), as.double(t(x)), as.integer(t(knn-1)), 
                             as.double(t(knnd)), as.integer(nLevels), 
			     as.logical(smooth))

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

    
    obj <- structure(list(call=call, y=y, x=x, knn = knn, level=mscl, 
            persistence=pers, nLevels=nLevels, predictLevel = nLevels), class="msc") 
    obj
}









