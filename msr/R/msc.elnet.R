msc.elnet <- function (ms, nfold = 10, blend=FALSE){

  buildelnet <- function(ms){
    msLevel <- ms$level[[ms$predictLevel]]
    nc <- ncol(ms$x)
    nr <- nrow(ms$x)
    

    #fit model to each crystal 
    lm <- vector("list", length(ms$mins))
    df <- data.frame(y=ms$y, ms$x) 
    cv <- 0;
    for(i in 1:length(msLevel$mins)){
      index = msc.level.ind(msLevel, i);
      x <- ms$x[index,]
      y <- ms$y[index]    
      cvm <- NULL
      cvsd <- NULL
      lambda <- NULL
      if(nfold > 0){
        cvg <- cv.glmnet(y = y, x=x, nfolds=nfold, standardize=FALSE)
        lambda <- cvg$lambda.min
        index <- which.min(cvg$cvm) 
        cvm <- cvg$cvm[index]
        cvsd <- cvg$cvsd[index]
        cv <- cv + nrow(x)/nr * cvm
      }

      df <- data.frame(y, x) 
      lm[[i]] = glmnet(y = y, x = x, lambda=lambda, standardize = FALSE)
    }

      obj <- structure(list(lm = lm, cv= cv), class =  "msc.elnet.level")
      obj
    }

    elnet <- c();
    minCV <- Inf
    for(i in 1:ms$nLevels){
      ms$predictLevel <- i
      elnet[[i]] <- buildelnet(ms)
      if(minCV > elnet[[i]]$cv){
        minCV <- elnet[[i]]$cv
        pLevel <- i
      }
    }
    ms$predictLevel <- pLevel 
    obj <- structure(list(ms=ms, elnet = elnet, blend=blend), class = "msc.elnet")
    obj
}

