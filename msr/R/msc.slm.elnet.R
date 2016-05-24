msc.slm.elnet <- function (ms, nfold = 10) {

  #build the model
  buildmodel <- function(ms, nfold){
    obj <- structure(list(ms = ms), class = "msc.slm.elnet" )
    mm <- model.matrix(obj, ms$x)

    cvg <- cv.glmnet(y = ms$y, x=mm, nfolds=nfold, standardize=FALSE)
    index <- which.min(cvg$cvm) 
    cvm <- cvg$cvm[index]
    cvsd <- cvg$cvsd[index]
    obj <- structure(list(ms = ms, elnet = cvg, nfold = nfold,
          cv = cvm, cvsd =cvsd), class = "msc.slm.elnet.level")
    obj
  }
  #end build model function

  nLevels <- ms$nLevels
  mincvm <- Inf
  pLevel =1
  slm <-c()
  for(i in 1:nLevels){
    ms$predictLevel <- i
    slm[[i]] <- buildmodel(ms, nfold=nfold)
        
    if(slm[[i]]$cv < mincvm){
     mincvm <- slm[[i]]$cv
     pLevel=i
    }
  }
  ms$predictLevel = pLevel
  obj <- structure(list(ms = ms, slm = slm), class = "msc.slm.elnet")
  obj
}

