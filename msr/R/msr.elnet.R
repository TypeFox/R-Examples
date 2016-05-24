msr.elnet <- function(ms, nfolds = 10) 
{

    
      buildlm <- function(ms, nfolds){
        nc <- ncol(ms$x)
        nr <- nrow(ms$x)
        lms <- vector("list", length(ms$mins))
        df <- data.frame(y=ms$y, ms$x) 
        if(length(ms$crystalsSize)  == 1) {
          df <- data.frame(y = ms$y, ms$x) 
          f <- formula(y ~ 0+.) #glment adds an intercept
          mm <- model.matrix(f, df);
        }
        else{
          df <- data.frame(y = ms$y, ms$x, cID = as.factor(ms$crystals)) 
          f <- formula(y ~ 0+. + .:cID) #glment adds an intercept
          mm <- model.matrix(f, df);
        }
        obj <- cv.glmnet(y = ms$y, x = mm, nfolds = nfolds, standardize = FALSE)
        obj 
      }


      buildslm <- function(ms, nfolds){
        obj <- structure(list(ms = ms), class = "msr")
        mm <- model.matrix(obj, ms$x)
        slm <- cv.glmnet(y = ms$y, x = mm, nfolds = nfolds, standardize = FALSE)
        slm
      }

      if(is.null(ms$nLevels)){ 
        l <- buildlm(ms, nfolds)
        if(class(ms) == "msc.svm"){
          if( ! is.null(ms$svm) ){
            df <- data.frame(cry = as.factor(ms$crystals), obj$x)
            ms$svm <- svm(cry ~ ., df, probability=TRUE, cost=ms$cost, scale = FALSE)
          }
        }
        obj <- structure(list(ms = ms, lm = l, slm = buildslm(ms, nfolds)), class = "msr.elnet")
      }
      else{      
        lms <- c()
        minCV <- Inf
        pLevel <- 0
        for(i in 1:ms$nLevels){
          lms[[i]] <- buildlm(ms$mscl[[i]], nfolds)
          cv <- min(lms[[i]]$cvm)
          if(minCV > cv){
            minBIC <- cv
            pLevel <- i
          }
        }
        if(class(ms) == "msc.svm"){
          if( ! is.null(ms$mscl[[pLevel]]$svm) ){
            df <- data.frame(cry = as.factor(ms$mscl[[pLevel]]$crystals), obj$mscl[[pLevel]]$x)
            ms$mscl[[pLevel]]$svm <- svm(cry ~ ., df, probability=TRUE, cost=ms$cost, scale = FALSE)
          }
        }
        slm <- buildslm(ms$mscl[[pLevel]], nfolds)
        
        obj <- structure(list(ms = ms, lms = lms, slm= slm, predictLevel = pLevel), class = "msr.elnet")
      }
      obj
}

