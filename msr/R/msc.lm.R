msc.lm <- function (ms, nfold = 10, modelSelect = FALSE, blend=FALSE,
                    verbose=FALSE) 
{

    
    #compute weights for each crystals    
#d <- predict(ms)
      
     #fit model to each crystals
      buildlm <- function(ms){
        msLevel <- ms$level[[ms$predictLevel]]
        nc <- ncol(ms$x)
        nr <- nrow(ms$x)
        lm <- vector("list", length(ms$mins))
        df <- data.frame(y=ms$y, ms$x) 
        r2 <- 0;
        cvAll <- 0 
        for(i in 1:length(msLevel$mins)){
          index = msc.level.ind(msLevel, i);
          x <- ms$x[index,]
          y <- ms$y[index]
          df <- data.frame(y, x) 
          lm[[i]] <- lm(y ~ ., data=df)
          if(modelSelect){
            lm[[i]] <- step(lm[[i]], trace = 0)
          }
          s <- summary(lm[[i]])
          r2 <- r2 + nrow(x)/nr * s$r.squared;
          cv <- cv.lm(y~., df, nfold)
          cvAll  <- cvAll + nrow(x)/nr * cv$meanSE
        }
     
        obj <- structure(list(cv = cvAll, r2=r2, lm = lm), class = "msc.lm.Level")
      }

      lms <- c()
      minCV <- Inf
      pLevel <- 0
      for(i in 1:ms$nLevels){
        ms$predictLevel <- i
        lms[[i]] <- buildlm(ms)
        if(verbose){
          print(paste("Level: ", i))
          print(paste("R^2: ", lms[[i]]$r2))
          print(paste("CV mean SE:", lms[[i]]$cv))
        }
        if(minCV > lms[[i]]$cv){
          minCV <- lms[[i]]$cv
          pLevel <- i
        }
      }        
      ms$predictLevel <- pLevel
      obj <- structure(list(lms=lms, ms = ms, blend=blend), class = "msc.lm")
      obj
}

