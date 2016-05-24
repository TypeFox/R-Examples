msc.slm <- function (ms, nfold = 10, modelSelect = FALSE) 
{
    obj <- NULL
    slm <- c()
    pLevel = 1
    if(nfold > 1){
      nLevels <- ms$nLevels
      minCV <- Inf
      minobj <- NULL
      nr <- nrow(ms$x)
      nc <- ncol(ms$x) +1
      split <- ceiling(nr/nfold)
      df <- data.frame(y = ms$y, ms$x)
      for(i in 1:nLevels){
        ms$predictLevel <- i
        obj <- structure(list(ms = ms), class = "msc.slm")
        cv = 0;
        mmdf <- data.frame( y = df[,1], model.matrix(obj, df[, 2:nc]) )
        for(n in 1:nfold){
          s <- sample(1:nr)
          va <- s[1:split]
          tr <- s[(split+1):nr]
          lm <- lm(y ~ 0+., mmdf[tr,])
          if(modelSelect){
            lm <- step(lm, trace = 0)
          }
          yp = predict(lm, mmdf[va, ])
          cv = cv + mean((yp - df[va, 1])^2)
        }
        cv = cv/nfold
        if(cv< minCV){
          pLevel = i 
          minCV <- cv
        }        
        mm <- model.matrix(obj, ms$x)
        mmdf <- data.frame(y = ms$y, mm)
        lm <- lm(y ~ 0+., mmdf)
        if(modelSelect){
          lm <- step(lm, trace = 0)
        }
        slm[[i]] <- structure(list(lm = lm, cv = cv), name="msc.slm.level")
      }
    } 
    else{
      nLevels <- ms$nLevels
      minBIC <- Inf
      minobj <- NULL
      for(i in 1:nLevels){
        ms$predictLevel <- i
        obj <- structure(list(ms = ms), class = "msc.slm")
        mm <- model.matrix(obj, ms$x)
        df <- data.frame(y = ms$y, mm)
        slm <- lm(y ~ 0+., df)
        bic <- BIC(slm)
        slm[[i]] <- structure(list(lm = lm, bic = bic), name="msc.slm.level")
        if(bic == -Inf){
          break;
        } 
        if(bic < minBIC){
          pLevel = i;
        }
      }
    }
    ms$predictLevel <- pLevel 
    obj <- structure(list(ms = ms, slm = slm), class = "msc.slm")
    obj
}

