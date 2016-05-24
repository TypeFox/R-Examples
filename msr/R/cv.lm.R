cv.lm <- function(formula, df, nfolds=10) 
{
  nr <- nrow(df)
  meanSE <- 0;
  medianSE <- 0
  if(nfolds < nr){
    split <- ceiling(nr/nfolds)
    for(i in 1:nfolds){          
      s <- sample(1:nr)
      l1 <- lm(formula, df[s[split:nr], ])
      yp <- predict(l1, df[s[1:split], ])
      meanSE <- meanSE + mean((yp - df[s[1:split], 1])^2)
      medianSE <- medianSE + median((yp - df[s[1:split], 1])^2)
    }
  }
  else{
    nfolds = nr;
    for(i in 1:nr){
      if(i == 1){          
        l1 <- lm(formula, df[2:nr, ])
      }
      else if( i==nr){
        l1 <- lm(formula, df[1:(nr-1), ])
      }
      else{
        l1 <- lm(formula, df[c(1:(i-1), (i+1):nr), ])
      }
      yp <- predict(l1, df[i, ])
      meanSE <- meanSE + (yp - df[i, 1])^2
    }    
    medianSE <- meanSE
  }
  cv <- list(meanSE = meanSE/nfolds, medianSE = medianSE/nfolds)
  cv
}

