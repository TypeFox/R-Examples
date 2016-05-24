lmcv <-
function(obj, ng = NULL, random = FALSE) {
  #
  # The Leave-One-Out (LOO) and/or Leave-Group-Out (LGO) Cross-Validation in R for (Multiple) Linear Regression.
  #
  # Input:
  #   obj: the model of MLR
  #   ng: number of group, if missing, do LOO
  #   random: logical, if TRUE, do random CV
  #
  # Output:
  #   q.squared: cross-validation relation coefficient.
  #   SDEP: Standard Deviation of Error of Predictions
  #   newsq: variance in Y explained only for LOO CV
  #
  # Usage:
  #   loo <- lmcv(obj)
  #   lgo <- lmcv(obj, ng = 5)
  data<-data.frame(obj$model)
  col.names <- colnames(data)
  colnames(data)[1]<-c("expr")
  N <- nrow(data)
  if (random == TRUE) data <- data
  if (missing(ng)) ng <- N 
  ytest <- data[,1]
  ypred <- numeric(N)
  newrsq <- numeric(ng)
  g <- N %/% ng
  for (i in 1:ng) {
    if (g == 1) {
      index <- i 
      newtrain <- data[-index,]
      newlm <- lm(expr ~., data = newtrain)
      newrsq[i] <- summary(newlm)$r.squared
      ypred[index] <- predict(newlm,newdata=data[index,])
    }
    else {
      index <- c(i, ng * seq(1, (g - 1)) + i)
      if (N %% ng != 0 & i <= N %% ng) index <- c(index, (g * ng + i))
      newtrain <- data[-index,]
      newlm <- lm(expr ~., data = newtrain)
      newrsq[i] <- summary(newlm)$r.squared
      ypred[index] <- predict(newlm,newdata=data[index,])
    }
   
  }  
  ypred<-as.data.frame(ypred)  
  ytest<-as.data.frame(ytest)
  rownames(ypred)<-rownames(data)
  rownames(ytest)<-rownames(data)
  q.squared <- 1 - sum((ytest[,1] - ypred[,1])^2) / sum((ytest[,1] - mean(ytest[,1]))^2)
  SDEP <- sqrt(sum((ytest[,1] - ypred[,1])^2) / N)
  return(list(q.squared = q.squared, SDEP = SDEP))
  
}
