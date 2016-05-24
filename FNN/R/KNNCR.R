################################################################################
# KNN Classification and Regression                                            #
# File:   KNNCR.R                                                              #
# Author: Shengqiao Li                                                         #
# Date:   December 12, 2008                                                    #
#                                                                              #
# for small dataset, nrow(train)<50 and nrow(test)<50,                         #
# ANN is slower than VR's method.                                              #
# for large dataset, ANN is faster                                             #
################################################################################
# note: data input for VR_knnc, VR_knnr in n x d format.
#
knn<- function(train, test, cl, k=1, prob=FALSE, 
              algorithm=c("kd_tree", "cover_tree", "brute")
              )
{
  algorithm<- match.arg(algorithm);
  
  train <- as.matrix(train)
  if (is.null(dim(test))) 
      dim(test) <- c(1, length(test))
  test <- as.matrix(test)
  if (any(is.na(train)) || any(is.na(test)) || any(is.na(cl))) 
      stop("no missing values are allowed")
  p <- ncol(train)
  ntr <- nrow(train)
  if (length(cl) != ntr) 
      stop("'train' and 'class' have different lengths")
  if (ntr < k) {
      warning(gettextf("k = %d exceeds number %d of patterns", k, ntr), domain = NA)
      k <- ntr
  }
  if (k < 1)
      stop(gettextf("k = %d must be at least 1", k), domain = NA)
  nte <- nrow(test)
  if (ncol(test) != p) 
      stop("dims of 'test' and 'train' differ")
  clf <- as.factor(cl)
  nc <- max(unclass(clf))
  switch(algorithm,
         cover_tree =,
         kd_tree =,
         brute = {Z<- get.knnx(data=train, query=test, k=k, algorithm);
                      nn.class<- matrix(cl[Z$nn.index], ncol=k); #factor levels are taken.      
                      pred_prob<- function(x) 
                      {
                        freq<- table(x)
                        prob<- freq/k;
                        max.ind<- which.max(freq)         
                        list(names(max.ind), prob[max.ind])
                      }      
                      res<- apply(nn.class, 1, pred_prob);      
                      res<- matrix(unlist(res), ncol=2, byrow=TRUE);
                      if(prob) pr<- as.numeric(res[,2]);
                      res<- as.factor(res[,1]);
         }
      )
      
  if(prob) attr(res, "prob") <- pr;
  attr(res, "nn.index")<- matrix(Z$nn.index, ncol=k);;
  attr(res, "nn.dist")<- matrix(Z$nn.dist, ncol=k);
 
  return(res);
}

knn.cv <- function(train, cl, k=1, prob=FALSE, 
          algorithm=c("kd_tree", "cover_tree", "brute")
          )
{
  algorithm<- match.arg(algorithm);
  
  train <- as.matrix(train)
  if(any(is.na(train)) || any(is.na(cl)))
     stop("no missing values are allowed")
  p <- ncol(train)
  ntr <- nrow(train)
  if(length(cl) != ntr) stop("'train' and 'class' have different lengths")
  if(ntr-1 < k) {
    warning(gettextf("k = %d exceeds number %d of patterns", k, ntr-1), domain = NA)
      k <- ntr - 1
    }
  if (k < 1)
            stop(gettextf("k = %d must be at least 1", k), domain = NA)
  clf <- as.factor(cl)
  nc <- max(unclass(clf))

  switch(algorithm,
     cover_tree =,
     kd_tree  =,
     brute = {Z<- get.knn(data=train, k=k, algorithm=algorithm);
                  nn.class<- matrix(cl[Z$nn.index], ncol=k); #factor levels are taken.      
                  pred_prob<- function(x) 
                  {
                    freq<- table(x)
                    prob<- freq/k;
                    max.ind<- which.max(freq)         
                    list(names(max.ind), prob[max.ind])
                  }
  
                  res<- apply(nn.class, 1, pred_prob);      
                  res<- matrix(unlist(res), ncol=2, byrow=TRUE);
                  if(prob) pr<- as.numeric(res[,2]);
                  res<- as.factor(res[,1]);
     }
    
  );
      
  if(prob) attr(res, "prob") <- pr;
  attr(res, "nn.index")<- matrix(Z$nn.index, ncol=k);;
  attr(res, "nn.dist")<- matrix(Z$nn.dist, ncol=k);

  return(res);
}

knn.reg<- function(train, test=NULL, y, k=3, algorithm=c("kd_tree", "cover_tree", "brute")
          )
{
  #KNN regression. If test is not supplied, LOOCV is performed and R2 is the predicted R2
  algorithm<- match.arg(algorithm); 
  
  train<- as.matrix(train);
  if(!is.null(test)){
    if (is.null(dim(test)))  dim(test) <- c(1, length(test))   #1 x p 
    test <- as.matrix(test)
  }
  ntr<- nrow(train); p<- ncol(train);
  n<- ifelse(is.null(test), nrow(train), nrow(test)); #number of samples to be predict
 
  pred<- switch(algorithm,
            cover_tree =, 
            kd_tree =,
            brute = {Z<- if(is.null(test)) get.knn(train, k, algorithm)
                           else get.knnx(train, test, k, algorithm);
                      rowMeans(matrix(y[Z$nn.index], ncol=k));
            }          
  )

  if(is.null(test)){
    residuals<- y-pred ;
    PRESS<- sum(residuals^2);
    R2<- 1-PRESS/sum((y-mean(y))^2);
  }
  else{
    residuals<-  PRESS<-  R2<- NULL;
  }
  
  res <- list(call = match.call(), k = k, n = n, pred = pred, 
      residuals = residuals, PRESS=PRESS, R2Pred = R2)

  class(res)<- if(!is.null(test)) "knnReg" else "knnRegCV";
  
  return(res);
}
print.knnRegCV <-function (x, ...) 
{      
    if (!inherits(x, "knnRegCV")) 
        stop(deparse(substitute(x)), " is not a knnRegCV object");
        
    cat("PRESS = ", x$PRESS, "\n")
    cat("R2-Predict = ", x$R2Pred, "\n")
}
print.knnReg <-function (x, ...)
{
    if (!inherits(x, "knnReg"))
        stop(deparse(substitute(x)), " is not a knnReg object");

    cat("Prediction:\n");
    
    print(x$pred);
}



