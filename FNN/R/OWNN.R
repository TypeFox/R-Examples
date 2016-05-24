#optimal weighted nn
ownn<- function(train, test, cl, testcl=NULL, k=NULL, prob=FALSE, algorithm=c("kd_tree", "cover_tree", "brute"))
{
  algorithm <- match.arg(algorithm)
  train <- as.matrix(train);
  if (is.null(dim(test)))
      dim(test) <- c(1, length(test))
      
  test <- as.matrix(test);
  if (any(is.na(train)) || any(is.na(test)) || any(is.na(cl)))
      stop("no missing values are allowed")
      
  d<- p<- ncol(train);
  
  n <- nrow(train);
  if (length(cl) != n)
      stop("'train' and 'class' have different lengths")
  nte <- nrow(test);
  if (ncol(test) != p)
      stop("dims of 'test' and 'train' differ")
      
  clf <- as.factor(cl);

  #Z<- knnx.index(data=train, query=test, k=n, algorithm='kd_tree');
  Z<- get.knnx(data=train, query=test, k=n, algorithm=algorithm)$nn.index;
  nn.class<- matrix(cl[Z], ncol=n); #factor levels are taken.
  
  pred_prob<- function(x, k, method=c("knn", "ownn", "bnn"))
  {
    method<- match.arg(method);
    x<- as.factor(x);
    cl.lab<- levels(x);

    nc <- max(unclass(x));
  
    weight<- switch(method,
                    knn = c(rep(1/k,k),rep(0,n-k)),
                    ownn = {
                        kstar <- floor((2*(d+4)/(d+2))^(d/(d+4))*k)
                        i<- 1:kstar
                        alpha <- i^(1+2/d) - (i-1)^(1+2/d)
                        c((1+d/2 - d/(2*kstar^(2/d))*alpha)/kstar,rep(0,n-kstar))
                        },
                    bnn= {
                        #q <- gamma(2+2/d)^(2*d/(d+4))/(2^(d/(d+4))*k)
                        q <- 2^(d/(d+4))*gamma(2+2/d)^(2*d/(d+4))*(1/k);
                        q*(1-q)^(1:n-1)/(1-(1-q)^n);
                        }
    );
    
    etaknn <- rep(0,nc);
    names(etaknn)<- cl.lab;
    for(r in 1:nc) etaknn[r] <- sum(weight[x==cl.lab[r]])

    prob<- etaknn/sum(etaknn);
    
    max.ind<- which.max(etaknn);
    list(names(max.ind), prob[max.ind]);
  }

  knn.predict<- function(nn.class, k, method)
  {
  
    res<- apply(nn.class, 1, pred_prob, k=k, method=method);
    res<- matrix(unlist(res), ncol=2, byrow=TRUE);
    if(prob) pr<- as.numeric(res[,2]);
    res<- as.factor(res[,1]);

    if(prob) attr(res, "prob") <- pr;

    return(res);
  }
  
  Choosek <- function(Y, d, nn.class, klower, kupper,gridpoints=20)
  {
    n <- length(Y)
    Y <- as.factor(Y)
    fold <- 0
    while(sum(fold==1) == 0 || sum(fold==2) == 0 || sum(fold==3) == 0 || sum(fold==4) == 0 || sum(fold==5) == 0)fold <- sample(1:5,n,replace=T)

    CVPred2 <- matrix(NA, nrow=n, ncol=gridpoints+1)

    k <- rep(0, gridpoints+1)

    for(g in 0:gridpoints){
      k[g+1] <- round((gridpoints-g)*klower/gridpoints + g*kupper/gridpoints)
      for(i in 1:5){
        CVPred2[fold==i,g+1] <- knn.predict(nn.class[fold==i, fold!=i], k=k[g+1], method="knn")
      }
    }
    CVRisk <- apply(CVPred2,2,function(x) sum(Y!=x))

    return(floor(k[which.min(CVRisk)]*(5/4)^(4/(d+4))))
  }

  if(is.null(k)) {
    Z_train<- get.knnx(data=train, query=train, k=n, algorithm=algorithm)$nn.index;
    nn.class_train<- matrix(cl[Z_train], ncol=n); #factor levels are taken.
    
    k<- Choosek(Y=clf, d=d, nn.class=nn.class_train, klower=5, kupper=n/2, gridpoints=20)
  }
  
  knnpred <- knn.predict(nn.class,  k, method="knn");
  ownnpred <- knn.predict(nn.class, k, method="ownn");
  bnnpred  <- knn.predict(nn.class, k, method="bnn");
  
  if(!is.null(testcl))
  {
    knnacc <- mean(testcl==knnpred);
    ownnacc <- mean(testcl==ownnpred);
    bnnacc <- mean(testcl==bnnpred);
    return(list(k=k, knnpred=knnpred, ownnpred=ownnpred, bnnpred=bnnpred, accuracy=c(knn=knnacc, ownn=ownnacc, bnn=bnnacc)))
 }
 else return(list(k=k, knnpred=knnpred, ownnpred=ownnpred, bnnpred=bnnpred))
}

