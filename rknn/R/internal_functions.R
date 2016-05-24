##########################################################################################
# Random KNN Internal Functions                                                
# These internal functions are not intent for end users                        
# File:   internal_functions.R                                                 
# Author: Shengqiao Li                                                         
# Date:   June 24, 2008 (initial)                                              
#         2013-08-05 replace choose.bigz by gmp::chooseZ, factorial.bigz by gmp::factorialZ
###########################################################################################

#factorial.bigz<- function(n)
#{
#    if(n==0) return(1);
#    prod.bigz(as.bigz(1:n))
#}
#
#chooses.bigz<- function(n, m)
#  {
#    if(m>n) return(0)
#    else if(m==0||m==n) return(1)
#
#    res<- prod.bigz(as.bigz(n:(n-m+1)))/prod.bigz(as.bigz(1:m));
#
#    as.bigz(res);
#}

#replaced by gmp::chooseZ
#choose.bigz<- function(n, m)
#{
#    ml<- length(n);
#    nl<- length(m);
#    vl<- max(ml, nl);

#    #recycle the short one
#    if(ml<vl) n<- rep(n, length=vl)
#    else if(nl<vl) m<- rep(m, length=vl)
#
#    res<- as.bigz(vl)
#
#    for(i in 1:vl) res[i]<- chooses.bigz(n[i], m[i])
#
#    return(res)
#}


splitR<- function(r, nc)
{
  r.min<- r%/%nc;
  remain<- r-r.min*nc;
  rep(c(r.min, r.min+1), times=c(nc-remain, remain))      
}


#
# disable set.return.seed since the specific "L'Ecuyer-CMRG" RNG is required for a cluster 
#
#set.return.seed<- function(Random.seed=NULL, seed=NULL)
#{
#  #set and record random seeds even no seed is supplied
#  #set.seed only accepts single integer seed and don't return seed
#  #Random.seed  -- a seed in the .Random.seed format
#  #seed  -- an integer seed
#                                            
# if(!is.null(Random.seed)) {
#    assign(".Random.seed", Random.seed, envir=.GlobalEnv)
#  }
#  else  if(!is.null(seed))  set.seed(seed)
#  else if(!exists(".Random.seed", envir=.GlobalEnv)) runif(1);
#
#  Random.seed<- get(".Random.seed", envir=.GlobalEnv);
#         
#  invisible(Random.seed);
#}


knn<- function (train, test, cl, k = 1) 
{

    train <- as.matrix(train)
    
    if (is.null(dim(test))) dim(test) <- c(1, length(test))
    test <- as.matrix(test)
    
    if (any(is.na(train)) || any(is.na(test)) || any(is.na(cl))) 
        stop("no missing values are allowed")
        
    p <- ncol(train)
    ntr <- nrow(train)
    
    if (length(cl) != ntr) 
        stop("'train' and 'class' have different lengths")
    if (ntr < k) {
        warning(gettextf("k = %d exceeds number %d of patterns", 
            k, ntr), domain = NA)
        k <- ntr
    }
    if (k < 1) stop(gettextf("k = %d must be at least 1", k), domain = NA)
    
    nte <- nrow(test)
    if (ncol(test) != p) stop("dims of 'test' and 'train' differ")
        
    clf <- as.factor(cl)
    nc <- max(unclass(clf))
  
    use.all = TRUE;
   
    Z <- .C(C_knnc, as.integer(k),  as.integer(ntr), 
        as.integer(nte), as.integer(p), as.double(train), 
        as.integer(unclass(clf)), as.double(test), res = integer(nte), 
        pr = double(nte), integer(nc + 1), as.integer(nc), 
        as.integer(FALSE), as.integer(use.all), nn.index = integer(nte * k), 
        nn.dist = double(nte * k)) 
    
    res <- factor(Z$res, levels = seq_along(levels(clf)), labels = levels(clf))

    return(res)
}
       
knn.cv<- function (train, cl, k = 1) 
{

    train <- as.matrix(train)
    if (any(is.na(train)) || any(is.na(cl))) 
        stop("no missing values are allowed")
    p <- ncol(train)
    ntr <- nrow(train)
    if (length(cl) != ntr) 
        stop("'train' and 'class' have different lengths")
    if (ntr - 1 < k) {
        warning(gettextf("k = %d exceeds number %d of patterns", 
            k, ntr - 1), domain = NA)
        k <- ntr - 1
    }
    if (k < 1) 
        stop(gettextf("k = %d must be at least 1", k), domain = NA)
        
    clf <- as.factor(cl)
    nc <- max(unclass(clf))
    
   use.all = TRUE;

    Z <- .C(C_knnc, as.integer(k),  as.integer(ntr), 
        as.integer(ntr), as.integer(p), as.double(train), 
        as.integer(unclass(clf)), as.double(train), res = integer(ntr), 
        pr = double(ntr), integer(nc + 1), as.integer(nc), 
        as.integer(TRUE), as.integer(use.all), nn.index = integer(ntr * 
            k), nn.dist = double(ntr * k))
    res <- factor(Z$res, levels = seq_along(levels(clf)), labels = levels(clf))

    return(res)
}

#KNN regression. If test is not supplied, LOOCV is performed and R2 is the predicted R2 
knn.reg<- function (train, test = NULL, y, k = 3) 
{
    train <- as.matrix(train)
    if (!is.null(test)) {
        if (is.null(dim(test))) 
            dim(test) <- c(1, length(test))
        test <- as.matrix(test)
    }
    ntr <- nrow(train)
    p <- ncol(train)
    n <- ifelse(is.null(test), nrow(train), nrow(test))
    
    use.all<-  FALSE;
    pred <- .C(C_knnr, as.integer(k), as.integer(ntr), as.integer(n), 
        as.integer(p), as.double(train), as.double(y), as.double(if (is.null(test)) train else test), 
        res = double(n), as.integer(is.null(test)), as.integer(use.all), 
        nn.index = integer(n * k), nn.dist = double(n * k))$res;
    if (is.null(test)) {
        residuals <- y - pred
        PRESS <- sum(residuals^2)
        R2 <- 1 - PRESS/sum((y - mean(y))^2)
    }
    else {
        residuals <- PRESS <- R2 <- NULL
    }
    res <- list(call = match.call(), k = k, n = n, pred = pred, 
        residuals = residuals, PRESS = PRESS, R2Pred = R2)
    class(res)<- if(!is.null(test)) "knnReg" else "knnRegCV";

    return(res)
}
  
#########################################################################
