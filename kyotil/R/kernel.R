getK=function(X,kernel,para=NULL,X2=NULL,C = NULL){    
    
    kernel <- substr(tolower(kernel[1]),1,1)
        
    if (is.null(C)) {
        # by default, we set C=TRUE for ibs and hamming kernel b/c the R implementation uses nested loop and can be slow
        if (kernel %in% c("i","h")) C=TRUE else C=FALSE
    }    
    if(C) return(getKernel(X=X,kernel=kernel,para = para,X2=X2))
    
    
    if(!any(kernel == c("l","e","i","h")) && is.null(para)) stop("kernel parameter is not set")
    if(any(kernel == c("i","h"))){
        if(is.null(para)) para <- 1
        para <- rep(para,length.out = ncol(X))
        if(any(para < 0)) stop("negative weights encountered for IBS or HAMMING kernel")
    }   
    
    if (kernel=="r" || kernel=="e") {
        if (!is.null(X2)) {
            aux = X[rep(1:nrow(X),nrow(X2)),,drop=FALSE] - X2[rep(1:nrow(X2),each=nrow(X)),,drop=FALSE]
            dist.mat = matrix(rowSums(aux^2), nrow=nrow(X))
#            aux=X2[rep(1:nrow(X2),nrow(X)),] - X[rep(1:nrow(X),each=nrow(X2)),]
#            dist.mat = matrix(rowSums(aux^2), nrow=nrow(X2))
        } else {
            dist.mat = as.matrix(dist(X))^2
        }
    }
    
    if (is.null(X2)) X2=X
    switch(kernel, 
        p=(tcrossprod(X,X2)+1)^para, # polynomial
        r=exp(-para*dist.mat), # rbf
        e=dist.mat, # Euclidean distance
        l=tcrossprod(X,X2), # linear
        i = R_ibs(X,X2,para), # IBS 
        h = R_hamming.sim(X,X2,para),
        stop(kernel %+% " kernel not supported")
    )
    
}


# Author: Krisztian Sebestyen <ksebestyen@gmail.com>
getKernel <- function(X,X2 = NULL,kernel = c("linear","euclidean","polynomial","rbf","ibs","hamming"),para = NULL){
    
    
    # note: the order of the kernels to be matched against must be exactly as below to match C
    kernel <- substr(tolower(kernel[1]),1,1)
    if(!any(kernel == c("l","e","i","h")) && is.null(para)) stop("kernel parameter is not set")
    .kernel <- kernel
    kernel <- match(kernel,c("l","e","p","r","i","h"),nomatch = 0)
    if(kernel == 0) stop("Invalid kernel selected.")
    kernel <- kernel - 1 # 0-based indexing in C
    
    if(is.null(X2))X2 <- X
    if(!is.matrix(X)) X <- as.matrix(X)
    d1 <- dim(X)
    if(!is.double(X)){
        X <- as.double(X)
        dim(X) <- c(d1[1] , d1[2])
    }
    if(!is.matrix(X2)) X2 <- as.matrix(X2)
    d2 <- dim(X2)
    if(!is.double(X2)){
        X2 <- as.double(X2)
        dim(X2) <- c(d2[1] , d2[2])
    }
    
    # check data integrity x \in{k,k+1,k+2} for IBS kernel  
    if(.kernel == "i"){ 
        U <- unique(c(X))
        if( !( (length(U) <= 3) || (setequal(abs(apply(expand.grid(U,U),1,diff)),c(0,1,2))) ) ) 
            stop("Values in 'X' should be in {0,1,2} or in {k,k+1,k+2} for some k in R")
        else if(!setequal(c(0,1,2),U))
            warning("Values in 'X' should be in {0,1,2} but are in {k,k+1,k+2} for some k in R")
        U2 <- unique(c(X2))
        if( !( (length(U2) <= 3) || (setequal(abs(apply(expand.grid(U2,U2),1,diff)),c(0,1,2))) ) ) 
            stop("Values in 'X2' should be in {0,1,2} or in {k,k+1,k+2} for some k in R")
        else if(!setequal(c(0,1,2),U2))
            warning("Values in 'X2' should be in {0,1,2} but are in {k,k+1,k+2} for some k in R")
    }   
    # para = weights for 'hamming' and 'ibs'
    if(!is.null(para)){  
        if(any(.kernel == c("i","h"))){
            para <- rep(para,length.out = min(c(d1[2],d2[2])))
            if(any(para < 0)) stop("negative weights encountered for IBS or HAMMING kernel")
        }
        if(!is.double(para)) para <- as.double(para)
    }
    
    K <- matrix(NaN,d1[1],d2[1])    
    .Call("Call_getKernel", X, X2, as.integer(kernel), para,K)
    K
}

# R's 'euclidean' dist()^2
edist2 <- function(X,X2 = NULL){
    if(is.null(X2)) X2 <- X
    
    if(!is.matrix(X)) X <- as.matrix(X)
    dx <- dim(X)
    if(!is.double(X)){
        X <- as.double(X)
        dim(X) <- c(dx[1] , dx[2])
    }
    if(!is.matrix(X2)) X2 <- as.matrix(X2)
    dy <- dim(X2)
    if(!is.double(X2)){ 
        X2 <- as.double(X2)
        dim(X2) <- c(dy[1] , dy[2])
    }
    
    K <- matrix(NaN,dx[1],dy[1])
    .Call("Call_edist2", X, X2,K)
    K
}


# if 'X' and 'X2' have differing column dimensions, the first 1..min(ncol(X),ncol(X2)) are used
ibs2 <- function(X,X2 = NULL,para = NULL){
    
    if(is.null(X2))X2 <- X
    
    # check data integrity x \in{k,k+1,k+2} for IBS kernel  
    U <- unique(c(X))
    if( !( (length(U) <= 3) || (setequal(abs(apply(expand.grid(U,U),1,diff)),c(0,1,2))) ) ) 
        stop("Values in 'X' should be in {0,1,2} or in {k,k+1,k+2} for some k in R")
    else if(!setequal(c(0,1,2),U))
        warning("Values in 'X' should be in {0,1,2} but are in {k,k+1,k+2} for some k in R")
    U2 <- unique(c(X2))
    if( !( (length(U2) <= 3) || (setequal(abs(apply(expand.grid(U2,U2),1,diff)),c(0,1,2))) ) ) 
        stop("Values in 'X2' should be in {0,1,2} or in {k,k+1,k+2} for some k in R")
    else if(!setequal(c(0,1,2),U2))
        warning("Values in 'X2' should be in {0,1,2} but are in {k,k+1,k+2} for some k in R")
    
    if(!is.matrix(X)) X <- as.matrix(X)
    d1 <- dim(X)
    if(!is.double(X)){
        X <- as.double(X)
        dim(X) <- c(d1[1] , d1[2])
    }
    
    if(!is.matrix(X2)) X2 <- as.matrix(X2)
    d2 <- dim(X2)
    if(!is.double(X2)){
        X2 <- as.double(X2)
        dim(X2) <- c(d2[1] , d2[2])
    }
    if(d1[2] != d2[2])stop("Column dimensions of 'X1' and 'X2' do not match.")
    
    K <- matrix(NaN,d1[1],d2[1])    
    
    if(!is.null(para)){
        para <- rep(para,length.out = min(c(d1[2],d2[2])))
        if(!is.double(para)) para <- as.double(para)
        if(any(para < 0)) stop("negative weights encountered for IBS kernel")       
    }
    .Call("Call_ibs2_kernel", X, X2, kernel = as.integer(kernel), para = para,K = K)
    K
    
}


R_ibsX <- function(X,para = NULL,order = 2){
    
    # check data integrity x \in{k,k+1,k+2} for IBS kernel  
    U <- unique(c(X))
    if( !( (length(U) <= 3) || (setequal(abs(apply(expand.grid(U,U),1,diff)),c(0,1,2))) ) ) 
        stop("Values in 'X' should be in {0,1,2} or in {k,k+1,k+2} for some k in R")
    else if(!setequal(c(0,1,2),U))
        warning("Values in 'X' should be in {0,1,2} but are in {k,k+1,k+2} for some k in R")
    
    n <- nrow(X)
    K <- matrix(NA,n,n)
    S <- ncol(X)
    
    if(is.null(para)){
        denom <- order * S
        for(i in 1:n)
            for(j in i:n)
                K[i,j] <- sum(order - abs(X[i,]-X[j,])) / denom
    }else{
        if(any(para < 0)) stop("Some weights are negative")
        para <- rep(para,length.out = S)
        w <- para / sum(para) / order
        for(i in 1:n)
            for(j in i:n)
                K[i,j] <- sum(w * (order - abs(X[i,]-X[j,])))       
    }
    
    for(i in 1:(n-1))
        for(j in (i+1):n)
            K[j,i] <- K[i,j]
    K   
}

R_ibs <- function(X,X2 = NULL,para = NULL,order = 2.0){
    if(is.null(X2))X2 <- X
    
    # check data integrity x \in{k,k+1,k+2} for IBS kernel  
    U <- unique(c(X))
    if( !( (length(U) <= 3) || (setequal(abs(apply(expand.grid(U,U),1,diff)),c(0,1,2))) ) ) 
        stop("Values in 'X' should be in {0,1,2} or in {k,k+1,k+2} for some k in R")
    else if(!setequal(c(0,1,2),U))
        warning("Values in 'X' should be in {0,1,2} but are in {k,k+1,k+2} for some k in R")
    U2 <- unique(c(X2))
    if( !( (length(U2) <= 3) || (setequal(abs(apply(expand.grid(U2,U2),1,diff)),c(0,1,2))) ) ) 
        stop("Values in 'X2' should be in {0,1,2} or in {k,k+1,k+2} for some k in R")
    else if(!setequal(c(0,1,2),U2))
        warning("Values in 'X2' should be in {0,1,2} but are in {k,k+1,k+2} for some k in R")
    
    S <- ncol(X)
    if(S != ncol(X2))stop("Column dimensions of 'X' and 'X2' do not match.")
    n1 <- nrow(X)
    n2 <- nrow(X2)
    K <- matrix(NaN,n1,n2)
    if(is.null(para)){
        denom <- order * S
        for(i in 1:n1)
            for(j in 1:n2)
                K[i,j] <- sum(order - abs(X[i,]-X2[j,]))/denom
    }else{
        if(any(para < 0)) stop("Some weights are negative")
        para <- rep(para,length.out = S)
        w <- para / sum(para) / order
        for(i in 1:n1)
            for(j in 1:n2)
                K[i,j] <- sum(w * (order - abs(X[i,]-X2[j,])))    
    }
    K   
}


R_ibs2 <- function(X,X2 = NULL,para = NULL){
    if(is.null(X2))X2 <- X
    
    # check data integrity x \in{k,k+1,k+2} for IBS kernel  
    U <- unique(X)
    if( !( (length(U) <= 3) || (setequal(abs(apply(expand.grid(U,U),1,diff)),c(0,1,2))) ) ) 
        stop("Values in 'X' should be in {0,1,2} or in {k,k+1,k+2} for some k in R")
    else if(!setequal(c(0,1,2),U))
        warning("Values in 'X' should be in {0,1,2} but are in {k,k+1,k+2} for some k in R")
    U2 <- unique(X2)
    if( !( (length(U2) <= 3) || (setequal(abs(apply(expand.grid(U2,U2),1,diff)),c(0,1,2))) ) ) 
        stop("Values in 'X2' should be in {0,1,2} or in {k,k+1,k+2} for some k in R")
    else if(!setequal(c(0,1,2),U2))
        warning("Values in 'X2' should be in {0,1,2} but are in {k,k+1,k+2} for some k in R")
    
    S <- ncol(X)
    if(S != ncol(X2))stop("Column dimensions of 'X' and 'X2' do not match.")
    n1 <- nrow(X)
    n2 <- nrow(X2)
    K <- matrix(NaN,n1,n2)
    if(is.null(para)){
        for(i in 1:n1)
            for(j in 1:n2)
                K[i,j] <- sum(2.0 - abs(X[i,]-X2[j,]))/(2.0 * S)
    }else{
        if(any(para < 0)) stop("Some weights are negative")
        para <- rep(para,length.out = S)
        w <- para / sum(para) / 2.0
        for(i in 1:n1)
            for(j in 1:n2)
                K[i,j] <- sum(w * (2 - abs(X[i,]-X2[j,])))    
    }
    K   
}


R_hamming.sim <- function(X,X2 = NULL,para = NULL){
    if(is.null(X2))X2 <- X
    S <- ncol(X)
    if(S != ncol(X2))stop("Column dimensions of 'X' and 'X2' do not match.")
    n1 <- nrow(X)
    n2 <- nrow(X2)
    K <- matrix(NaN,n1,n2)
    
    if(is.null(para)){
        for(i in 1:n1)
            for(j in 1:n2)
                K[i,j] <- sum(X[i,] == X2[j,]) / S
    }else{
        if(any(para < 0)) stop("Some weights are negative")
        para <- rep(para,length.out = ncol(X))
        w <- para / sum(para)
        for(i in 1:n1)
            for(j in 1:n2)
                K[i,j] <- sum(w * c(X[i,] == X2[j,]))   
    }
    K   
}

# hamming similarity (as opposed to dissimilarity)
hamming.sim <- function(X,X2 = NULL,para = NULL){
    
    if(is.null(X2))X2 <- X
    
    if(!is.matrix(X)) X <- as.matrix(X)
    d1 <- dim(X)
    if(!is.double(X)){
        X <- as.double(X)
        dim(X) <- c(d1[1] , d1[2])
    }
    
    if(!is.matrix(X2)) X2 <- as.matrix(X2)
    d2 <- dim(X2)
    if(!is.double(X2)){
        X2 <- as.double(X2)
        dim(X2) <- c(d2[1] , d2[2])
    }
    
    if(!is.null(para)){
        para <- rep(para,length.out = min(c(d1[2],d2[2])))
        if(!is.double(para)) para <- as.double(para)
        if(any(para < 0)) stop("negative weights encountered for IBS or HAMMING kernel")        
    }
    
    if(d1[2] != d2[2])stop("Column dimensions of 'x1' and 'X2' do not match.")
    
    K <- matrix(NaN,d1[1],d2[1])
    .Call("Call_hammingSim_kernel",X,X2,para,K)
    K
}

 
