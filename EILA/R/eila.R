
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Subset AIMs
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

comp.anc.dist <- function(anc1,anc2){
    m <- nrow(anc1)
    res <- rep(NA,m)
    for (i in 1:m){
        res[i] <- wilcox.test(anc1[i,],anc2[i,])$p.value
    }
    return(res)
}

get.thin.idx <- function(pvalues, position, V=50000, C=0.2){
    pos.int <- floor((position-min(position)+0.5)/V)

    res <- rep(FALSE,length(position))
    for (i in unique(pos.int)){
        idx <- pos.int != i
        p <- pvalues
        p[idx] <-  1
        res[which.min(p)] <- TRUE
    }

    res[pvalues > C] <- FALSE
    return(res)
}

subset.data <- function(obj,idx){
    m <- length(obj)
    res <- obj
    for (i in 1:m){
        x <- obj[[i]]
        if(is.vector(x)){
            if(length(x) > 1){
                res[[i]] <- x[idx]
            }else{
                res[[i]] <- x
            }
        }else if(is.matrix(x)){
            res[[i]] <- x[idx,]
        }
    }
    return(res)
}

get.AIMs <- function(obj){

    pval <- comp.anc.dist(obj$anc1,obj$anc2)
    idx <- get.thin.idx(pvalues=pval,position=obj$position)
    res <- subset.data(obj,idx)

    return(res)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## dec genotype to probability
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get.prob.score <- function(obj){
    anc1 <- obj$anc1
    anc2 <- obj$anc2
    admixed <- obj$admixed

    old.o <- options(warn= -1)
    r <- nrow(admixed)
    c <- ncol(admixed)
    prob <- matrix(NA, nrow=r, ncol=c)
    y <- c(rep(1,ncol(anc1)),rep(0,ncol(anc2)))
    ANC <- cbind(anc1,anc2)
    for (i in 1:r){
        x <- ANC[i,]
        res.fit <- glm(y~x,family=binomial(link=logit))
        prob[i,] <- predict(res.fit, newdata=data.frame(x=admixed[i,]),
                           type="response")
    }
    options(old.o)

    res <- list(prob=prob,
                position=obj$position)
    return(res)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## partition chromosome
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

partition.chr <- function(pos){

    EPI <- 0.1
    b0 <- min(pos)
    b1 <- max(pos)
    res <- c(b0-EPI,b1+EPI)
    D <- 300
    if(length(pos)> D){
        M <- ceiling(length(pos)/D)
        K <- (b1-b0)/M
        res <- c(res,seq(1,M-1)*K)
        res <- sort(res)
    }
    return(res)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate dec QR est
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


quantreg.est <- function(y, lambda, eps=1e-3, tau=0.5) {
    m <- length(y)
    E <- diag(m)
    D <- diff(E)
    B <- rbind(E, lambda * D)

    uniq.y <- unique(y)
    if(length(uniq.y) < length(y)){
        d0 <- min(diff(sort(uniq.y))) / 100
        y <- y + rnorm(length(y), mean=0, sd=d0)
    }
    yhat <- c(y,rep(0,m-1))
    res.rq.fit <- rq.fit.fnb(B,yhat, eps=eps, tau=tau)
    rq.est <- res.rq.fit$coef
    return(rq.est)
}



prob.to.qr <- function(prob, lambda){
    R <- nrow(prob)
    C <- ncol(prob)
    res <- matrix(NA, nrow=R, ncol=C)

    for (i in 1:C){
        res[,i] <- quantreg.est(prob[,i], lambda=lambda)
    }
    return(res)
}



get.fused.qr <- function(obj,lambda){

    partition.points <- partition.chr(obj$position)

    N <- nrow(obj$prob)
    C <- ncol(obj$prob)
    dec.qr <- matrix(NA,nrow=N,ncol=C)
    for (k in 1:(length(partition.points)-1)){
        ##cat(k," out of ",length(partition.points)," ......\n")
        idx <- obj$position > partition.points[k] &
            obj$position <= partition.points[k+1]
        dec.qr[idx,] <- prob.to.qr(obj$prob[idx,],lambda)

    }
    res <- list(admixed.qr = dec.qr,
                position = obj$position)

    return(res)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get breakpoints(jumps)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.breakpoints <- function(obj){
    dec.qr <- obj$admixed.qr
    position <- obj$position
    EPS <- 10^(-2)
    N <- ncol(dec.qr)
    res <- list()
    for (i in 1:N){
        left.idx <- which(abs(diff(dec.qr[,i])) > EPS)
        res[[i]] <- (position[left.idx] + position[left.idx + 1])/2
    }
    return(res)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## cluster each segment
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simulate.admixed <- function(obj1,obj2){
    n1 <- nrow(obj1)
    n2 <- nrow(obj2)
    n <- min(c(n1,n2))
    idx1 <- sample(n1,n)
    idx2 <- sample(n2,n)
    res <- matrix(NA, nrow=n, ncol=ncol(obj1))
    for (i in 1:n){
        x1 <- obj1[idx1[i],]
        x2 <- obj2[idx2[i],]
        idx1.0 <- x1 == 0
        idx1.1 <- x1 == 1
        idx1.2 <- x1 == 2
        idx2.0 <- x2 == 0
        idx2.1 <- x2 == 1
        idx2.2 <- x2 == 2
        z <- (x1 + x2)/2

        idx <- (idx1.0 & idx2.1) | (idx1.1 & idx2.0)
        if(any(idx)) z[idx] <- sample(c(0,1),sum(idx), replace=T)

        idx <- (idx1.1 & idx2.2) | (idx1.2 & idx2.1)
        if(any(idx)) z[idx] <- sample(c(1,2),sum(idx), replace=T)

        idx <- idx1.1 & idx2.1
        if(any(idx)) z[idx] <- sample(c(0,1,2), size=sum(idx),
                                      replace=T, prob=c(0.25,0.5,0.25))
        res[i,] <- z
    }
    return(res)
}

cluster.segment <- function(obj){

    test <- obj$test
    train <- obj$train
    train.lab <- obj$train.lab
    ini <- obj$ini

    get.mode <- function(x){
        y <- table(x)
        names(y[which.max(y)])
    }

    cluster <- try(kmeans(rbind(train,test),centers=ini)$cluster,T)
    if(class(cluster) == "try-error"){
        res <- knn(train, test, train.lab,k=7)
        res <- as.character(res)
    }else{
        test.cluster <- cluster[length(cluster)]
        train.cluster <- cluster[-length(cluster)]
        names(train.cluster) <- train.lab
        res <- get.mode(names(train.cluster[train.cluster ==
                                            test.cluster]))
        if(length(res) == 0){
            res <- knn(train, test, train.lab,k=7)
            res <- as.character(res)
        }
    }
    res <- as.integer(res)
    return(res)
}


get.train.test.data <- function(obj){
    dec <- obj$admixed
    anc1 <- obj$anc1
    anc2 <- obj$anc2
    anc3.NULL <- obj$anc3.NULL

    admix12 <- simulate.admixed(anc1,anc2)

    if(!anc3.NULL){
        anc3 <- obj$anc3
        admix13 <- simulate.admixed(anc1,anc3)
        admix23 <- simulate.admixed(anc2,anc3)
        n3 <- nrow(anc3)
        n13 <- nrow(admix13)
        n23 <- nrow(admix23)
        ini <- rbind(apply(anc1,2,mean),
                     apply(anc2,2,mean),
                     apply(anc3,2,mean),
                     apply(admix12,2,mean),
                     apply(admix13,2,mean),
                     apply(admix23,2,mean))
    }else{
        anc3 <- NULL
        admix13 <- NULL
        admix23 <- NULL
        n3 <- n13 <- n23 <- 0
        ini <- rbind(apply(anc1,2,mean),
                     apply(anc2,2,mean),
                     apply(admix12,2,mean)
                     )
    }
    train <- rbind(anc1,anc2,anc3,admix12,admix13,admix23)
    train.lab <- c(rep(11,nrow(anc1)),
                   rep(22,nrow(anc2)),
                   rep(33,n3),
                   rep(12,nrow(admix12)),
                   rep(13,n13),
                   rep(23,n23))

    return(list(test=dec,
                train=train,
                train.lab=train.lab,
                ini=ini)
           )
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Infer local ancestry
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subset.segment <- function(obj,i,idx){
    list(test      = obj$test[i,idx,drop=F],
         train     = obj$train[,idx],
         train.lab = obj$train.lab,
         ini       = obj$ini[,idx]
         )
}
transpose.matrix <- function(obj){
    m <- length(obj)
    res <- obj
    for (i in 1:m){
        x <- obj[[i]]
        if(is.matrix(x)) res[[i]] <- t(x)
    }
    return(res)
}

get.local.ancestry <- function(obj, breakpoints){


    obj.t <- transpose.matrix(obj) ## From GxN to NxG
    my.data <- get.train.test.data(obj.t)

    pos <- obj$position
    b0 <- min(pos) - 1
    bm <- max(pos) + 1
    R <- nrow(my.data$test)
    C <- ncol(my.data$test)
    res <- matrix(NA, nrow=R, ncol=C)

    for (i in 1:R){
        ##cat("\ni=",i,"\n")
        bps <- c(b0, breakpoints[[i]], bm)
        m <- length(bps)
        for (b in 1:(m-1)){
            ##cat("\tb=",b,"\n")
            idx <- pos > bps[b]  & pos <= bps[b+1]
            ##cat("\tidx length =",sum(idx),"\n")
            if(sum(idx) > 1){
                obj.ib <- subset.segment(my.data,i,idx)
                res.cluster <- cluster.segment(obj.ib)
                res[i,idx] <- res.cluster
            }
        }
    }
    res <- t(res) ## Transpose back to match the input data format:
    return(res)
}



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MAIN program
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

infer.breakpoints <- function(obj,lambda){
#    cat("\n\nget AIMs ......")
    AIM <- get.AIMs(obj)

#    cat("\n\nget prob score ......")
    PROB <- get.prob.score(AIM)

#    cat("\n\nget fused quantile regression ......\n")
    FQR <- get.fused.qr(PROB, lambda=lambda)

#    cat("\n\nget breakpoints ......")
    res <- get.breakpoints(FQR)
    return(res)
}

infer.local.ancestry1 <- function(obj,lambda){

    Breakpoints <- infer.breakpoints(obj,lambda)

#    cat("\n\ninfer local ancestry ......")
    res <- get.local.ancestry(obj, Breakpoints)
    return(res)

}

infer.local.ancestry2 <- function(obj,obj1,obj2,obj3,lambda){

    Breakpoints1 <- infer.breakpoints(obj1,lambda)
    Breakpoints2 <- infer.breakpoints(obj2,lambda)
    Breakpoints3 <- infer.breakpoints(obj3,lambda)

    m <- length(Breakpoints1)
    Breakpoints <- list()
    for (i in 1:m){
        z <- c(Breakpoints1[[i]],
               Breakpoints2[[i]],
               Breakpoints3[[i]])
        Breakpoints[[i]] <- sort(z)
    }

    cat("infer local ancestry ......\n")
    res <- get.local.ancestry(obj, Breakpoints)
    return(res)
}


eila <- function(admixed,position,anc1,anc2,anc3=NULL,
                 lambda=15,rng.seed=172719943){
    
    rng.state <- NULL
    if(exists(".Random.seed")) rng.state <- .Random.seed
    set.seed(rng.seed)

    anc3.NULL <- ifelse(is.null(anc3), TRUE, FALSE)

    old.op <- options(error = NULL)
    if(!hasArg(admixed))  stop("Argument admixed is missing")
    if(!hasArg(position)) stop("Argument position is missing")
    if(!hasArg(anc1))     stop("Argument anc1 is missing")
    if(!hasArg(anc2))     stop("Argument anc2 is missing")

    if(!is.vector(position)) position <- as.vector(position)
    if(!is.matrix(admixed)) admixed <- as.matrix(admixed)
    if(!is.matrix(anc1)) anc1 <- as.matrix(anc1)
    if(!is.matrix(anc2)) anc2 <- as.matrix(anc2)

    G0 <- nrow(admixed)
    G1 <- nrow(anc1)
    G2 <- nrow(anc2)

    if(length(unique(c(G0,G1,G2))) > 1 )
        stop("Numbers of rows between admxied, anc1, or anc2 do not match")

    if(any(c(ncol(admixed )<5, ncol(anc1)<5, ncol(anc2)<5)))
        stop("Current verison of eila requires at least five columns in admixed, anc1, and anc2")

    if(any(is.na(position))) stop("There is a missing value in position")
    if(any(is.na(admixed)))
        stop("There is a missing value in admixed. Please impute missing (eg. impute package)")
    if(any(is.na(anc1)))
       stop("There is a missing value in anc1. Please impute missing (eg. impute pacakge)")
    if(any(is.na(anc2)))
       stop("There is a missing value in anc2. Please impute missing (eg. impute pacakge)")

    if(any(admixed < 0) || any(admixed > 2))
        stop("The values in anc1 are not valid")
    if(any(anc1 < 0) || any(anc1 > 2))
        stop("The values in anc1 are not valid")
    if(any(anc2 < 0) || any(anc2 > 2))
        stop("The values in anc1 are not valid")

    if(lambda < 0) stop("lambda has to be a positive number")
    

    obj <- list(admixed=admixed,
                position=position,
                anc1=anc1,
                anc2=anc2,
                anc3.NULL=anc3.NULL)

    if(!anc3.NULL){
        if(!is.matrix(anc3)) anc3 <- as.matrix(anc3)
        G3 <- nrow(anc3)
        if(G0 != G3) stop("Numbers of rows between admxied and anc3 do not match")
        if(any(is.na(anc3)))
            stop("There is a missing value in anc3. Please impute missing (eg. impute pacakge)")

        obj <- c(obj, list(anc3=anc3))
        obj1 <- list(admixed=admixed,position=position,anc1=anc1,anc2=anc2)
        obj2 <- list(admixed=admixed,position=position,anc1=anc1,anc2=anc3)
        obj3 <- list(admixed=admixed,position=position,anc1=anc2,anc2=anc3)
    }

    if(anc3.NULL){
        local.ancestry <- infer.local.ancestry1(obj,lambda)
    }else{
        local.ancestry <- infer.local.ancestry2(obj,obj1,obj2,obj3,lambda)
    }
    options(old.op)

    return(list(local.ancestry=local.ancestry,
                rng.seed=rng.seed, rng.state=rng.state))
}
