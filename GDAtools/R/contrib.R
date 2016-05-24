contrib <-
 function(resmca) {
    s <- vector()
    for(i in 1:ncol(resmca$call$X)) s <- c(s,rep(i,times=length(levels(resmca$call$X[,i]))))
    e <- 99999
    classe <- class(resmca)[1] # new
    #if(classe=='stMCA') classe=resmca$call$input.mca # new
    if(classe %in% c('speMCA','csMCA')) e <- resmca$call$excl # new
    s <- s[-e]
    dims <- paste('dim',1:resmca$call$ncp,sep='.')
    x <- aggregate(resmca$var$contrib,list(s),sum)[,-1]
    dimnames(x) <- list(colnames(resmca$call$X),dims)
    Z <- dichotom(resmca$call$X,out='numeric')[,-e]
    fK <- colSums(resmca$call$row.w*Z)/nrow(Z)
    Q <- ncol(resmca$call$X)
    ctr.cloud <- data.frame(100*(1-fK)/(ncol(Z)-Q))
    colnames(ctr.cloud) <- 'ctr.cloud'
    rownames(ctr.cloud) <- rownames(resmca$var$contrib)
    vctr.cloud <- aggregate(ctr.cloud,list(s),FUN=sum)[-1]
    colnames(vctr.cloud) <- 'vctr.cloud'
    rownames(vctr.cloud) <- colnames(resmca$call$X)
    list(ctr=round(resmca$var$contrib,6),var.ctr=round(x,6),ctr.cloud=round(ctr.cloud,6),vctr.cloud=round(vctr.cloud,6))
    }
