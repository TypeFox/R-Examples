"crossval"<- function(x,y,theta.fit,theta.predict,...,ngroup=n){
    call <- match.call()
    x <- as.matrix(x)
    n <- length(y)
    ngroup <- trunc(ngroup)
    if( ngroup < 2){
        stop ("ngroup should be greater than or equal to 2")
    }
    if(ngroup > n){
        stop ("ngroup should be less than or equal to the number of observations")
    }
  
    if(ngroup==n) {groups <- 1:n; leave.out <- 1}
    if(ngroup<n){
        leave.out <- trunc(n/ngroup);
        o <- sample(1:n)
        groups <- vector("list",ngroup)
        for(j in 1:(ngroup-1)){
            jj <- (1+(j-1)*leave.out)
            groups[[j]] <- (o[jj:(jj+leave.out-1)])
        }
        groups[[ngroup]] <- o[(1+(ngroup-1)*leave.out):n]
    }
    u <- vector("list",ngroup)
    cv.fit <- rep(NA,n)
    for(j in 1:ngroup){
        u <- theta.fit(x[-groups[[j]], ],y[-groups[[j]]],...)
        cv.fit[groups[[j]]] <-  theta.predict(u,x[groups[[j]],])
        
    }
  
    if(leave.out==1) groups <- NULL
    return(list(cv.fit=cv.fit, 
                ngroup=ngroup, 
                leave.out=leave.out,
                groups=groups, 
                call=call)) 
}
