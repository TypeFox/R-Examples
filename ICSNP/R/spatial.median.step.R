`spatial.median.step` <-
function(y,datas)
    { 
    datas.y <- sweep(datas,2,y) 
    eta.y <- 1-as.numeric(apply(datas.y,1,setequal,rep(0,length(y))))
    if (sum(eta.y)==0) stop("all observations are equal")
    
    d.y <- norm(datas.y)
    T.tilde.y <- (sum(1/d.y[eta.y==1]))^(-1)*colSums((sweep(datas[eta.y==1,],1,d.y[eta.y==1],"/")))

    if (length(eta.y)==length(eta.y[eta.y==1])) {
    y.new <- T.tilde.y 
    } else {
    R.tilde.y <- colSums((sweep(datas.y[eta.y==1,],1,d.y[eta.y==1],"/")))
    r.y <- sqrt(sum(R.tilde.y^2))
    s.eta.y.r <- (length(eta.y) -sum(eta.y))/r.y
    y.new <- max(c(0,(1-s.eta.y.r)))*T.tilde.y+ min(c(1,s.eta.y.r))*y
    }
    return(y.new)
    }
