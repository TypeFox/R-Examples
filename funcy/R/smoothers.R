#  Original sm1 sm2: Copyright (C) Chiou J-M and Li P-L.2007
#  Modified Code: Copyright (C) 2011-2015 Christina Yassouridis
# 
#

sm1 <- function(x, y, h, eval.points, weights, poly.index,
                eval.grid=FALSE, display=NULL){
    m <- length(eval.points)
    n <- length(x)
   
    aa <- 1
    poly.index2 <- max(poly.index,1)
    x_sort <- rep(0,n)
    you <- rep(0,m)
    x_sort[1:n-1] <- x[2:n]
  
    if(any(x_sort[1:n-1] < x[1:n-1])==1){
        srt <- sort(x, index.return=TRUE)
        b1 <- srt$x
        ix1 <-srt$ix
        x[1:n] <- b1
        y[1:n] <- y[ix1[1:n]]
        weights[1:n] <- weights[ix1[1:n]]
    }

    id <- weights > 0
    tempx <- x[id]
    tempy <- y[id]
    tempw <- weights[id]

    for(j in 1:m){
        ##calculate bandwidth 
        uu <- eval.points[j]
        ul <- uu - aa*h
        uh <- uu + aa*h
        you[j] <- -99
        count <- 0
        ##find all x,y values within bandwidth
        id <- tempx >= ul
        aux1 <- tempx[id]
        aux2 <- tempy[id]
        aux3 <- tempw[id]
        id <-  aux1 <= uh
        aux1 <- (aux1[id]-uu)/h  #((x_ij - x)/h)
        aux2 <- aux2[id]
        aux3 <- aux3[id]
        
        count <- length(aux2)
        if(count < (2 + poly.index2)){
            if ( poly.index2 == 1){
                if(count == 1)  you[j] <- aux2[1]
                if(count == 2){
                    if(aux1[1] >= 0)  you[j] <- aux2[1]
                    if(aux1[2] <= 0)  you[j] <- aux2[2]
                    if(aux1[1] < 0 & aux1[2] > 0){
                        xh <- aux1[2] - aux1[1]
                        if(aux2[2] > aux2[1]){
                            you[j] <- ( -aux1[1] / xh * ( aux2[2] - aux2[1] ) ) + aux2[1]
                        } else if(aux2[2] < aux2[1]){
                            you[j] <- ( aux1[2] / xh * ( aux2[1] - aux2[2] ) ) + aux2[2]
                        } else
                            you[j] <- (aux2[1] + aux2[2])/2
                        
                    }
                }
            }
            next
        }

        aux3= aux3*(1-aux1*aux1)
        aux3[aux3<0] <- 0
        aux8 <- matrix(0,count,poly.index2+1)
        aux8[,1] <- 1
        aux1 <- aux1*h
        
        for(l in 1:poly.index2)
            aux8[1:length(aux1),l+1] <- as.numeric(aux1^l)
        
        nxmat <- poly.index2 + 1
        aux7 <- rep(0,nxmat)  
        swt <- sqrt(aux3[1:count])
        b <- aux2[1:count]*swt
        SW <- matrix(swt,nrow= length(swt),ncol=nxmat)
        A <- aux8*SW
        aux7 <- lm.fit(A, b)$coefficients
        you[j] <- aux7[1] 
    }

    list(estimate=you)
}

##a 2-dimensional smoother
sm2 <- function(x, y, h, eval.points, weights, poly.index,
                eval.grid=FALSE, display=NULL){

    
    m <- dim(eval.points)[1]

    spoly <- 0
    for(l in 1:2){
        poly.index[l] <- max(poly.index[l],1)
        spoly <- spoly + poly.index[l]
    }
    
    you <- rep(0,m)
    id <- weights > 0
    
    tempx <- x[id,];
    tempy <- y[id];
    tempw <- weights[id]
    for(j in 1:m){
        you[j] <- 0 
        count <- 0  
        x <- tempx
        y <- tempy
        w <- tempw
        aa <- 1
        
        ##extract x and y points within certain bandwidth
        for(l in 1:2){
            if(is.null(dim(x)))
                x <- t(x)
            id1 <- (x[,l] >= (eval.points[j,l] - aa*h[l]))
            x <- x[id1,];
            y <- y[id1];
            w <- w[id1]
            if(is.null(dim(x)))
                x <- t(x)
            id2 <- (x[,l] <= (eval.points[j,l] + aa*h[l]))
            x <- x[id2,];
            y <- y[id2];
            w <- w[id2]
        }

        
        count <- length(y)
        
        if(count==0)
            next
        temp <- as.numeric(eval.points[j,])
        temp <- t(matrix(temp, nrow=length(temp), ncol=count))
        tempd <- t(h)           
        tempd <- t(matrix(tempd, nrow=length(tempd), ncol=count))
        x <- (x - temp)/tempd
        
        we <- rep(1,count)
        we <- we*apply(1-x^2,1,prod) #1-x^2 = Kernel function
        
        we <- w*we
        we[we<0] <- 0 
        temp <- t(h)
        x <- x*matrix(temp,nrow=count,ncol=dim(temp)[2])
        
        xmat <- rep(1,count)
        for(l in 1:2){
            temp1 <- x[,l]
            temp2 <- NULL
            for(i1 in 1:poly.index[l]){
                temp2 <- cbind(temp2,temp1^i1)
            }
            xmat <- cbind(xmat,temp2)
        }
        nxmat <- dim(xmat)[2]
        auy <- rep(0,nxmat) 
        swt <- sqrt(we[1:count])
        b <- y[1:count]*swt
        SW <-matrix(swt,nrow=length(swt),ncol=nxmat)
        A <- SW*xmat
        
        auy <- lm.fit(A, b)$coefficients
        you[j] <- auy[1]
    }
    return(list(estimate=you))
}


##select BW
selBw <- function(data, reg){
    ##Reformat data
    if(!is.null(data)){
        if(reg==0){
            res <- formatFuncy(data=data, format="Format3")
            time <- data[,3]
            y <- data[,2] 
        }else{
            res <- formatFuncy(data=data, format="Format3")
            time <-
                rep(seq(0,dim(data)[1],length=dim(data)[1]),dim(data)[2])
            y <- as.vector(data)
        }
        Tin <- res$Tin; Yin <- res$Yin; isobs <- res$isobs; t_all <- res$t_all
    }
    hmu.cv <- h.select(time, y, method="cv")
    
    ##smoother for mu
    fitmu <- sm.regression(time, y, h=hmu.cv, poly.index=1,
                           eval.points=t_all, display="none")$estimate
    ctrY <- longCtrY(Yin, Tin, t_all, isobs, fitmu)
    conver <- longCov(ctrY, Tin, isobs, t_all)
    Cova <-conver$cov.raw
    CovT <-conver$cov.time
    
    ##select smoothing parameter for covariance
    hcov.cv <- h.select(CovT, Cova, method="cv")
    
    return(list(h1Dim=hmu.cv, h2Dim=hcov.cv))
}
