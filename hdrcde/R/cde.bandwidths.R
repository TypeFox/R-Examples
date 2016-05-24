cde.bandwidths <- function(x,y,deg=0,link="identity",method=1,y.margin,passes=2,
                   ngrid=8,min.a=NULL,ny=25,use.sample=FALSE,GCV=TRUE,b=NULL,...)
{
    if(deg>0 & method>2)
        stop("Method unavailable for degree > 0")

    if(missing(y.margin))
        y.margin <- seq(min(y),max(y),l=ny)
    else
        ny <- length(y.margin)

    if(deg==0)
        return(cde.bandwidths0(x,y,method,y.margin,
                ny=ny,passes=passes,use.sample=use.sample,ngrid=ngrid,...))
    else if(method == 2)
        return(cdeband.rules(x,y,deg=deg,link=link,...))
    else # method = 1
    {
        bands <- cdeband.rules(x,y,deg=deg,link=link,...)
        if(!is.null(b))
            bands$b <- b
        print(paste("Initial values: a=",round(bands$a,3),"  b=",round(bands$b,3)),sep="")
        if(is.null(min.a))
            min.a <- 0.2*bands$a
        a.grid <- bands$a*seq(min.a/bands$a,2,l=ngrid)
        regout <- cdeband.Mbh(x,y,a.grid,bands$b,y.margin=y.margin,use.sample=use.sample,deg=deg,link=link,GCV=GCV,passes=1,ngrid=ngrid,usequad=FALSE,...)
        return(list(a=regout$a,b=bands$b,a.grid=regout$a.grid,q=regout$q))
    }
}


"cde.bandwidths0" <-
function(x,y,method=1,y.margin,sdlinear=FALSE,xden="normal",penalty=4,
    ngrid=8,modified=FALSE,k=3, m=25,nx=30,ny=25,passes=2,tol=0.99999,use.sample=FALSE,usequad=TRUE)
{
    if(missing(y.margin))
        y.margin <- seq(min(y),max(y),l=ny)
    else
        ny <- length(y.margin)

    if(method==2)
        return(cdeband.rules0(x,y,sdlinear,xden,k,modified))
    else if(method==3)
        return(cdeband.regress(x,y,use.sample=use.sample,y.margin=y.margin,
                            penalty=penalty,tol=tol,usequad=usequad))
    else if(method==4)
        return(cdeband.bootstrap(x,y,m=m,nx=nx,y.margin=y.margin))
    else
    {
        firstbands <- bands <- cdeband.rules0(x,y,sdlinear=sdlinear,xden=xden,modified=modified,k=k)
        print(paste("Initial values: a=",round(bands$a,3),"  b=",round(bands$b,3)),sep="")

        for(i in 1:passes)
        {
            if(i==1)
            {
                a.grid <- bands$a*seq(0.1,3,l=ngrid)
                b.grid <- bands$b*seq(0.1,3,l=ngrid)
            }
            else # Choose values around previous minimums
            {
                na <- length(regout$a.grid)
                nb <- length(bootout$b.grid)
                adiff <- abs(bands$a-regout$a.grid)
                idx <- (1:na)[adiff == min(adiff)]
                a.grid <- seq(regout$a.grid[max(idx-2,1)],regout$a.grid[min(idx+2,na)],l=ngrid)
                bdiff <- abs(bands$b-bootout$b.grid)
                idx <- (1:nb)[bdiff == min(bdiff)]
                b.grid <- seq(bootout$b.grid[max(idx-2,1)],bootout$b.grid[min(idx+2,nb)],l=min(8,ngrid))
            }
            regout <- cdeband.regress(x,y,a.grid,bands$b,y.margin=y.margin,tol=tol,use.sample=use.sample,usequad=usequad)
            regna <- is.na(regout$a)
            if(regna)
                bands$a <- firstbands$a
            else
                bands$a <- regout$a
            regna <- is.na(regout$a)
            bootout <- cdeband.bootstrap(x,y,bands$a,b.grid,m,nx,y.margin=y.margin)
            bands$b <- bootout$b
        }
        return(list(a=bands$a,b=bands$b,a.grid=regout$a.grid,b.grid=bootout$b.grid,q=regout$q,imse=bootout$imse,regna=regna))
    }
}

"cdeband.rules0" <- function(x,y,sdlinear=FALSE,xden="normal",k=3,modified=FALSE)
{
    if(xden=="normal" & sdlinear)
        stop("Option not available")
    if(modified & xden!="normal")
        stop("Modified estimator requires a normal marginal density")
    Rk <- 0.5/sqrt(pi)
    sk <- 1
    n <- length(x)
    junk <- Lsquare(x,y,sdlinear)
    p <- junk$pc
    d <- junk$dx
    pl <- NA
    qx <- NA
    if(xden=="uniform")
    {
        u <- max(x)
        l <- min(x)
        if(!sdlinear)
        {
            a <- ((4*sqrt(pi)*Rk^2*(u-l)*abs(p)^5)/(3*n*sk^4*abs(d)^5))^(1/6)
            b <- abs(d)*a
        }
        else
        {
            pl <- junk$pl
            qx <- junk$qx
            z <- ((pl+qx*u)^4 - (pl+qx*l)^4)/((pl+qx*u)^4 *(pl+qx*l)^4)
            w <- 19*qx^4 + 4*d^4 + 28*d^2*qx^2
            a <- ((2^(7.5)*sqrt(pi)*Rk^2*(u-l)^2*qx)/(3*n*sk^4*z*w^(0.75)*(sqrt(w)+2*d^2-3*qx^2)))^(1/6)
            b <- a*w^(0.25)/sqrt(2)
        }
    }
    else if(xden=="normal")
    {
        erf <- 0.954499736
        sm <- sqrt(var(x))
        if(!modified)
        {
            v <- 3*erf*d^2*sm^3*pi - 8*sqrt(2*pi)*p^2*k*exp((-k^2)/2) + 8*pi*p^2*erf
            anum <- ((16*Rk^2*k*(pi)^(1.25)*p^5*sm^(2.5)) / (n*abs(d)^(5/2)*sk^4))^(1/6)
            adem <- (((v^5)/(3*(pi)^2*sm^4*erf))^(1/4) + 3*d*((v*erf^(1/3))/3)^(3/4))^(1/6)
            a <- anum/adem
            b <- ((d^2*v)/(3*pi*sm*erf))^(1/4) * a
        }
        else
        {
            u <- 3*d^2*sm^2 + 4*p^2
            anum <- 16*pi*Rk^2*sm^2.5*p^5
            adem <- n * sk^4 * abs(d)^2.5 * ((u^5/(3*sm^4))^0.25 + (3*d^4*u^3)^0.25)
            a <- (anum/adem)^(1/6)
            b <- a *(d^2*u/3/sm^2)^0.25
        }
    }
    else
        stop("Unknown distribution")
    return(list(a=a,b=b,p=p,d=d,pl=pl,q=qx))
}

"cdeband.rules" <- function(x,y,deg,link,mean.order,...)
{
    if(missing(mean.order))
    {
        if(deg==1)
            mean.order <- 2
        else
            mean.order <- 1
    }
    if(mean.order<1 | mean.order>2 | deg<0 | deg>2)
        stop("Not implemented")
    if(deg==2 & mean.order!=1)
        stop("Not implemented")
    if(deg==0)
    {
        if(mean.order!=1 | link!="identity")
            stop("Not implemented")
        else
            return(cdeband.rules0(x,y,...))
    }

    n <- length(x)
    if(mean.order==1)
    {
        model <- lm(y~x)
            d1 <- abs(model$coef[2])
        d2 <- 0
    }
    else
    {
        mu <- mean(x)
        model <- lm(y~ I(x-mu) + I((x-mu)^2))
        d1 <- abs(model$coef[2])
        d2 <- model$coef[3]
    }
        res <- residuals(model)
        sigma <- sqrt(sum(res^2)/(length(x)-2))

    v <- sqrt(var(x))
    gamma <- 3/(64*pi*v*sigma^5)
    if(deg==1)
    {
        eta <- 0.5
        tau2 <- 0.25/pi
        r1 <- 2
        if(link=="identity")
        {
            alpha <- (3*d1^4+36*d2^2*v^2*(d1^2+v^2*d2^2)+8*d2^2*sigma^2)/(64*pi*sigma^5*v)
            beta <- 3*(d1^2+2*d2^2*v^2)/(32*pi*sigma^5*v)
        }
        else
        {
            alpha <- (d1^4+12*d2^2*v^2*(d1^2+v^2*d2^2)+2*d2^2*sigma^2)/(16*pi*sigma^5*v)
            beta <- (d1^2+2*d2^2*v^2)/(16*pi*sigma^5*v)
        }
        cc <- (alpha/gamma)^0.25
    }
    else
    {
        if(d2!=0)
            stop("Not implemented")
        if(link=="identity")
        {
            alpha <- 105*d1^8/(256*pi*sigma^9*v)
            beta <- 15*d1^4/(64*pi*sigma^7*v)
            cc <- d1^2*sqrt(sqrt(305)+5)/(2*sigma)
            eta <- -0.5
            tau2 <- 27/(64*pi)
            r1 <- 4
        }
        else
            stop("Not implemented")
    }
    power <- 2/(5*r1+2)
    hstar <- (tau2/(n*cc*r1*(2*alpha + beta*cc^2)))^power
    bstar <- cc*(hstar)^(r1/2)

    return(list(a=hstar,b=bstar,sigma=sigma,d1=d1,d2=d2,v=v))
}

cdeband.Mbh <- function(x,y,a.grid,b,ny=50,use.sample=FALSE,nx=100,y.margin,passes=2,deg=deg,
    link="identity",usequad=TRUE,GCV=TRUE,ngrid=10)
{
    ## Initialization

    bands <- cdeband.rules(x,y,deg=deg,link=link)
    if(missing(a.grid))
        a.grid <- bands$a*seq(0.1,5,l=ngrid)
    if(missing(b))
        b <- bands$b
    a.grid <- a.grid[a.grid>0]
    na <- length(a.grid)

    n <- length(x)
    if(use.sample & n > nx)
        idx <- sample(1:n,nx,replace=F)
    else
        idx <- 1:n
    xx <- x[idx]
    yy <- y[idx]
    if(missing(y.margin))
        y.margin <- seq(min(y),max(y),l=ny)

    ## First pass

    first <- CDEband.Mbh(xx,yy,a.grid,b,y.margin,deg,link,GCV=GCV)
    if(is.na(first$a)) # Expand grid
    {
        a.grid <- seq(a.grid[1]/20,1.1*a.grid[na],l=ngrid)
        firstb <- CDEband.Mbh(xx,yy,a.grid,b,y.margin,deg,link,GCV=GCV)
        first$a.grid <- c(first$a.grid,firstb$a.grid)
        idx <- order(first$a.grid)
        first$q <- c(first$q,firstb$q)[idx]
        first$a.grid <- first$a.grid[idx]
        first$a <- firstb$a
    }

    ## Second pass

    if(passes>1 & !is.na(first$a))
    {
        na <- length(first$a.grid)
        idx <- c(1:na)[first$a.grid==first$a]
        if(na>4)
            idx <- idx + (-2:2)
        else
            idx <- idx + (-1:1)
        idx <- idx[idx>=1 & idx<=na]
        newa.grid <- seq(0.95*first$a.grid[idx[1]],1.05*first$a.grid[idx[length(idx)]],l=ngrid)
        second <- CDEband.Mbh(xx,yy,newa.grid,b,y.margin,deg=deg,link=link,GCV=GCV)
        fulla.grid <- c(first$a.grid,second$a.grid)
        idx <- order(fulla.grid)
        fullq <- c(first$q,second$q)[idx]
        fulla.grid <- fulla.grid[idx]
        a <- second$a
    }
    else
    {
        a <- first$a
        fulla.grid <- first$a.grid
        fullq <- first$q
    }

    ## Quadratic solution

    if(!is.na(a) & usequad)
    {
        na <- length(fulla.grid)
        idx <- (c(1:na)[fulla.grid==a])[1]
        if(na>4)
            idx <- idx + (-2:2)
        else
            idx <- idx + (-1:1)
        idx <- idx[idx>=1 & idx<=na & fullq[idx] != Inf]
        if(length(idx)>3)
        {
            fit <- lm(q ~ a+I(a^2),data=data.frame(q=fullq[idx],a=fulla.grid[idx]))
            a <- -0.5*fit$coef[2]/fit$coef[3]
        }
    }
    else
    {
        a <- fulla.grid[fullq==min(fullq,na.rm=T)]
        a <- a[!is.na(a)]
    }

    return(list(a=a,b=b,a.grid=fulla.grid,q=fullq))
}

CDEband.Mbh <- function(x,y,a.grid,b,y.margin,deg=deg,link=link,GCV=TRUE)
{
    na <- length(a.grid)
    q <- numeric(na)
    cat("\n Trying a=")
    for(i in 1:na)
    {
        cat(round(a.grid[i],3)," ")
        junk <- cde(x,y,a=a.grid[i],b=b,x.margin=0,y.margin=y.margin,deg=deg,link=link)
        if(GCV)
            q[i] <- junk$GCV
        else
            q[i] <- junk$AIC
    }
    cat("\n")
    idx <- (1:na)[q==min(q,na.rm=TRUE)]
    idx <- idx[!is.na(idx)]
    if(idx==1 | idx==na)
    {
        warning("No minimum found")
        a <- NA
    }
    else
        a <- a.grid[idx]
    return(list(a=a,a.grid=a.grid,q=q))
}

"cdeband.bootstrap" <- function(x,y,a.grid,b.grid,m=25,nx=30,ny=25,y.margin)
{
    ## Fit parametric model and calculate parametric conditional density
    df <- data.frame(x=x,y=y)
    aic <- numeric(4)
    n <- length(x)
    fit1 <- lm(y ~ 1,data=df)
    aic[1] <- deviance(fit1) + 2*(n-fit1$df.residual)*deviance(fit1)/fit1$df.resid
    fit2 <- lm(y ~ x,data=df)
    aic[2] <- deviance(fit2) + 2*(n-fit2$df.residual)*deviance(fit2)/fit2$df.resid
    fit3 <- lm(y ~ x + I(x^2),data=df)
    aic[3] <- deviance(fit3) + 2*(n-fit3$df.residual)*deviance(fit3)/fit3$df.resid
    fit4 <- lm(y~ x+I(x^2)+I(x^3),data=df)
    aic[4] <- deviance(fit4) + 2*(n-fit4$df.residual)*deviance(fit4)/fit4$df.resid
    fit <- switch((1:4)[aic==min(aic)],fit1,fit2,fit3,fit4)
    rse <- summary(fit)$sigma
    fits <- fitted(fit)
    n <- length(x)
    if(missing(y.margin))
        y.margin <- seq(min(y),max(y),l=ny)
    else
        ny <- length(y.margin)
    if(length(x)<nx)
    {
        x.margin <- sort(x)
        nx <- length(x)
    }
    else
        x.margin <- sort(sample(x,nx))
    fit.grid <- approx(x,fits,xout=x.margin,rule=2)$y  # Quicker than predict and works when constant model used.
    truecde <- list(x=x.margin,y=y.margin,z=matrix(NA,nx,ny))
    for(i in 1:nx)
        truecde$z[i,] <- dnorm(y.margin,fit.grid[i],rse)

    ## Get bandwidth grids

    bands <- cdeband.rules0(x,y)
    if(missing(a.grid))
        a.grid <- bands$a*seq(0.4,2,0.2)
    if(missing(b.grid))
        b.grid <- bands$b*seq(0.4,2,0.2)

    ## Simulate bootstrap samples

    m <- max(m,5)
    bootsamples <- matrix(NA,nrow=n,ncol=m)
    for(i in 1:m)
        bootsamples[,i] <- fits + rnorm(n,0,rse)

    ## Calculate IMSE

    diff.cde <- numeric(m)
    na <- length(a.grid)
    nb <- length(b.grid)
    delta <- y.margin[2]-y.margin[1]
    imse <- matrix(NA,na,nb)
    cat("\n Trying (a,b)=")

    ## First pass

    for(i in 1:na)
    {
        a <- a.grid[i]
        for(j in 1:nb)
        {
            b <- b.grid[j]
            cat("(",round(a,3),",",round(b,3),") ",sep="")
            for(k in 1:3)
            {
                bootcde <- cde(x,bootsamples[,k],deg=0,link="identity",a=a,b=b,x.margin=x.margin,y.margin=y.margin)
                diff.cde[k] <- sum((bootcde$z - truecde$z)^2)/(m*nx)
            }
            imse[i,j] <- mean(diff.cde[1:3])
        }
    }
    cat("\n")

    ## Second pass near minimum

    idx <- (imse==min(imse))
    idx.a <- max(3,(1:na)[apply(idx,1,sum)==1]) + c(-2:2)
    idx.b <- max(3,(1:nb)[apply(idx,2,sum)==1]) + c(-2:2)
    idx.a <- idx.a[idx.a>0 & idx.a <=na]
    idx.b <- idx.b[idx.b>0 & idx.b <=nb]

    for(i in idx.a)
    {
        a <- a.grid[i]
        for(j in idx.b)
        {
            b <- b.grid[j]
            cat("(",round(a,3),",",round(b,3),") ",sep="")
            for(k in 4:m)
            {
                bootcde <- cde(x,bootsamples[,k],deg=0,link="identity",a=a,b=b,x.margin=x.margin,y.margin=y.margin)
                diff.cde[k] <- sum((bootcde$z - truecde$z)^2)/(m*nx)
            }
            imse[i,j] <- (3*imse[i,j] + (m-3)*mean(diff.cde[4:m]))/m
        }
    }
    cat("\n")

    # Find minimum on grid
    idx <- (imse==min(imse))
    idx.a <- (1:na)[apply(idx,1,sum)>0]
    idx.b <- (1:nb)[apply(idx,2,sum)>0]
    best.a <- a.grid[idx.a[1]]
    best.b <- b.grid[idx.b[1]]
    bestmse <- imse[idx.a[1],idx.b[1]]

    ## Fit quadratic surface to points near minimum and find minimum of surface
    na <- length(idx.a)
    nb <- length(idx.b)
    if(na>2 & nb>2)
    {
        fit <- lm(imse ~ a + b + I(a^2) + I(b^2) + I(a*b),
                data=data.frame(imse=c(imse[idx.a,idx.b]),a=rep(a.grid[idx.a],nb),b=rep(b.grid[idx.b],rep(na,nb))))
        best.b <- (fit$coef[2]*fit$coef[6] - 2*fit$coef[3]*fit$coef[4])/(4*fit$coef[4]*fit$coef[5] - fit$coef[6]*fit$coef[6])
        best.a <- (best.b*fit$coef[6] + fit$coef[2])/(-2*fit$coef[4])
    }
    else if(nb>2) # na=1 or 2
    {
        for(i in 1:na)
        {
            fit <- lm(imse ~ b+I(b^2),data=data.frame(imse=imse[idx.a[i],idx.b],b=b.grid[idx.b]))
            bestb <- -0.5*fit$coef[2]/fit$coef[3]
            minimse <- predict(fit,newdata=data.frame(b=bestb))
            if(minimse<bestmse & bestb >= min(b.grid) & bestb <= max(b.grid))
            {
                best.b <- bestb
                best.a <- a.grid[idx.a[i]]
                bestmse <- minimse
            }
        }
    }
    else if(na>2) # nb=1 or 2
    {
        for(i in 1:nb)
        {
            fit <- lm(imse ~ a+I(a^2),data=data.frame(imse=imse[idx.a,idx.b[i]],a=a.grid[idx.a]))
            besta <- -0.5*fit$coef[2]/fit$coef[3]
            minimse <- predict(fit,newdata=data.frame(a=besta))
            if(minimse<bestmse  & besta >= min(a.grid) & besta <= max(a.grid))
            {
                best.a <- besta
                best.b <- b.grid[idx.b[i]]
                bestmse <- minimse
            }
        }
    }

    return(list(a=best.a,b=best.b,a.grid=a.grid,b.grid=b.grid,imse=imse*delta))
}

## This function called by cdeband.regress.

CDEband.regress <- function(x, y, a.grid, b, y.margin, penalty=4, tol=0.999)
{
    # Calculate v
    n <- length(x)
    ny <- length(y.margin)
    if(ny>1)
        delta <- y.margin[2]-y.margin[1]
    else
        delta <- 1
    v <- Kernel(y,y.margin,b)
    na <- length(a.grid)
    q <- numeric(na)
    cat("\n Trying a=")
    for(i in 1:na)
    {
        cat(round(a.grid[i],3)," ")
        w <- Kernel(x,x,a.grid[i])
        row.sum <- apply(w,1,sum)
        w <- w  / matrix(rep(row.sum,n),ncol=n)
        junk <- (v - (t(w) %*% v))^2
        diag.w <- diag(w)
        pen <- switch(penalty,1+2*diag.w,1/(1-diag.w)^2,exp(2*diag.w),(1+diag.w)/(1-diag.w),1/(1-2*diag.w))
        tmp <- apply(junk,1,sum)
        q[i] <- mean(tmp*pen,na.rm=TRUE)
#        if(sum(diag.w>tol)>0)
#            q[i] <- Inf
#        else
#        {
#            pen <- switch(penalty,1+2*diag.w,1/(1-diag.w)^2,exp(2*diag.w),(1+diag.w)/(1-diag.w),1/(1-2*diag.w))
#            q[i] <- mean(apply(junk,1,sum)*pen)
#        }
    }
    cat("\n")
    idx <- (1:na)[q==min(q)]
    if(idx==1 | idx==na)
    {
#        warning("No minimum found")
        a <- NA
    }
    else
        a <- a.grid[idx]
#    browser()
    return(list(a=a,a.grid=a.grid,q=q*delta))
}

cdeband.regress <- function(x,y,a.grid,b,ny=25,use.sample=FALSE,nx=100,y.margin,passes=2,usequad=TRUE,penalty=4,tol=0.999)
{
    ## Initialization

    bands <- cdeband.rules0(x,y)
    if(missing(a.grid))
        a.grid <- bands$a*seq(0.01,3,l=10)
    if(missing(b))
        b <- bands$b
    a.grid <- a.grid[a.grid>0]
    na <- length(a.grid)

    n <- length(x)
    if(use.sample & n > nx)
        idx <- sample(1:n,nx,replace=FALSE)
    else
        idx <- 1:n
    xx <- x[idx]
    yy <- y[idx]
    if(missing(y.margin))
        y.margin <- seq(min(y),max(y),l=ny)

    ## First pass

    first <- CDEband.regress(xx,yy,a.grid,b,y.margin,penalty=penalty,tol=tol)
    if(is.na(first$a)) # Expand grid
    {
        a.grid <- seq(a.grid[1]/20,1.1*a.grid[na],l=50)
        firstb <- CDEband.regress(xx,yy,a.grid,b,y.margin,penalty=penalty,tol=tol)
        first$a.grid <- c(first$a.grid,firstb$a.grid)
        idx <- order(first$a.grid)
        first$q <- c(first$q,firstb$q)[idx]
        first$a.grid <- first$a.grid[idx]
        first$a <- firstb$a
    }

    ## Second pass

    if(passes>1 & !is.na(first$a))
    {
        na <- length(first$a.grid)
        idx <- c(1:na)[first$a.grid==first$a]
        if(na>4)
            idx <- idx + (-2:2)
        else
            idx <- idx + (-1:1)
        idx <- idx[idx>=1 & idx<=na]
        newa.grid <- seq(0.95*first$a.grid[idx[1]],1.05*first$a.grid[idx[length(idx)]],l=length(a.grid))
        second <- CDEband.regress(xx,yy,newa.grid,b,y.margin,tol=tol,penalty=penalty)
        fulla.grid <- c(first$a.grid,second$a.grid)
        idx <- order(fulla.grid)
        fullq <- c(first$q,second$q)[idx]
        fulla.grid <- fulla.grid[idx]
        a <- second$a
    }
    else
    {
        a <- first$a
        fulla.grid <- first$a.grid
        fullq <- first$q
    }

    ## Quadratic solution

    if(!is.na(a) & usequad)
    {
        na <- length(fulla.grid)
        idx <- (c(1:na)[fulla.grid==a])[1]
        if(na>4)
            idx <- idx + (-2:2)
        else
            idx <- idx + (-1:1)
        idx <- idx[idx>=1 & idx<=na & fullq[idx] != Inf]
        if(length(idx)>3)
        {
            fit <- lm(q ~ a+I(a^2),data=data.frame(q=fullq[idx],a=fulla.grid[idx]))
            besta <- -0.5*fit$coef[2]/fit$coef[3]
            if(besta >0)
                a <- besta
        }
    }

    return(list(a=a,b=b,a.grid=fulla.grid,q=fullq))
}

"Lsquare" <- function(x,y,sdlinear=TRUE)
{
    # Fits linear model.
    # If sdlinear=TRUE, allows heteroskedastic errors.

    x <- c(x)
    y <- c(y)
    model <- lm(y~x)
    res <- residuals(model)
    pc <- sqrt(sum(res^2)/(length(x)-2))
    if(!sdlinear)
        return(list(pc=pc,dx=model$coef[2]))

    ## Initial estimate of heteroskedasticity
    model.res <- lm(abs(res)~x)

    ## Recalculate linear fit
    weights <- fitted(model.res)^(-2)
    model <- lm(y~x,weights=weights)

    ## Recalculate heteroskedasticity
    res <- residuals(model)
    f1 <- function(co,res,x)
    {
        sum((res * res - (co[1] + co[2] * x)^2)^2)
    }
    junk <- nlm(f1,model.res$coef,res=res,x=x)
    return(list(pc=pc,pl=junk$estimate[1],qx=junk$estimate[2],dx=model$coef[2]))
}
