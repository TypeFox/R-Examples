
#get the auc through Mann-Whitney statistics 
#References:
#Comparing the Areas under Two or More 
#Correlated Receiver Operating Characteristic Curves: 
#A Nonparametric Approach
#Author(s): Elizabeth R. DeLong, David M. DeLong and Daniel L. Clarke-Pearson
#Source: Biometrics, Vol. 44, No. 3 (Sep., 1988), pp. 837-845

#STATISTICS IN MEDICINE
#Statist. Med. 2006; 25:559-573
#Confidence intervals for an effect size measure based on
#the Mann-Whitney statistic. Part 2: Asymptotic methods
#and evaluation
#Robert G. Newcombe

#x: the scores of subjects from class P
#y: the scores of subjects from class N
#alpha: type I error

#NOTE: the larger the score, the more likely a subject is from class P


getAUCmw <- function(x, y){
    
    xy <- expand.grid(x, y)
    mean(ifelse(xy[,2] < xy[,1], 1, (ifelse(xy[,2] == xy[,1], 1/2, 0)))) 
        
}

#get the variance V_2(\theta) from method 4 of Newcombe 2006
getVARmw  <- function(theta, m, n){
    nstar <- mstar <- (m+n)/2  - 1
    theta * 
        (1-theta) * 
        (1 + nstar*(1-theta)/(2-theta) + mstar*theta/(1+theta)) /
        (m*n)    
    
}

#note that here x < y
#m is the length of x and n is the length of y
getS2mw <- function(x, y, m, n){
    
    xi <- sort(x)
    yj <- sort(y)
    rxy <- rank(c(xi, yj))
    Ri <- rxy[1:m]
    Sj <- rxy[m+1:n]
    S102 <- 1/((m-1)*n^2) * (sum((Ri-1:m)^2) - m*(mean(Ri)-(m+1)/2)^2)
    S012 <- 1/((n-1)*m^2) * (sum((Sj-1:n)^2) - n*(mean(Sj)-(n+1)/2)^2)
    S2 <- (m*S012 + n*S102) / (m + n)
    S2
    
}


#get the CI by method 5 of Newcombe 2006
getAUCCImwu  <- function(hat.theta, zalpha, m, n){
    
    nstar <- mstar <- (m+n)/2  - 1
    
    a <- -1 - nstar - mstar
    b <- 1 + 2*mstar
    c <- 2 + nstar
    d <- 1 + 2*hat.theta
    ht2 <- hat.theta^2
    e <- 2 - 2*hat.theta - ht2
    f <- ht2 - 4*hat.theta
    g <- 2*ht2
    
    z2 <- zalpha^2
    mn <- m*n
    z5 <- -mn  + z2*a
    z4 <- mn*d - z2*(a-b)
    z3 <- mn*e - z2*(b-c)
    z2 <- mn*f - z2*c
    z1 <- mn*g
    
    roots <- polyroot(c(z1, z2, z3, z4, z5))
    real <- Re(roots[sapply(1:4, function(i) all.equal(Im(roots[i]), 0))])
    real <- real[real > 0 & real < 1]
    ci <- sort(real)
    if(length(ci) > 2) 
        warning("There are three roots meet the requirement when computing
                the confidence interval.")
    else{
        if(length(real) == 1){
            if(real <= hat.theta){
                ci <- c(real, 1)
            }
            else{
                ci <- c(0, real)
            }
        }
        if(length(real) == 0){
            ci <- c(0, 1)
        }
    }
    ci
}



#get the CI by method 5 of Newcombe 2006
auc.mw.newcombe <- function(x, y, alpha){
    
    point <- getAUCmw(x, y)
    nx <- length(x)
    ny <- length(y)
    zalpha <- qnorm(1-alpha/2)
    ci <- getAUCCImwu(point, zalpha, ny, nx)
    
    c(point, ci)
}

auc.mw.zhou <- function(x, y, alpha){
    if(max(y) < min(x)){
        c(1, 1, 1)
    }
    else{
        point <- getAUCmw(x, y)
        nx <- length(x)
        ny <- length(y)
        zalpha <- qnorm(1-alpha/2)
        varHatTheta <- getVARmw(point,nx, ny)
        
        Z <- 1/2 * log((1+point)/(1-point))
        varZ <- 4 / (1-point^2)^2 * varHatTheta
        LL <- Z - zalpha*sqrt(varZ)
        UL <- Z + zalpha*sqrt(varZ)
        ci <- c((exp(2*LL) - 1)/ (exp(2*LL) + 1), (exp(2*UL) - 1)/ (exp(2*UL) + 1))
        
        c(point,ci)
    }
}

auc.mw.pepe <- function(x, y, alpha){

    if(max(y) < min(x)){
        c(1, 1, 1)
    }
    else{
        point <- getAUCmw(x, y)
        nx <- length(x)
        ny <- length(y)
        zalpha <- qnorm(1-alpha/2)
        
        #varHatTheta <- (nx+ny) * getS2mw(y, x, ny, nx) / (nx*ny)
        varHatTheta <- getVARmw(point,nx, ny)
        LL <- log(point/(1-point)) - zalpha*sqrt(varHatTheta)/(point*(1-point))
        UL <- log(point/(1-point)) + zalpha*sqrt(varHatTheta)/(point*(1-point))
        ci <- c(exp(LL) / (1+exp(LL)), exp(UL) / (1+exp(UL)))
        
        c(point,ci)
    }
}

auc.mw.delong <- function(x, y, alpha){
    
    point <- getAUCmw(x, y)
    nx <- length(x)
    ny <- length(y)
    zalpha <- qnorm(1-alpha/2)
    
    D10 <- sapply(1:ny, function(i)
                  mean(ifelse(x > y[i], 1, ifelse(x == y[i], 1/2, 0))))
    D01 <- sapply(1:nx, function(i)
                  mean(ifelse(x[i] > y, 1, ifelse(x[i] == y, 1/2, 0))))
    varDhatTheta <- 1/(ny*(ny-1))*sum((D10-point)^2) + 
                    1/(nx*(nx-1))*sum((D01-point)^2)
    ci <- c(point - zalpha*sqrt(varDhatTheta),
            point + zalpha*sqrt(varDhatTheta))
    
    c(point, ci)
}

mw.jackknife <- function(x, y){
    nx <- length(x)
    ny <- length(y)
    n <- nx + ny
    
    hatThetaPartial <- rep(0, n)
    
    for(i in 1:nx){
        hatThetaPartial[i] <- getAUCmw(x[-i], y)
    }
    for(i in 1:ny){
        hatThetaPartial[i+nx] <- getAUCmw(x, y[-i])
    }

    hatThetaPartial

}

auc.mw.jackknife <- function(x, y, alpha){
  
    nx <- length(x)
    ny <- length(y)
    n <- nx + ny
    
    hatTheta <- getAUCmw(x, y)
    
    hatThetaPseudo <- rep(0, n)

    hatThetaPartial <- mw.jackknife(x, y)
    
    for(i in 1:n){
        hatThetaPseudo[i] <- n*hatTheta - (n-1)*hatThetaPartial[i]
    }

    point <- mean(hatThetaPseudo)
    ST2 <- mean((hatThetaPseudo - point)^2) / (n-1)
    ST <- sqrt(ST2)
    z.alpha2 <- qt(1 - alpha/2, df=n-1)
    ci <- c(point - z.alpha2*ST, point + z.alpha2*ST)
   
    c(point, ci)
}

auc.mw.boot <- function(x, y, alpha, nboot=1000, method){
    if(max(y) < min(x)){
        c(1, 1, 1)
    }
    else{
        nx <- length(x)
        ny <- length(y)
        
        point <- getAUCmw(x, y)
        index.x <- matrix(sample.int(nx, size = nx*nboot, replace = TRUE),
                          nboot, nx)
        index.y <- matrix(sample.int(ny, size = ny*nboot, replace = TRUE),
                          nboot, ny)
        mw.boot <- sapply(1:nboot, function(i) getAUCmw(x[index.x[i,]], 
                                                        y[index.y[i,]]))
        if(method=="P"){
            ci <- as.vector(quantile(mw.boot, c(alpha/2, 1-alpha/2), type=6))
        }
        else{
            hatZ0 <- qnorm(mean(mw.boot < point))
            
            partial <- mw.jackknife(x, y)
            mpartial <- mean(partial)
            hatA <- sum((mpartial - partial)^3) /
              (6 * (sum((mpartial - partial)^2))^(3/2))
            
            alpha1 <- pnorm(hatZ0 + (hatZ0 + qnorm(alpha/2)) /
                            (1 - hatA*(hatZ0 + qnorm(alpha/2))))
            alpha2 <- pnorm(hatZ0 + (hatZ0 + qnorm(1-alpha/2)) /
                            (1 - hatA*(hatZ0 + qnorm(1-alpha/2))))
            
            ci <- as.vector(quantile(mw.boot, c(alpha1, alpha2), type=6))
        }
        c(point, ci)
    }
}

auc.nonpara.mw <- function(x, y, conf.level=0.95, 
                           method=c("newcombe", "pepe", "delong", "jackknife", "bootstrapP", "bootstrapBCa"), 
                           nboot){

    alpha <- 1 - conf.level
    method <- match.arg(method)
    estimate <- switch(method,
                       newcombe=auc.mw.newcombe(x, y, alpha),
                       pepe=auc.mw.pepe(x, y, alpha),
                       delong=auc.mw.delong(x, y, alpha),
                       jackknife=auc.mw.jackknife(x, y, alpha),
                       bootstrapP=auc.mw.boot(x, y, alpha, nboot, method="P"),
                       bootstrapBCa=auc.mw.boot(x, y, alpha, nboot, method="BCa"))
    estimate
}
