
#get the auc through kernel smoothing
#References: 

#For the CDF method
#Using Smoothed Receiver Operating Characteristic Curves to 
#Summarize and Compare Diagnostic Systems
#Chris J. Lloyd
#Journal of the American Statistical Association, 
#Vol. 93, No. 444 (Dec., 1998), pp. 1356-1364

#For the PDF method
#STATISTICS IN MEDICINE, VOL. 16, 2143-2156 (1997)
#SMOOTH NON-PARAMETRIC RECEIVER OPERATING
#CHARACTERISTIC (ROC) CURVES FOR CONTINUOUS DIAGNOSTIC TESTS
#KELLY H. ZOU, W. J. HALL AND DAVID E. SHAPIRO

#For the coomparison of bandwidth
#Statist. Med. 2002; 21:2045-2055 (DOI: 10.1002/sim.1156)
#Comparison of bandwidth selection methods for kernel
#smoothing of ROC curves
#Xiao-Hua Zhou and Jaroslaw Harezlak

#The band width "nrd0" is from Zou et al. (1997)
#The band width "al" is from Zhou and Harezlak (2002)

#x: the scores of subjects from class P
#y: the scores of subjects from class N
#alpha: type I error

#NOTE: the larger the score, the more likely a subject is from class P

getBandwidth.nrd0 <- function(x){
    bw.nrd0(x)
}

getBandwidth.sj <- function(x){
    bw <- try(bw.SJ(x), silent=TRUE)
    if(inherits(bw, "try-error")){
        bw <- bw.nrd0(x)
    }
    bw
}

getBandwidth <- function(x, bw=c("nrd0", "sj")){
    bw <- match.arg(bw)
    switch(bw,
           nrd0 = getBandwidth.nrd0(x),
           sj = getBandwidth.sj(x))
}

#get the variance V_2(\theta) from method 4 of Newcombe 2006
getVARmw  <- function(theta, m, n){
    nstar <- mstar <- (m+n)/2  - 1
    theta * 
        (1-theta) * 
        (1 + nstar*(1-theta)/(2-theta) + mstar*theta/(1+theta)) /
        (m*n)    
    
}

getAUCkernelCDF <- function(x, y, nx, ny, bw){
    hx <- getBandwidth(x, bw)
    hy <- getBandwidth(y, bw)
    h <- sqrt(hx^2 + hy^2)
    xy <- expand.grid(x, y)
    mean(pnorm((xy[,1] - xy[,2]) / h))
}

#get the roc curve and then the auc using trapezoidal rule
#nint: the number of equally spaced points at which the density is to be estimated. 

getAUCkernelPDF <- function(x, y, nx, ny, bw, nint=512){
    
    #get bandwidth
    hx <- getBandwidth(x, bw)
    hy <- getBandwidth(y, bw)
    
    #get range of thresholds
    thr <- rep(0, 2)
    xys <- sort(unique(c(x, y)))
    thr[1] <- min(xys) - max(c(3*hx, 3*hy))
    thr[2] <- max(xys) + max(c(3*hx, 3*hy))
    
    #get first the pdf 
    #second the cdf by numerical intergration using trapezoidal rule
    #third true positive (tp) and false positive (fp). They are equal to 1-cdf
    pdf.x <- density(x, bw=hx, kernel="gaussian", n=nint, from=thr[1], to=thr[2])
    pdf.y <- density(y, bw=hy, kernel="gaussian", n=nint, from=thr[1], to=thr[2])
    dint <- pdf.x$x[2] - pdf.x$x[1]
    cdf.x <- c(0, cumsum((pdf.x$y[1:(nint-1)] + pdf.x$y[2:nint])*dint/2))
    cdf.x <- cdf.x / max(cdf.x)
    tp <- 1- cdf.x
    
    cdf.y <- c(0, cumsum((pdf.y$y[1:(nint-1)] + pdf.y$y[2:nint])*dint/2))
    cdf.y <- cdf.y / max(cdf.y)
    fp <- 1-cdf.y
    
    
    #get the auc
    diffs.x <- fp[-1] - fp[-nint]
    means.vert <- (tp[-1] + tp[-nint])/2
    auc <- sum(means.vert * diffs.x)
    -1*auc
}

getAUCkernel <- function(x, y, nx, ny, method=c("CDF", "PDF"), bw, nint=512){
    
    method <- match.arg(method)
    point <- switch(method,
                    CDF = getAUCkernelCDF(x, y, nx, ny, bw),
                    PDF = getAUCkernelPDF(x, y, nx, ny, bw, nint))
    point
}

auc.kernel.mw <- function(x, y, alpha, method, 
                         bw, nint){

    nx <- length(x)
    ny <- length(y)
    
    point <- getAUCkernel(x, y, nx, ny, method, bw, nint)
    theta <- getAUCmw(x, y)
    se <- sqrt(getVARmw(theta, ny, nx))
    zalpha <- qnorm(alpha/2)
    ci <- c(1-(1-point)*exp(-zalpha * se / (1-point)), 
            1-(1-point)*exp(zalpha * se / (1-point)))
    
    c(point, ci)
    
}

kernel.jackknife <- function(x, y, method, bw, nint){

    nx <- length(x)
    ny <- length(y)
    n <- nx + ny

    hatThetaPartial <- hatThetaPseudo <- rep(0, n)
    for(i in 1:nx){
        hatThetaPartial[i] <- getAUCkernel(x[-i], y, nx-1, ny, method, bw, nint)
    }
    for(i in 1:ny){
        hatThetaPartial[i+nx] <- getAUCkernel(x, y[-i], nx, ny-1, method, bw, nint)
    }

    hatThetaPartial
}

auc.kernel.jackknife <- function(x, y, alpha, method, 
                                 bw, nint){
    
    nx <- length(x)
    ny <- length(y)
    n <- nx + ny
    
    hatTheta <- getAUCkernel(x, y, nx, ny, method, bw, nint)

    hatThetaPseudo <- rep(0, n)

    hatThetaPartial <- kernel.jackknife(x, y, method, bw, nint)
    
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

auc.kernel.boot <- function(x, y, alpha, method, 
                            bw, nint, nboot, cimethod){
    
    nx <- length(x)
    ny <- length(y)
    
    point <- getAUCkernel(x, y, nx, ny, method, bw, nint)
    index.x <- matrix(sample.int(nx, size = nx*nboot, replace = TRUE),
                      nboot, nx)
    index.y <- matrix(sample.int(ny, size = ny*nboot, replace = TRUE),
                      nboot, ny)
    kernel.boot <- sapply(1:nboot, function(i) getAUCkernel(x[index.x[i,]], 
                                                     y[index.y[i,]],
                                                     nx, ny, method, bw, nint))

    if(cimethod=="P"){
        ci <- as.vector(quantile(kernel.boot, c(alpha/2, 1-alpha/2), type=6, na.rm=TRUE))
    }
    else{
        hatZ0 <- qnorm(mean(kernel.boot < point))
            
        partial <- kernel.jackknife(x, y, method, bw, nint)
        mpartial <- mean(partial)
        hatA <- sum((mpartial - partial)^3) /
          (6 * (sum((mpartial - partial)^2))^(3/2))
        
        alpha1 <- pnorm(hatZ0 + (hatZ0 + qnorm(alpha/2)) /
                        (1 - hatA*(hatZ0 + qnorm(alpha/2))))
        alpha2 <- pnorm(hatZ0 + (hatZ0 + qnorm(1-alpha/2)) /
                        (1 - hatA*(hatZ0 + qnorm(1-alpha/2))))
        
        ci <- as.vector(quantile(kernel.boot, c(alpha1, alpha2), type=6))
        
    }
    c(point, ci)
    
}

auc.nonpara.kernel <- function(x, y, conf.level=0.95,
                               integration=c("FALSE", "TRUE"),
                               bw=c("nrd0", "sj"), nint=512, 
                               method=c("mw", "jackknife", "bootstrapP", "bootstrapBCa"), 
                               nboot){

    alpha <- 1 - conf.level
    infer <- match.arg(method)
    method <- match.arg(integration)
    method <- ifelse(method, "PDF", "CDF")
    estimate <- switch(infer,
                       mw=auc.kernel.mw(x, y, alpha, method, bw, nint),
                       jackknife=auc.kernel.jackknife(x, y, alpha, method, bw, nint),
                       bootstrapP=auc.kernel.boot(x, y, alpha, method, bw, nint, nboot, cimethod="P"),
                       bootstrapBCa=auc.kernel.boot(x, y, alpha, method, bw, nint, nboot, cimethod="BCa"))
    estimate
}
