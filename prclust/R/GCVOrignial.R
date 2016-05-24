###################################
## general degree of freedom and GCV
###################################

GDF.Step23 <- function(seed,data,lambda1,lambda2,tau,mumethods, methods,sigma,algorithm,epsilon)
{
    set.seed(seed) ## Set the random seed for this simulation
    # peturbation of data
    deltaB = matrix(rnorm(dim(data)[1]*dim(data)[2],0,sigma),dim(data)[1],dim(data)[2])
    data1 = data + deltaB
    if (algorithm == 1){
        rho = lambda1
        a = .Call('prclust_PRclustADMM', PACKAGE = 'prclust', data1, rho, lambda2, tau, mumethods, methods,epsilon)
    } else {
        a =    .Call('prclust_PRclustOriginal', PACKAGE = 'prclust', data1, lambda1, lambda2, tau,mumethods, methods)

    }
    out = list()
    out[[1]] = deltaB
    out[[2]] = a$mu
    out[[3]] = a$group
    out
}

##########################################

##########################################

GCV <- function(data,lambda1,lambda2,tau,sigma,B=100,loss.method = c("quadratic","lasso"),group.method = c("gtlp","lasso","SCAD","MCP"), algorithm = c("ADMM","Quadratic"), epsilon = 0.001)
{   ## judge for different situation
    mumethods = switch(match.arg(loss.method), `quadratic` = 0,lasso = 1)
    methods = switch(match.arg(group.method), `gtlp` = 0,lasso = 1, MCP = 2, SCAD = 3)
    nalgorithm = switch(match.arg(algorithm), `ADMM` = 1,Quadratic = 2)
    
    if(is.character(lambda1))
    stop("lambda1 must be a number")
    if(is.character(sigma))
    stop("sigma must be a number")
    if(is.character(B))
    stop("B must be a number")
    if(is.character(lambda2))
    stop("lambda2 must be a number")
    if(is.character(tau))
    stop("tau must be a number")
    
    if(lambda1<0 | is.na(lambda1))
    stop("lambda1 must be a postive number.")
    if(lambda2<0 | is.na(lambda2))
    stop("lambda2 must be a postive number.")
    if(tau<0 | is.na(tau))
    stop("tau must be a postive number.")
    if(sigma<0 | is.na(sigma))
    stop("sigma must be a postive number.")
    if(B<0 | is.na(B))
    stop("B must be a postive integer.")
    
    B = as.integer(B)
    data = as.matrix(data)
    if(sum(is.na(data)))
    stop("Clustering data contains NA or character value. The current version does not support missing data situation.")
    
    ##require("multicore")
    if( nalgorithm ==2) {
        if (mumethods!= 0 || methods >=2) {
            stop("Quadtraic penalty based algorithm cannot deal with the selected objective function. You can try ADMM instead.")
        }
    }


    res = mclapply(1:B,GDF.Step23,data = data,lambda1 = lambda1,lambda2 = lambda2,
    tau = tau, mumethods = mumethods, methods = methods,sigma = sigma,algorithm = nalgorithm,epsilon = epsilon)
    nrows = dim(data)[1]
    ncols = dim(data)[2]
    num = nrows * ncols
    slope = matrix(NA,1,num)
    for(ii in 1:num)
    {
        x = matrix(NA,B,1)
        y = matrix(NA,B,1)
        row = floor((ii-1)/ncols)+1
        col = ii %% ncols
        if (col == 0)
        {
            col = ncols
        }
        for(i in 1:B){
            
            x[i] = res[[i]][[1]][row,col]
            y[i] = res[[i]][[2]][row,col]
        }
        slope[ii] = coefficients(lm(y~x))[2]
    }
    
    GDF = sum(slope)
    
    ## calculate the original mu
    if (nalgorithm == 1){
        rho = lambda1
        a = .Call('prclust_PRclustADMM', PACKAGE = 'prclust', data, rho, lambda2, tau,mumethods, methods,epsilon)
    } else {
        a =    .Call('prclust_PRclustOriginal', PACKAGE = 'prclust', data, lambda1, lambda2, tau, mumethods, methods)
        
    }
    
    GCV = sum((data-a$mu)^2)/(nrows*ncols- GDF)^2
    groupNum = length(unique(a$group))
    
    ## esitmated sigma
    estSigma = sum((data-a$mu)^2)/(nrows*ncols- GDF)
    #  sum((data[1:2,]-a$mu[1:2,])^2)/(2 * ncols- GDF)^2
    out = t(as.matrix(c(GDF,GCV,groupNum,estSigma)))
    colnames(out) = c("GDF","GCV","groupNum","estSigmaSquare")
    out
}