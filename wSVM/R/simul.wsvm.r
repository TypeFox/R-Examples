simul.wsvm <- function(set.seeds = 123, obs.num = 1000, mu = c(0,0), cov.mat = 0.2*diag(2), X = NULL, Y = NULL){
# generate simulation data set 
    require(MASS)
    set.seed(set.seeds)
    # generate train data set. 
    # We borrow the train from mixture example by Elements of Statistical Learning(Freidman et al. 2000)
    if(is.null(X)) X <- mixture.example$x 
    if(is.null(Y)) Y <- mixture.example$y

    # For SVM, reponse value 0 replace 1
    Y <- ifelse(Y == 0, -1 , 1)
    N <- length(Y)
    
    #plot(x[,2] ~ x[,1], col = ifelse(y == 1,'black','green'))
    
    centers <- c(sample(1:10, obs.num / 2, replace=TRUE), 
                 sample(11:20, obs.num / 2, replace=TRUE))
    means <- mixture.example$means
    means <- means[centers, ]
    mix.test <- mvrnorm(obs.num, mu, cov.mat)
    mix.test <- mix.test + means
    cltest <- c(rep(-1, obs.num/2), rep(1, obs.num/2))
    
    new.X <- as.data.frame(mix.test)
    new.Y <- as.data.frame(cltest)
    colnames(new.Y)[1] <- "z"
    colnames(new.X)[1] <- "x"
    colnames(new.X)[2] <- "y"
    res <- list(X = X, Y = Y, new.X = new.X, new.Y = new.Y)
    return(res)
}
