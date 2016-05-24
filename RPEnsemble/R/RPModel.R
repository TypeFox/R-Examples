RPModel <-
function(Model.No #Model number
                , n # sample size
                , p # dimension
                , Pi = 1/2 # class 1 prior
                    )
{
    if (Model.No == 1)
        {
            Y1 <- rmultinom(1,n,c(Pi,1-Pi))
            Y <- c(rep(1,Y1[1,1]),rep(2,Y1[2,1]))
            mu <- rep(1/8,p)
            U <- runif(Y1[1,1]*p)
            X1 <- matrix(log(2*U)*(U < 1/2) - log(2-2*U) * (U >= 1/2),Y1[1,1],p)
            X2 <- mvrnorm(Y1[2,1],mu,diag(p))
            X  <- rbind(X1,X2)
        }
    if (Model.No == 2)
        {
            if (p <= 5) stop("model 2 requires p > 5")
            Y1 <- rmultinom(1,n,c(Pi,1-Pi))
            Y <- c(rep(1,Y1[1,1]),rep(2,Y1[2,1]))
            mu <- c(rep(2,5),rep(0,p-5))
            U1 <- rchisq(Y1[1,1],1)
            U2 <- rchisq(Y1[2,1],2)
            Sigma1 <- diag(p)
            Sigma2 <- 0.5*diag(p)+0.5*c(rep(1,5),rep(0,p-5))%*%t(c(rep(1,5),rep(0,p-5))) + 0.5*diag(c(rep(0,5),rep(1,p-5)))
            X1 <- mvrnorm(Y1[1,1],rep(0,p),Sigma1)/sqrt(U1/1)
            X2 <- t(mu + t(mvrnorm(Y1[2,1],rep(0,p),Sigma2)/sqrt(U2/2)))
            X  <- rbind(X1,X2)
        }
    if (Model.No == 3)
        {
        if (p <= 5) stop("model 3 requires p > 5")
            Y1 <- rmultinom(1,n,c(Pi,1-Pi))
            Y <- c(rep(1,Y1[1,1]),rep(2,Y1[2,1]))
            Y11 <- rmultinom(1,Y1[1,1],c(1/2,1/2))
            mu <- c(rep(1,5), rep(0,p-5))
            Sigma <- diag(p)
            X1 <- rbind(t(matrix(mu/2,p,Y11[1,1])),t(matrix(mu/2,p,Y11[2,1]))) + mvrnorm(Y1[1,1],rep(0,p),Sigma)
            X2 <- cbind(matrix(rcauchy(Y1[2,1]*5),Y1[2,1],5), matrix(rnorm(Y1[2,1]*(p-5),0,1),Y1[2,1],p-5))
            X <- rbind(X1,X2)
        }
    if (Model.No == 4)
        {
    if (p != 50) stop("model 4 requires p = 50")
    R <- NULL
    load("R.RData")
    Y1 <- rmultinom(1,n,c(Pi,1-Pi))
    Y <- c(rep(1,Y1[1,1]),rep(2,Y1[2,1]))
    mu <- c(rep(1,3), rep(0,p-3))
    Sigma1 <- 0.5*diag(c(rep(1,3),rep(0,p-3)))+0.5*c(rep(1,3),rep(0,p-3))%*%t(c(rep(1,3),rep(0,p-3))) + 0.5*diag(c(rep(0,3),rep(1,p-3)))+0.5*c(rep(0,3),rep(1,p-3))%*%t(c(rep(0,3),rep(1,p-3)))
    Sigma2 <- 1.5*diag(c(rep(1,3),rep(0,p-3)))+0.5*c(rep(1,3),rep(0,p-3))%*%t(c(rep(1,3),rep(0,p-3))) + 0.5*diag(c(rep(0,3),rep(1,p-3)))+0.5*c(rep(0,3),rep(1,p-3))%*%t(c(rep(0,3),rep(1,p-3)))
    X1 <- mvrnorm(Y1[1,1],R%*%rep(0,p),R%*%Sigma1%*%t(R))
    X2 <- mvrnorm(Y1[2,1],R%*%mu,R%*%Sigma2%*%t(R))
    X <- rbind(X1,X2)
}
    return(list(x=X,y=Y))
}
