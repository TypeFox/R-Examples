cox.ridge <-
function(formula = formula(data),lambda=1, lambdaFixed=FALSE, eps=10e-06, data = sys.parent(),iter.max=200,mon=FALSE){
    Sumj <- function(x) { sum(x)-cumsum(x)+x}
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data")
    m <- m[match(temp, names(m), nomatch = 0)]
    Terms <- terms(formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    Y <- model.extract(m, "response")
    nm <- attr(Terms, "term.labels")
    if(!inherits(Y, "Surv"))
        stop("Response must be a survival object")
    X <- model.matrix(Terms, m)[, -1, drop = F]
    X <- apply(X,2,function(x) {x-mean(x)})
    death <- Y[,2]
    n <-nrow(X); p <- ncol(X); b <-rep(0,p);Ip <- diag(rep(1,p))
    it <- 0; itl <- 0; dl <- 5
    while(dl>eps & itl <iter.max){   
        dz <- 5
        while (dz> eps & it < iter.max){
            H0 <- cumsum(death/(Sumj(exp(X %*% b))))
            h0 <- 1/(Sumj(exp(X %*% b)))
            lpen <- sum((X%*%b)*death - death*log(Sumj(exp(X%*%b))))-0.5 * lambda * t(b) %*% b
            Diaf <- death- H0 * exp(X %*%b) 
            sc <- t(X) %*% Diaf-lambda*b
            Dm <- diag(as.vector(H0 * exp(X %*%b)))
            Info <-  -(t(X) %*% Dm %*% X + lambda* Ip)
            sInfo <- solve(Info)
            b <- b-sInfo%*%sc
            Hat <- t(X) %*% Dm %*% X %*% sInfo
            H0 <- cumsum(death/(Sumj(exp(X %*% b))))   
            h0 <- 1/(Sumj(exp(X %*% b)))
            lpenn <- sum((X%*%b)*death - death*log(Sumj(exp(X%*%b))))-0.5 * lambda * t(b) %*% b
            dz <- abs(lpen -lpenn)
            it <- it+1}                      
        lambda2 <- -as.numeric(sum(diag(Hat))/(t(b)%*%b))
        if(lambda2>1e+6){lambda2=1e+06}
        dl <- abs(lambda-lambda2)
        if(lambdaFixed==TRUE){dl= eps}
        if(lambdaFixed==FALSE){lambda <- lambda2}
        itl <- itl +1
        if(mon==TRUE){cat("iteration",itl,". Penalty weight=",lambda,"\n")}
    }
    lpenn <- sum((X%*%b)*death - death*log(Sumj(exp(X%*%b))))-0.5 * lambda * t(b) %*% b 
    fit <- list(call = call, coef=b, loglik=lpen, time = Y[,1],death=Y[,2], X=X,iter =it,inter.it =itl, lambda=lambda, Hat=hat, hess=Info)
    class(fit) <- "cox.ridge"
    fit}
