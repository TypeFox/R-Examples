bridge.EM.R = function(y, X, alpha, ratio=1.0, lambda.max=1e9*ratio, tol=1e-9, max.iter=30, init=NULL){
    X <- as.matrix(X)
    xx <- t(X)%*%X
    xy <- t(X)%*%y
    ixx <- chol2inv(chol(xx))

    p <- ncol(X)

    bhat <- drop(ixx%*%xy)
    Beta = bhat
    if (!is.null(init)) Beta = init

    diff = 1

    # tau = (Nu)^{-1/alpha}
    sigma = 1;
    tau = ratio;

    iter = 0

    while(diff > tol && iter < max.iter)
    {
        # YHat = X %*% Beta
        # sigma = sqrt(sum( (Y-YHat)^2 )/(n-p))
        # sigma=1
        # EXPECTATION STEP
        Lambda = pmin( alpha*(tau^(2-alpha)) * abs(Beta)^(alpha-2), lambda.max)
        #OmegaInv = as.numeric((d+1)/(d*sigma^2+(Y-YHat)^2))
        # H = solve((1/tau^2)*diag(as.numeric(Lambda))+t(X) %*% X) %*% t(X)
        # BetaNew = H %*% y
        BetaNew = solve((1/tau^2)*diag(as.numeric(Lambda))+t(X) %*% X, xy);
        # S = X %*% H
        diff = sum(abs(Beta - BetaNew))
        Beta = BetaNew
        # print(Beta);
        #Nu = (b.nu + sum(abs(Beta)/sigma))/(p + a.nu - 1)
        iter = iter + 1;
    }

    out = list("beta"=drop(Beta), "iter"=iter, "diff"=diff)
    out
}
