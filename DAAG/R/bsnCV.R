bsnCV <-
function (m = 100, n = 40, method = "exhaustive", nvmax = 3,
              X = NULL, y=NULL, intercept=TRUE, nfolds = 2,
              print.summary = TRUE, really.big = FALSE)
{
    leaps.out <- try(requireNamespace("leaps"), silent = TRUE)
    if (!is.logical(leaps.out) | (leaps.out == FALSE)) {
        print("Error: package leaps is not installed properly")
        return()
    }
        if (is.null(X)) {
            X <- matrix(rnorm(m * n), ncol = n)
            colnames(X) <- paste("V", 1:n, sep = "")
        }
        else {
            if(is.data.frame(X)){
            if(intercept) X <- model.matrix(~., data=X)[,-1] else
            X <- model.matrix(~-1+., data=X)
        }
            m <- dim(X)[1]
            n <- dim(X)[2]
        }
        if (is.null(colnames(X)))
                colnames(X) <- paste("V", 1:n, sep = "")
    if(is.null(y))y <- rnorm(m)
    foldid <- sample(1:nfolds, m, replace = TRUE)
    objlist <- vector("list", length = nfolds)
    for (i in 1:nfolds) {
        train <- foldid != i
        test <- !train
        xxi <- X[train, ]
        yi <- y[train]
        u <- leaps::regsubsets(xxi, yi, method = method, nvmax = nvmax,
                        nbest = 1,  intercept=intercept, really.big = really.big)
        if(intercept){
        x <- X[test, summary(u)$which[nvmax, -1]]
        objlist[[i]] <- lm(y[test] ~ x)} else {
        x <- X[test, summary(u)$which[nvmax, ]]
        objlist[[i]] <- lm(y[test] ~ -1+x)}
    }
    if (print.summary)
        for (i in 1:nfolds) print(summary(objlist[[i]]))
    invisible(objlist)
}
