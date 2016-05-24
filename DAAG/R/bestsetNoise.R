bestsetNoise <-
function (m = 100, n = 40, method = "exhaustive", nvmax = 3,
              X = NULL, y=NULL, intercept=TRUE,
              print.summary = TRUE, really.big = FALSE, ...)
{
    leaps.out <- try(requireNamespace("leaps"), silent = TRUE)
    if ((is.logical(leaps.out) == TRUE) & (leaps.out == TRUE)) {
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
        u <- leaps::regsubsets(X, y, method = method, nvmax = nvmax,
                        nbest = 1, intercept=intercept, really.big = really.big,
                        ...)
        if(is.null(intercept))intercept <- TRUE
        if(intercept){
        x <- X[, summary(u)$which[nvmax, -1]]
        u1 <- lm(y ~ x)} else {
        x <- X[, summary(u)$which[nvmax, ]]
        u1 <- lm(y ~ -1+x)}
        if (print.summary)
            print(summary(u1, corr = FALSE))
        invisible(list(best=u1, regsubsets_obj=u))
    }
    else {
        print("Error: package leaps is not installed properly")
    }
}
