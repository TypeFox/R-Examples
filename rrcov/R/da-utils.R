.AER <- function(tab)
{
    1 - sum(tab[row(tab) == col(tab)])/sum(tab)
}

mtxconfusion <- function(actual, predicted, prior = NULL, printit=FALSE) {
    .confusion(actual, predicted, prior, printit)
}

.confusion <- function(actual, predicted, prior = NULL, printit=FALSE) {

    names <- levels(actual)
    tab <- table(actual, predicted)
    acctab <- t(apply(tab, 1, function(x) x/sum(x)))
    dimnames(acctab) <- list(Actual = names, "Predicted" = names)
    dimnames(tab) <- list(Actual = names, "Predicted" = names)

    if(is.null(prior))
    {
        cnt <- table(actual)
        prior <- cnt/sum(cnt)
    }
    else
        names(prior) <- names

    AER <- 1 - sum(tab[row(tab) == col(tab)])/sum(tab)

    if(printit)
    {
        prt <- as.matrix(round(c("Apparent error rate" = AER, "Prior frequency" = prior),4))
        colnames(prt) <- ""
        print(prt)

        cat("\nClassification table", "\n")
        print(tab)
        cat("\nConfusion matrix", "\n")
        print(round(acctab, 3))
    }

    invisible(tab)
}

## Internal function to perform leaving-one-out cross validation by brute force -
##  recalculates the estimator n times, excluding each observation in turn.
##
##  - The discriminantion method (QDA or LDA) is selected according to the type of
##      the object.
##  - The method for computing the covariance matrices (or the common
##      cov matrix in LDA) is selected according the slot methodID.
##
.CV <- function(obj){

    if(!is(obj, "Lda") && !is(obj, "Qda"))
        stop("The object must be an Lda or Qda object")

    classic <- is(obj, "LdaClassic") || is(obj, "QdaClassic")
    ret <- predict(obj)

    X <- obj@X
    grp <- obj@grp
    ng <- length(levels(grp))
    method <- obj@method

    ptm <- proc.time()

    n <- nrow(X)
    p <- ncol(X)

    if(!classic && n*p > 500)
        warning("This could take some time!")

    for(i in 1:n)
    {
        cat("i=",i,"\n")

        ll <- if(is(obj, "LdaClassic")) {
                LdaClassic(X[-i,], grouping=grp[-i])
            } else if(is(obj, "Linda")){
                Linda(X[-i,], grouping=grp[-i], method=method)
            } else if(is(obj, "QdaClassic")){
                QdaClassic(X[-i,], grouping=grp[-i])
            } else if(is(obj, "QdaCov")){
                QdaCov(X[-i,], grouping=grp[-i], method=obj@control)
            } else {
                stop("ERROR: unknown class")
            }

        pp <- predict(ll, newdata=t(X[i,]))

        ret@classification[i] <- pp@classification[1]
        ret@posterior[i,] <- pp@posterior[1,]
    }

    ret@ct <- mtxconfusion(grp, ret@classification)

##    cat("\nElapsed time (loo): ",(proc.time() - ptm)[1],"\n")
    ret
}
