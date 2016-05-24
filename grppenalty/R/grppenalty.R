## concave 1-norm and 2-norm group penalties
## Dingfeng Jiang
## Version 2.0-0
##dyn.load("../src-x64/grppenalty.dll")

grppenalty <- function(y, x, index, family = "gaussian",
                       type = "l1", penalty = "mcp", kappa = 1/2.7,
                       nlambda = 100, lambda.min = 0.01,
                       epsilon = 1e-3, maxit = 1e+3){
    ## check the consistency of the argument
    if (length(y) != nrow(x)) stop("Dimension of x does not match that of y! \n")
    if (length(index) != ncol(x)) stop("The # of columns in x does not match the group index! \n")
    if (!any(penalty == c("scad", "mcp"))) stop("Please specify 'scad' or 'mcp'! \n")
    if (!any(type == c("l1", "l2"))) stop("Please specify 'l1' or 'l2'! \n")
    if (!any(family == c("gaussian", "binomial"))) stop("Please specify 'gaussian' or 'binomial'! \n")
    if (max(kappa) >= 1.0 | min(kappa) < 0) stop("Please specify a kappa in [0, 1)!\n")
    ## get the dimension
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1, n)
    qp <- qq + p
    ## re-organize x such that variables are grouped
    ordix <- order(index)
    if (is.null(colnames(x))){
        names <- paste("x", 1:p, sep = "")
    } else {
        names <- colnames(x)
    }
    xnames <- c("intercept", names[ordix])
    xx <- x[, ordix]
    nidx <- index[ordix]
    dfgrp <- as.vector(table(nidx))
    ngrp <- length(dfgrp)
    ## working kappa
    if (any(kappa == 0)){
        ukas <- kappa
    }else{
        ukas <- c(0,kappa)
    }
    ookas <- ukas
    if (family == "gaussian"){
        if (penalty == "scad") {
            ukas <- ukas/(1+ukas)
        } else if (penalty == "mcp") {
            ukas <- ukas
        }
    } else if (family == "binomial"){
        if (penalty == "scad") {
            ukas <- ukas/(4+ukas)
        } else if (penalty == "mcp") {
            ukas <- ukas/4
        }
    }
    nkap <- length(ukas)
    ## assign memory for FORTRAN
    olmdas <- rep(0, nkap*nlambda)
    okas <- rep(0, nkap*nlambda)
    ocoef <- matrix(0, qp, nkap*nlambda)
    odf <- rep(0, nkap*nlambda)
    oevidx <- rep(0, nkap*nlambda)
    ## hard work by FORTRAN
    if (family == "gaussian" & type == "l1"){
        if (penalty == "mcp"){
            out <- try(.Fortran("bfl1gmcpga",as.double(olmdas),as.double(okas),
                                as.double(ocoef),as.integer(odf),as.integer(oevidx),
                                as.double(y),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else {
            out <- try(.Fortran("bfl1gscadga",as.double(olmdas),as.double(okas),
                                as.double(ocoef),as.integer(odf),as.integer(oevidx),
                                as.double(y),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else if (family == "binomial" & type == "l1"){
        if (penalty == "mcp"){
            out <- try(.Fortran("bfl1gmcpbi",as.double(olmdas),as.double(okas),
                                as.double(ocoef),as.integer(odf),as.integer(oevidx),
                                as.double(y),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else{
            out <- try(.Fortran("bfl1gscadbi",as.double(olmdas),as.double(okas),
                                as.double(ocoef),as.integer(odf),as.integer(oevidx),
                                as.double(y),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else if (family == "gaussian" & type == "l2"){
        if (penalty == "mcp"){
            out <- try(.Fortran("bfl2gmcpga",as.double(olmdas),as.double(okas),
                                as.double(ocoef),as.integer(odf),as.integer(oevidx),
                                as.double(y),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else {
            out <- try(.Fortran("bfl2gscadga",as.double(olmdas),as.double(okas),
                                as.double(ocoef),as.integer(odf),as.integer(oevidx),
                                as.double(y),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else if (family == "binomial" & type == "l2"){
        if (penalty == "mcp"){
            out <- try(.Fortran("bfl2gmcpbi",as.double(olmdas),as.double(okas),
                                as.double(ocoef),as.integer(odf),as.integer(oevidx),
                                as.double(y),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else{
            out <- try(.Fortran("bfl2gscadbi",as.double(olmdas),as.double(okas),
                                as.double(ocoef),as.integer(odf),as.integer(oevidx),
                                as.double(y),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    }
    ## organize the output
    if (!(inherits(out,"try-error"))) {
        lambda <- unique(out[[1]])
        ocoef <- matrix(out[[3]], nrow=qp)
        coef <- lapply(1:nkap, function(v){
            idx <- (0:(nlambda-1))*nkap + v
            outm <- ocoef[ ,idx]
            rownames(outm) <- xnames
            outm
        })
        ##print(kappa)
        ##print(ookas)
        if (length(kappa) + 1 == length(ookas)){
            kappa <- ookas[-1]
            coef.beta <- coef[-1]
        } else {
            kappa <- ookas
            coef.beta <- coef
        }
    }
    list(coef.beta, kappa, lambda)
}

path.plot <- function(out){
    coef.beta <- out[[1]]
    kappa <- out[[2]]
    nkap <- length(coef.beta)
    qp <- nrow(coef.beta[[1]])
    nlmda <- ncol(coef.beta[[1]])
    ymax <- max(as.vector(unlist(coef.beta)))
    ymin <- min(as.vector(unlist(coef.beta)))
    if (nkap != 1){
        par(mfrow = c(ceiling(nkap/2),2))
    }
    for(i in 1:nkap){
        plot(1:nlmda, coef.beta[[i]][2,], type = "l",
             ylim=c(ymin,ymax), xlab="Grids of lambda",
             ylab="Profile",
             main=paste("kappa=",format(kappa[i],digits=3)))
        for(j in 3:qp){
            lines(1:nlmda, coef.beta[[i]][j,])
        }
    }
}

cv.plot <- function(cv.out){
    cvval <- cv.out[[2]]
    kappa <- unique(as.vector(cv.out[[3]]))
    nkap <- nrow(cvval)
    nlmda <- ncol(cvval)
    ymin <- min(cvval)
    ymax <- max(cvval)
    plot(1:nlmda, cvval[1,], ylim=c(ymin, ymax),
         type="l", xlab="Grids of lambda",
         ylab="CV performance", main = "Cross validation process" )
    if (nkap >= 2){
        for(i in 2:nkap){
            lines(1:nlmda, cvval[i,], col = i)
        }
    }
    legend("bottomright", col=1:nkap, lwd=1.2,
           legend=paste(format(kappa,digits=3)))
}



## tuning parameter selection using cross validation
## use PMSE for linear model, PAUC for logistic model
cv.grppenalty <- function(y, x, index, family = "gaussian",
                          type = "l1", penalty = "mcp", kappa = 1/2.7,
                          nfold = 5, regular.only = TRUE,
                          nlambda = 100, lambda.min = 0.01,
                          epsilon = 1e-3, maxit = 1e+3, seed = 1000){
    ## check the consistency of the argument
    if (length(y) != nrow(x)) stop("Dimension of x does not match that of y! \n")
    if (length(index) != ncol(x)) stop("The # of columns in x does not match the group index! \n")
    if (!any(penalty == c("scad", "mcp"))) stop("Please specify 'scad' or 'mcp'! \n")
    if (!any(type == c("l1", "l2"))) stop("Please specify 'l1' or 'l2'! \n")
    if (!any(family == c("gaussian", "binomial"))) stop("Please specify 'gaussian' or 'binomial'! \n")
    if (max(kappa) >= 1.0 | min(kappa) < 0) stop("Please specify a kappa in [0, 1)!\n")
    ## get the dimension
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1, n)
    qp <- qq + p
    ## working kappa
    if (any(kappa == 0)){
        ukas <- kappa
    }else{
        ukas <- c(0,kappa)
    }
    ookas <- ukas
    if (family == "gaussian"){
        if (penalty == "scad") {
            ukas <- ukas/(1+ukas)
        } else if (penalty == "mcp") {
            ukas <- ukas
        }
    } else if (family == "binomial"){
        if (penalty == "scad") {
            ukas <- ukas/(4+ukas)
        } else if (penalty == "mcp") {
            ukas <- ukas/4
        }
    }
    nkap <- length(ukas)
    ## cross validation index
    cvk <- nfold
    set.seed(seed)
    if (family == "bionomial") {
        ## strafified cross validation for binomial outcome
        n0 <- length(y[y==0])
        n1 <- n-n0
        cvpool0 <- rep(rep(1:cvk), length = n0)
        cvpool1 <- rep(rep(1:cvk), length = n1)
        nidx0 <- sample(cvpool0, replace = FALSE)
        nidx1 <- sample(cvpool1, replace = FALSE)
        nindex <- c(nidx0, nidx1)
        ## reorder the data
        odridx <- order(y)
        oy <- y[odridx]
        ox <- x[odridx,]
    } else{
        ## cross validation without stratification
        cvpool <- rep(rep(1:cvk), length = n)
        nindex <- sample(cvpool, replace = FALSE)
        oy <- y
        ox <- x
    }
    ## re-organize x such that variables are grouped
    ordix <- order(index)
    if (is.null(colnames(x))){
        names <- paste("x", 1:p, sep = "")
    } else {
        names <- colnames(x)
    }
    xnames <- c("intercept", names[ordix])
    xx <- ox[, ordix]
    nidx <- index[ordix]
    dfgrp <- as.vector(table(nidx))
    ngrp <- length(dfgrp)
    ## assign membery for FORTRAN
    oout <- rep(0,3 + qp)
    opmse <- rep(0, nkap*nlambda)
    opauc <- rep(0, nkap*nlambda)
    olmdas <- rep(0, nkap*nlambda)
    okas <- rep(0, nkap*nlambda)
    ocoef <- matrix(0, qp, nkap*nlambda)
    ofull <- rep(0, nkap*nlambda)
    odf <- rep(0, nkap*nlambda)
    ocvx <- rep(0, nkap*nlambda)
    ocvfull <- rep(0, nkap*nlambda)
    ocvcvx <- rep(0, nkap*nlambda)
    ## use fortran to do the computation
    if (family == "gaussian" & type == "l1"){
        if (penalty == "mcp") {
            out <- try(.Fortran("cvl1gmcpga",as.double(oout),as.double(opmse),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.integer(ofull),as.integer(odf),as.integer(ocvx),
                                as.integer(ocvfull),as.integer(ocvcvx),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else {
            out <- try(.Fortran("cvl1gscadga",as.double(oout),as.double(opmse),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.integer(ofull),as.integer(odf),as.integer(ocvx),
                                as.integer(ocvfull),as.integer(ocvcvx),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else if (family == "binomial" & type =="l1") {
        if (penalty == "mcp"){
            out <- try(.Fortran("cvl1gmcpbi",as.double(oout),as.double(opmse),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.integer(ofull),as.integer(odf),as.integer(ocvx),
                                as.integer(ocvfull),as.integer(ocvcvx),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else {
            out <- try(.Fortran("cvl1gscadbi",as.double(oout),as.double(opmse),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.integer(ofull),as.integer(odf),as.integer(ocvx),
                                as.integer(ocvfull),as.integer(ocvcvx),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else if (family == "gaussian" & type == "l2"){
        if (penalty == "mcp"){
            out <- try(.Fortran("cvl2gmcpga",as.double(oout),as.double(opmse),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.integer(ofull),as.integer(odf),as.integer(ocvx),
                                as.integer(ocvfull),as.integer(ocvcvx),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else {
            out <- try(.Fortran("cvl2gscadga",as.double(oout),as.double(opmse),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.integer(ofull),as.integer(odf),as.integer(ocvx),
                                as.integer(ocvfull),as.integer(ocvcvx),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else if (family == "binomial" & type == "l2") {
        if (penalty == "mcp"){
            out <- try(.Fortran("cvl2gmcpbi",as.double(oout),as.double(opmse),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.integer(ofull),as.integer(odf),as.integer(ocvx),
                                as.integer(ocvfull),as.integer(ocvcvx),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else {
            out <- try(.Fortran("cvl2gscadbi",as.double(oout),as.double(opmse),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.integer(ofull),as.integer(odf),as.integer(ocvx),
                                as.integer(ocvfull),as.integer(ocvcvx),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy),as.double(int),as.double(xx),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(ngrp),as.integer(dfgrp),
                                as.integer(nkap),as.double(ukas),
                                as.integer(nlambda),as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    }
    ## re-organize output
    if (!(inherits(out,"try-error"))) {
        cvval <- out[[2]]
        olmdas <- out[[3]]
        okas <- rep(ookas, nlambda)
        ocoef <- matrix(out[[5]], nrow=qp)
        ofull <- out[[6]]
        odf <- out[[7]]
        if (length(kappa) + 1 == length(ookas)){
            sidx <- (0:(nlambda-1))*nkap +1
            cvval <- cvval[-sidx]
            olmdas <- olmdas[-sidx]
            okas <- okas[-sidx]
            ocoef <- ocoef[,-sidx]
            ofull <- ofull[-sidx]
            odf <- odf[-sidx]
        }
    }
    ## get the regular solution with ms<n
    if (regular.only){
        useidx <- ifelse(ofull == 0, TRUE, FALSE)
        ucv <- cvval[useidx]
        ulmda <- olmdas[useidx]
        uka <- okas[useidx]
        ucoef <- ocoef[, useidx]
    } else{
        ucv <- cvval
        ulmda <- olmdas
        uka <- okas
        ucoef <- ocoef
    }
    ## find the minimum for gaussian and maximum for binomial
    if (family == "gaussian"){
        idx <- ucv == min(ucv)
    } else {
        idx <- ucv == max(ucv)
    }
    ## output
    cv.value <- ucv[idx][1]
    lambda <- ulmda[idx][1]
    kappa <- uka[idx][1]
    coefficients <- ucoef[, idx]
    if (is.matrix(coefficients)) coefficients <- coefficients[, 1]
    names(coefficients) <- xnames
    selected <- list(cv.value, coefficients, kappa, lambda)
    cv.values <- matrix(cvval, ncol=nlambda)
    kappas <- matrix(okas, ncol=nlambda)
    lambdas <- matrix(olmdas, ncol=nlambda)
    list(selected, cv.values, kappas, lambdas)
}


