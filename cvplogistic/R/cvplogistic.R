## Compute solution surface for a concave penalized logistic regression model
## using majorization minimization by coordinate descent (MMCD) algorithm.
## Two concave penalties are considered: SCAD and MCP.
## For MCP, the local linear approximation (LLA-CD) and adaptive rescaling
## algorithms are also implemented.
## For all the algorithms, the solution surface is computed along kappa.
## Lasso-MCP hybrid penalty is also provided.
## Jan 3, 2012
## dyn.load("../src/cvplogistic.dll")


cvplogistic <- function(y, x, penalty = "mcp", approach = "mmcd",
                        kappa = 1/2.7, nlambda = 100, lambda.min = 0.01,
                        epsilon = 1e-3, maxit = 1e+3){
    ## error checking
    if (nrow(x) != length(y)) stop("# of rows in X does not match the length of Y! \n")
    if ( kappa >= 1 | kappa < 0) {
        stop("Regulation parameter kappa should be in [0, 1.0)!\n")
    }
    ## penalty
    pen <- pmatch(penalty, c("mcp", "scad"))
    app <- pmatch(approach, c("mmcd", "adaptive", "llacda"))
    ## space assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1, n)
    qp <- qq+p
    ## kappa used in the computation
    if (pen == 1) {
        if (app == 1) {
            maxkappa <- kappa/4
        } else {
            maxkappa <- kappa
        }
    } else if ( pen == 2){
        maxkappa <- kappa/(kappa + 4)
    }
    nkappa <- ifelse(maxkappa == 0, 1, 2)
    ## create output space
    olmdas <- rep(0, nkappa*nlambda)
    okas <- rep(0, nkappa*nlambda)
    ocoef <- matrix(0, qp, nkappa*nlambda)
    oaic <- rep(0, nkappa*nlambda)
    obic <- rep(0, nkappa*nlambda)
    oobj <- rep(0, nkappa*nlambda)
    odf <- rep(0, nkappa*nlambda)
    ocvx <- rep(0, nkappa*nlambda)
    ## fit through Fortran depend on pen, app and ss
    if (pen == 1){
        if (app == 1){
            out <- try(.Fortran("mcpkapa",
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx),
                                as.double(y), as.double(int), as.double(x),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(lambda.min),
                                as.double(epsilon), as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        } else if (app == 2) {
            out <- try(.Fortran("adpmcpkp",
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx),
                                as.double(y), as.double(int), as.double(x),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        } else if (app == 3) {
            out <- try(.Fortran("fllabi",
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx),
                                as.double(y), as.double(int), as.double(x),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        }
    } else if (pen == 2) {
        out <- try(.Fortran("scadkapa",
                            as.double(olmdas), as.double(okas), as.double(ocoef),
                            as.double(oaic), as.double(obic), as.double(oobj),
                            as.integer(odf), as.integer(ocvx),
                            as.double(y), as.double(int), as.double(x),
                            as.integer(n), as.integer(qq), as.integer(p),
                            as.integer(nkappa), as.double(maxkappa),
                            as.integer(nlambda), as.double(lambda.min),
                            as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
    }
    ## organize output
    if (!(inherits(out, 'try-error'))){
        lambdas <- unique(out[[1]])
        coef <- matrix(out[[3]], qp, nkappa*nlambda)[,(1:nlambda)*nkappa]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept", paste("x", 1:p, sep = ""))
        } else {
            rownames(coef) <- c("intercept", colnames(x))
        }
    }
    list(coef, lambdas)
}

hybrid.logistic <- function(y, x, penalty = "mcp",
                            kappa = 1/2.7, nlambda = 100, lambda.min = 0.01,
                            epsilon = 1e-3, maxit = 1e+3){
    ## error checking
    if (nrow(x) != length(y)) stop("# of rows in X does not match the length of Y! \n")
    if ( kappa >= 1 | kappa < 0) {
        stop("Regulation parameter kappa should be in [0, 1.0)!\n")
    }
    ## penalty
    pen <- pmatch(penalty, c("mcp", "scad"))
    ## space assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1, n)
    qp <- qq+p
    ## kappa used in the computation
    if (pen == 1) {
        maxkappa <- kappa/4
    } else if ( pen == 2){
        maxkappa <- kappa/(kappa + 4)
    }
    nkappa <- ifelse(maxkappa == 0, 1, 2)
    ## create output space
    olmdas <- rep(0, nkappa*nlambda)
    okas <- rep(0, nkappa*nlambda)
    ocoef <- matrix(0, qp, nkappa*nlambda)
    oaic <- rep(0, nkappa*nlambda)
    obic <- rep(0, nkappa*nlambda)
    oobj <- rep(0, nkappa*nlambda)
    odf <- rep(0, nkappa*nlambda)
    ocvx <- rep(0, nkappa*nlambda)
    ## fit through Fortran depend on pen, app and ss
    if (pen == 1){
        out <- try(.Fortran("mcpkapa2",
                            as.double(olmdas), as.double(okas), as.double(ocoef),
                            as.double(oaic), as.double(obic), as.double(oobj),
                            as.integer(odf), as.integer(ocvx),
                            as.double(y), as.double(int), as.double(x),
                            as.integer(n), as.integer(qq), as.integer(p),
                            as.integer(nkappa), as.double(maxkappa),
                            as.integer(nlambda), as.double(lambda.min),
                            as.double(epsilon), as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
    } else if (pen == 2) {
        out <- try(.Fortran("scadkapa2",
                            as.double(olmdas), as.double(okas), as.double(ocoef),
                            as.double(oaic), as.double(obic), as.double(oobj),
                            as.integer(odf), as.integer(ocvx),
                            as.double(y), as.double(int), as.double(x),
                            as.integer(n), as.integer(qq), as.integer(p),
                            as.integer(nkappa), as.double(maxkappa),
                            as.integer(nlambda), as.double(lambda.min),
                            as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
    }
    ## organize output
    if (!(inherits(out, 'try-error'))){
        lambdas <- unique(out[[1]])
        coef <- matrix(out[[3]], qp, nkappa*nlambda)[,(1:nlambda)*nkappa]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept", paste("x", 1:p, sep = ""))
        } else {
            rownames(coef) <- c("intercept", colnames(x))
        }
    }
    list(coef, lambdas)
}





## Tuning parameter selection using CV-AUC
cv.cvplogistic <- function(y, x, penalty = "mcp", approach = "mmcd",
                           nfold = 5, kappa = 1/2.7,
                           nlambda = 100, lambda.min = 0.01,
                           epsilon = 1e-3, maxit = 1e+3, seed = 1000){
    ## error checking
    if (nrow(x) != length(y)) stop("# of rows in X does not match the length of Y! \n")
    if ( kappa >= 1 | kappa < 0) {
        stop("Regulation parameter kappa should be in [0, 1.0)!\n")
    }
    ## penalty
    pen <- pmatch(penalty, c("mcp", "scad"))
    app <- pmatch(approach, c("mmcd", "adaptive", "llacda"))
    ## space assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1, n)
    qp <- qq+p
    ## kappa used in the computation
    if (pen == 1) {
        if (app == 1) {
            maxkappa <- kappa/4
        } else {
            maxkappa <- kappa
        }
    } else if ( pen == 2){
        maxkappa <- kappa/(kappa + 4)
    }
    nkappa <- ifelse(maxkappa == 0, 1, 2)
    ## cross validation index
    cvk <- nfold
    set.seed(seed)
    ## strafified cross validation only
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
    ## create output space
    oout <- rep(0, 3+qp)
    opauc <- rep(0, nkappa*nlambda)
    olmdas <- rep(0, nkappa*nlambda)
    okas <- rep(0, nkappa*nlambda)
    ocoef <- matrix(0, qp, nkappa*nlambda)
    oaic <- rep(0, nkappa*nlambda)
    obic <- rep(0, nkappa*nlambda)
    oobj <- rep(0, nkappa*nlambda)
    odf <- rep(0, nkappa*nlambda)
    ocvx <- rep(0, nkappa*nlambda)
    ofull <- rep(0, nkappa*nlambda)
    cvcvx <- rep(0, nkappa*nlambda)
    cvfull <- rep(0, nkappa*nlambda)
    ## cross validation process by Fortran
    if (pen == 1){
        if (app == 1){
            out <- try(.Fortran("cvauckapa", as.double(oout), as.double(opauc),
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                as.integer(cvcvx), as.integer(cvfull),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy), as.double(int), as.double(ox),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        } else if (app == 2){
            out <- try(.Fortran("adpcvauckp", as.double(oout), as.double(opauc),
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                as.integer(cvcvx), as.integer(cvfull),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy), as.double(int), as.double(ox),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        } else if (app == 3){
            out <- try(.Fortran("fcvllabi", as.double(oout), as.double(opauc),
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                as.integer(cvcvx), as.integer(cvfull),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy), as.double(int), as.double(ox),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(lambda.min),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        }
    } else if (pen == 2) {
        out <- try(.Fortran("cvaucsdka", as.double(oout), as.double(opauc),
                            as.double(olmdas), as.double(okas), as.double(ocoef),
                            as.double(oaic), as.double(obic), as.double(oobj),
                            as.integer(odf), as.integer(ocvx), as.integer(ofull),
                            as.integer(cvcvx), as.integer(cvfull),
                            as.integer(nindex), as.integer(cvk),
                            as.double(oy), as.double(int), as.double(ox),
                            as.integer(n), as.integer(qq), as.integer(p),
                            as.integer(nkappa), as.double(maxkappa),
                            as.integer(nlambda), as.double(lambda.min),
                            as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
    }
    ## selection based on CV-AUC,  df
    idx <- (1:nlambda)*nkappa
    if (!(inherits(out, 'try-error'))){
        cvauc <- out[[2]][idx]
        lambda <- out[[3]][idx]
        kappa <- out[[4]][idx]
        coef <- matrix(out[[5]], qp, nkappa*nlambda)[,idx]
        df <- out[[9]][idx]
        fullmodel <- out[[11]][idx]
        cvfullmodel <- out[[13]][idx]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept", paste("x", 1:p, sep = ""))
        } else {
            rownames(coef) <- c("intercept", colnames(x))
        }
        ## regular solution if n > df
        uidx <- n > df
        ucvauc <- cvauc[uidx]
        ulambda <- lambda[uidx]
        ukappa <- kappa[uidx]
        ucoef <- coef[, uidx]
        ## maximum cvauc
        aucidx <- ucvauc ==  max(ucvauc)
        scvauc <- ucvauc[aucidx][[1]]
        slambda <- ulambda[aucidx][[1]]
        skappa <- ukappa[aucidx][[1]]
        scoef <- ucoef[, aucidx]
        if (is.matrix(scoef)) scoef <- scoef[, 1]
        list(scvauc, slambda, scoef)
    } else {
         stop("Model fitting fails, double check!\n")
         NULL
    }
}


## Tuning parameter selection using CV-AUC
cv.hybrid <- function(y, x, penalty = "mcp", nfold = 5, kappa = 1/2.7,
                      nlambda = 100, lambda.min = 0.01,
                      epsilon = 1e-3, maxit = 1e+3, seed = 1000){
    ## error checking
    if (nrow(x) != length(y)) stop("# of rows in X does not match the length of Y! \n")
    if ( kappa >= 1 | kappa < 0) {
        stop("Regulation parameter kappa should be in [0, 1.0)!\n")
    }
    ## penalty
    pen <- pmatch(penalty, c("mcp", "scad"))
    ## space assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1, n)
    qp <- qq+p
    ## kappa used in the computation
    if (pen == 1) {
        maxkappa <- kappa/4
    } else if ( pen == 2){
        maxkappa <- kappa/(kappa + 4)
    }
    nkappa <- ifelse(maxkappa == 0, 1, 2)
    ## cross validation index
    cvk <- nfold
    set.seed(seed)
    ## strafified cross validation only
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
    ## create output space
    oout <- rep(0, 3+qp)
    opauc <- rep(0, nkappa*nlambda)
    olmdas <- rep(0, nkappa*nlambda)
    okas <- rep(0, nkappa*nlambda)
    ocoef <- matrix(0, qp, nkappa*nlambda)
    oaic <- rep(0, nkappa*nlambda)
    obic <- rep(0, nkappa*nlambda)
    oobj <- rep(0, nkappa*nlambda)
    odf <- rep(0, nkappa*nlambda)
    ocvx <- rep(0, nkappa*nlambda)
    ofull <- rep(0, nkappa*nlambda)
    cvcvx <- rep(0, nkappa*nlambda)
    cvfull <- rep(0, nkappa*nlambda)
    ## cross validation process by Fortran
    if (pen == 1){
        out <- try(.Fortran("cvauckapa2", as.double(oout), as.double(opauc),
                            as.double(olmdas), as.double(okas), as.double(ocoef),
                            as.double(oaic), as.double(obic), as.double(oobj),
                            as.integer(odf), as.integer(ocvx), as.integer(ofull),
                            as.integer(cvcvx), as.integer(cvfull),
                            as.integer(nindex), as.integer(cvk),
                            as.double(oy), as.double(int), as.double(ox),
                            as.integer(n), as.integer(qq), as.integer(p),
                            as.integer(nkappa), as.double(maxkappa),
                            as.integer(nlambda), as.double(lambda.min),
                            as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
    } else if (pen == 2) {
        out <- try(.Fortran("cvaucsdka2", as.double(oout), as.double(opauc),
                            as.double(olmdas), as.double(okas), as.double(ocoef),
                            as.double(oaic), as.double(obic), as.double(oobj),
                            as.integer(odf), as.integer(ocvx), as.integer(ofull),
                            as.integer(cvcvx), as.integer(cvfull),
                            as.integer(nindex), as.integer(cvk),
                            as.double(oy), as.double(int), as.double(ox),
                            as.integer(n), as.integer(qq), as.integer(p),
                            as.integer(nkappa), as.double(maxkappa),
                            as.integer(nlambda), as.double(lambda.min),
                            as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
    }
    ## selection based on CV-AUC,  df
    idx <- (1:nlambda)*nkappa
    if (!(inherits(out, 'try-error'))){
        cvauc <- out[[2]][idx]
        lambda <- out[[3]][idx]
        kappa <- out[[4]][idx]
        coef <- matrix(out[[5]], qp, nkappa*nlambda)[,idx]
        df <- out[[9]][idx]
        fullmodel <- out[[11]][idx]
        cvfullmodel <- out[[13]][idx]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept", paste("x", 1:p, sep = ""))
        } else {
            rownames(coef) <- c("intercept", colnames(x))
        }
        ## regular solution if n > df
        uidx <- n > df
        ucvauc <- cvauc[uidx]
        ulambda <- lambda[uidx]
        ukappa <- kappa[uidx]
        ucoef <- coef[, uidx]
        ## maximum cvauc
        aucidx <- ucvauc ==  max(ucvauc)
        scvauc <- ucvauc[aucidx][[1]]
        slambda <- ulambda[aucidx][[1]]
        skappa <- ukappa[aucidx][[1]]
        scoef <- ucoef[, aucidx]
        if (is.matrix(scoef)) scoef <- scoef[, 1]
        list(scvauc, slambda, scoef)
    } else {
         stop("Model fitting fails, double check!\n")
         NULL
    }
}

path.plot <- function(out){
    coef <- out[[1]]
    nlambda <- ncol(coef)
    p <- nrow(coef)-1
    ymin <- min(coef[-1,])
    ymax <- max(coef[-1,])
    plot(1:nlambda, coef[2, ], ylim=c(ymin, ymax), type = "l",
         xlab = "Grids of lambda", ylab = "Coefficient profile")
    for(i in 3:p){
        lines(1:nlambda, coef[i,])
    }
}
