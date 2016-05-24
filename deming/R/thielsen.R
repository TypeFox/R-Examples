# Function for Theil-Sen regression
thielsen <- function(formula, data, subset, weights, na.action, conf=.95,
                  nboot=0, symmetric=FALSE, 
                  eps=sqrt(.Machine$double.eps),
                  x=FALSE, y=FALSE, model=TRUE) {
    Call <- match.call()

    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    # then evaluate it in the proper frame
    indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp[[1]] <- as.name('model.frame')  # change the function called
    mf <- eval(temp, parent.frame())
    Terms <- terms(mf)

    X <- model.matrix(Terms, mf)
    if (!attr(Terms, "intercept")) stop ("an intercept is required")
    if (ncol(X) != 2) 
        stop("Thiel-Sen regression requires a single predictor variable")

    Y <- model.response(mf, type="numeric")
    n <- length(Y)
    if (is.null(Y))
        stop ("a response variable is required")
    if (n <3) stop("less than 3 non-missing observations in the data")

    wt <- model.weights(mf)
    if (is.null(wt)) wt <- rep(1.0, n)
    else {
        if (conf >0 && nboot==0) 
            stop("confidence intervals for weighted data require a bootstrap")
    }

    if (!is.numeric(eps) || eps <=0) stop("invalid value for eps")
    if (!is.logical(symmetric)) stop("invalid value for symmetric argument")
    if (conf <0 || conf >=1) stop("invalid confidence interval limit")
    
    fit <- theilsen.fit(X[,2], Y, wt, symmetric, conf*(nboot==0), eps)
    names(fit$coefficients) <-  dimnames(X)[[2]]
    yhat <- fit$coefficients[1] + X[,2]*fit$coefficients[2]
    fit$residuals <- Y-yhat

    if (nboot > 0) {
        # Set up for the wild bootstrap, which add random noise to
        #  the original x,y values
        fcoef <- fit$coefficients
        resid <- list(y= Y - yhat, x = X[,2] - (Y-fcoef[1])/fcoef[2])
        wild <- function(data, resid) {
            n <- length(data$x)
            rb <- matrix(rbinom(2*n, 1, prob=(sqrt(5)+1)/sqrt(20)), ncol=2)
            newx <- data$x + resid$x*ifelse(rb[,1]==1, -(sqrt(5)-1)/2,
                                             (sqrt(5)+1)/2)
            newy <- data$y + resid$y*ifelse(rb[,2]==1, -(sqrt(5)-1)/2,
                                             (sqrt(5)+1)/2)
            list(x=newx, y=newy, w=data$w)
        }
        dfun <- function(data) {  #will be called by boot function
            theilsen.fit(data$x, data$y, data$w, symmetric, conf=0, 
                           eps)$coefficients
         }
        boot.out <- boot(list(x=X[,2], y=Y, w= wt), dfun, nboot, 
                         sim="parametric",
                         ran.gen=wild, mle=resid)
        temp1 <- boot.ci(boot.out, type="perc", index=1, conf=conf)
        temp2 <- boot.ci(boot.out, type="perc", index=2, conf=conf)
        fit$boot <-boot.out
        fit$ci   <- rbind(temp1$perc[,4:5], temp2$perc[,4:5])
        fit$var <- var(boot.out$t)
    }

 #   dimnames(fit$ci) <- list(names(fit$coef), 
#                             paste(c("lower", "upper"), format(conf)))
    fit$n <- n
    fit$conf <- conf
    if (x) fit$x <- X
    if (y) fit$y <- Y
    if (model) fit$model <- mf
    
    na.action <- attr(mf, "na.action")
    if (length(na.action)) fit$na.action <- na.action
    fit$terms <- Terms
    class(fit) <- c("thielsen") 
    fit$call <- Call
    fit
}


print.thielsen <- function(x, ...) {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("n=", x$n)
    if (length(x$na.action))
        cat("  (", naprint(x$na.action), ")\n", sep='')
    else cat("\n")

    temp <- matrix(x$coefficients, 2,1)
    tname <- "Coefficient"
    if (!is.null(x$var)) {
        temp <- cbind(temp, sqrt(diag(x$var)))
        tname <- c(tname, "Std err")
    }
    if (!is.null(x$ci)){
        temp <- cbind(temp, x$ci)
        tname <- c(tname, paste(c("lower", "upper"), format(signif(x$conf,2))))
    }
    dimnames(temp) <- list(names(x$coefficients), tname)

    print(temp, ...)
    invisible(x)
    }
   
