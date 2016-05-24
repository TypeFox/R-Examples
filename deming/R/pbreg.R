pbreg <- function(formula, data, subset, weights, na.action, conf=.95,
                  nboot=0, method=1,
                  eps=sqrt(.Machine$double.eps),
                  x=FALSE, y=FALSE, model=TRUE) {
    Call <- match.call()

    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    # then evaluate it
    indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp[[1]] <- as.name('model.frame')  # change the function called
    mf <- eval(temp, parent.frame())
    Terms <- terms(mf)
    n <- nrow(mf)
    if (n <3) stop("less than 3 non-missing observations in the data")

    X <- model.matrix(Terms, mf)
    if (!attr(Terms, "intercept")) stop ("an intercept is required")
    if (ncol(X) != 2) stop("P-B regression requires a single predictor variable")
    xx <- X[,2]

    Y <- model.response(mf, type="numeric")
    if (is.null(Y))
        stop ("a response variable is required")
    wt <- model.weights(mf)
    if (is.null(wt)) wt <- rep(1.0, n)

    if (floor(method) != method || method <1 || method >3) 
        stop("method must be an integer between 1 and 3")
    if (!is.numeric(eps) || eps <=0) stop("invalid value for eps")

    if (conf <0 || conf >=1) stop("invalid confidence interval limit")

    if (conf !=0 && nboot==0 && any(wt != wt[1]))
        stop("with weighted data, a confidence interval requires bootstrap")
    
    fit <- pbreg.fit(xx, Y, wt, method, conf*(nboot==0), eps)
    names(fit$coefficients) <-  dimnames(X)[[2]]
    yhat <- fit$coefficients[1] + xx*fit$coefficients[2]
    fit$residuals <- Y-yhat
    
    if (nboot > 0) {
        bfun <- function(n) {
            temp <- rbinom(n, 1, (1+sqrt(5))/sqrt(20))
            ifelse(temp==1, 1-sqrt(5), 1+sqrt(5))/2
        }

        # Set up for the wild bootstrap, which adds random noise to
        #  the original x,y values
        fcoef <- fit$coefficients
#        u <- (xx + fcoef[2]*(Y - fcoef[1]))/(1+ fcoef[2]^2) #projection onto fit
#        resid <- list(x= u-xx, y=(fcoef[1] + fcoef[2]*u) - Y)
        resid <- list(y=fit$residuals, x= xx - (Y-fcoef[1])/fcoef[2])

        wild <- function(data, resid) {
            rb <- bfun(length(data$x))
            newx <- data$x + resid$x*rb
            newy <- data$y + resid$y*rb
            list(x=newx, y=newy, w=data$w)
        }
        dfun <- function(data) {  #will be called by boot function
            pbreg.fit(data$x, data$y, data$w, method, conf=0, eps)$coefficients
        }
        boot.out <- boot(list(x=xx, y=Y, w= wt), dfun, nboot, sim="parametric",
                         ran.gen=wild, mle=resid)
        temp1 <- boot.ci(boot.out, type="perc", index=1, conf=conf)
        temp2 <- boot.ci(boot.out, type="perc", index=2, conf=conf)
        fit$boot <-boot.out
        fit$ci   <- rbind(temp1$perc[,4:5], temp2$perc[,4:5])
        fit$var <- var(boot.out$t)
    }

    if (!is.null(fit$ci)) 
        dimnames(fit$ci) <- list(names(fit$coef), 
                                 paste(c("lower", "upper"), format(conf)))
    fit$n <- n
    if (x) fit$x <- X
    if (y) fit$y <- Y
    if (model) fit$model <- mf
    fit$method <- method
    
    na.action <- attr(mf, "na.action")
    if (length(na.action)) fit$na.action <- na.action
    fit$terms <- Terms
    fit$conf <- conf
    fit$call <- Call
    class(fit) <- "pbreg"
    fit
}

print.pbreg <- function(x, ...) {
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
        tname <- c(tname, paste(c("lower", "upper"), format(x$conf)))
    }
    dimnames(temp) <- list(names(x$coefficients), tname)

    print(temp, ...)
    invisible(x)
    }
   
