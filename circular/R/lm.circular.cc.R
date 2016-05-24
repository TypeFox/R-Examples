
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   lm.circular.cc function                                 #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: August, 10, 2006                                  #
#   Version: 0.2-3                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

lm.circular.cc <- function(y, x, order = 1, level = 0.05, control.circular=list()) {
    # Handling missing values
    ok <- complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]
    if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }

    if (is.circular(y)) {
       datacircularp <- circularp(y)     
    } else if (is.circular(x)) {
              datacircularp <- circularp(x)     
    } else {
       datacircularp <- list(type="angles", units="radians", template="none", modulo="2pi", zero=0, rotation="counter")
    }

    dc <- control.circular
    if (is.null(dc$type))
       dc$type <- datacircularp$type
    if (is.null(dc$units))
       dc$units <- datacircularp$units
    if (is.null(dc$template))
       dc$template <- datacircularp$template
    if (is.null(dc$modulo))
       dc$modulo <- datacircularp$modulo
    if (is.null(dc$zero))
       dc$zero <- datacircularp$zero
    if (is.null(dc$rotation))
       dc$rotation <- datacircularp$rotation
    
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "circularp") <- attr(x, "class") <- NULL
    y <- conversion.circular(y, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(y, "circularp") <- attr(y, "class") <- NULL
    
    circ.lm <- LmCircularccRad(y, x, order)
    
    circ.lm$call <- match.call()
    circ.lm$fitted <- conversion.circular(circular(circ.lm$fitted), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    circ.lm$residuals <- conversion.circular(circular(circ.lm$residuals), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)    
    if (circ.lm$p.values[1] > level & circ.lm$p.values[2] > level)
        circ.lm$message <- paste("Higher order terms are not significant at the ", level, " level", sep = "")
    else circ.lm$message <- paste("Higher order terms are significant at the ", level, " level", sep = "")
    class(circ.lm) <- "lm.circular.cc"
    return(circ.lm)
}

LmCircularccRad <- function(y, x, order) {
    n <- length(x)
    cy <- cos(y)
    sy <- sin(y)
    order.matrix <- t(matrix(rep(c(1:order), n), ncol = n))
    cos.x <- cos(x * order.matrix)
    sin.x <- sin(x * order.matrix)
    cos.lm <- lm(cy ~ cos.x + sin.x)
    sin.lm <- lm(sy ~ cos.x + sin.x)
    cos.fit <- cos.lm$fitted
    sin.fit <- sin.lm$fitted
    g1.sq <- t(cos.fit) %*% cos.fit
    g2.sq <- t(sin.fit) %*% sin.fit
    rho <- sqrt((g1.sq + g2.sq)/n)
    y.fitted <- atan2(sin.fit, cos.fit)
    Y1 <- cy
    Y2 <- sy
    ones <- matrix(1, n, 1)
    X <- cbind(ones, cos.x, sin.x)
    W <- cbind(cos((order + 1) * x), sin((order + 1) * x))
    M <- X %*% solve(t(X) %*% X) %*% t(X)
    I <- diag(n)
    H <- t(W) %*% (I - M) %*% W
    N <- W %*% solve(H) %*% t(W)
    cc <- n - (2 * order + 1)
    N1 <- t(Y1) %*% (I - M) %*% N %*% (I - M) %*% Y1
    D1 <- t(Y1) %*% (I - M) %*% Y1
    T1 <- cc * (N1/D1)
    N2 <- t(Y2) %*% (I - M) %*% N %*% (I - M) %*% Y2
    D2 <- t(Y2) %*% (I - M) %*% Y2
    T2 <- cc * (N2/D2)
    p1 <- 1 - pchisq(T1, 2)
    p2 <- 1 - pchisq(T2, 2)
    pvalues <- cbind(p1, p2)
    circ.lm <- list()
    circ.lm$rho <- rho
    circ.lm$fitted <- y.fitted %% (2 * pi)
    circ.lm$x <- cbind(x, y)
    circ.lm$residuals <- (y - y.fitted) %% (2 * pi)
    circ.lm$coefficients <- cbind(cos.lm$coefficients, sin.lm$coefficients)
    circ.lm$p.values <- pvalues
    circ.lm$A.k <- mean(cos(circ.lm$residuals))
    circ.lm$kappa <- A1inv(circ.lm$A.k)
    return(circ.lm)
}
