senmv <-
function (y, gamma = 1, method = NULL, inner = 0, trim = 2.5, lambda = 1/2,
              tau = 0, TonT = FALSE)
    {
        if (is.vector(y)) {
            y <- y[!is.na(y)]
            treat <- y/2
            cont <- (-y/2)
            y <- cbind(treat, cont)
        }
        stopifnot((0 <= inner) & (inner <= trim))
        stopifnot((lambda > 0) & (lambda < 1))
        stopifnot(gamma >= 1)
        stopifnot(min(apply(!is.na(y), 1, sum)) >= 2)
        m <- 1
        m1 <- 1
        m2 <- 1
        if (!is.null(method))
            stopifnot(is.element(method, c("h", "i", "w", "t","q","s","l")))
        wgt <- (m > 1)
        if (!is.null(method)) {
            if (is.element(method, c("w","q","s","l")))
                wgt <- TRUE
        }
        vc <- (sum(is.na(as.vector(y)))) > 0
        if (wgt & vc)
            warning("Weighting of sets is not available for matched sets of variable sizes.")
        stopifnot(FALSE == (wgt & vc))
        if (!is.null(method)) {
            if (method == "h") {
                inner <- 0
                trim <- 2.5
                lambda <- 1/2
                m1 <- 1
                m2 <- 1
                m <- 1
                TonT <- FALSE
            }
            if (method == "i") {
                inner <- 1/2
                trim <- 2.5
                lambda <- 1/2
                m1 <- 1
                m2 <- 1
                m <- 1
                TonT <- FALSE
            }
            if (method == "w") {
                inner <- 0
                trim <- 3
                lambda <- 1/2
                m1 <- 12
                m2 <- 20
                m <- 20
                TonT <- FALSE
            }
            if (method == "t") {
                inner <- 0
                trim <- Inf
                lambda <- 1/2
                m1 <- 1
                m2 <- 1
                m <- 1
                TonT <- TRUE
            }
            if (method == "q") {
                inner <- 0
                trim <- 3
                lambda <- 1/2
                m1 <- 2
                m2 <- 2
                m <- 2
                TonT <- FALSE
            }
            if (method == "s") {
                inner <- 0
                trim <- 3
                lambda <- 1/2
                m1 <- 16
                m2 <- 20
                m <- 20
                TonT <- FALSE
            }
            if (method == "l") {
                inner <- 0
                trim <- 2
                lambda <- 1/2
                m1 <- 12
                m2 <- 19
                m <- 20
                TonT <- FALSE
            }
        }
        if (!(tau == 0))
            y[, 1] <- y[, 1] - tau
        ms <- mscorev(y, inner = inner, trim = trim, qu = lambda,
                      TonT = TonT)
        if (m > 1)
            separable1v(newurks(ms, m1 = m1, m2 = m2, m = m), gamma = gamma)
        else if (m == 1)
            separable1v(ms, gamma = gamma)
    }
