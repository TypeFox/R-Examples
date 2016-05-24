senmw <-
function (y, gamma = 1, method = NULL, inner = 0, trim = 3, lambda = 1/2, 
              tau = 0, m1=1, m2=1, m=1) 
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
        stopifnot((m1>=1)&(m2>=m1)&(m>=m2))
		TonT<-FALSE
        if (!is.null(method)) 
            stopifnot(is.element(method, c("h", "t", "w", "f", "q","s","l","p")))
        if (!is.null(method)) {
            if (is.element(method, c("w","f","q","s","l"))) 
                wgt <- TRUE
        }
        vc <- (sum(is.na(as.vector(y)))) > 0
		if (vc) warning("For senmw, y cannot include NAs.  All matched sets must have the same number of controls. 
						If your study has variable numbers of controls, consider using senmv in the sensitivitymv package.")
		stopifnot(!vc)
        if (!is.null(method)) {
            if (method == "h") {
                inner <- 0
                trim <- 3
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
			if (method == "f") {
                inner <- 0
                trim <- 3
                lambda <- 1/2
                m1 <- 14
                m2 <- 20
                m <- 20
                TonT <- FALSE
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
                trim <- 3
                lambda <- 1/2
                m1 <- 12
                m2 <- 19
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
                TonT <- FALSE
            }
			if (method == "p") {
                inner <- 0.5
                trim <- 2
                lambda <- 1/2
                m1 <- 1
                m2 <- 1
                m <- 1
                TonT <- FALSE
            }
        }
        if (!(tau == 0)) 
            y[, 1] <- y[, 1] - tau
        ms <- mscorev(y, inner = inner, trim = trim, qu = lambda, TonT=FALSE)
        if (m > 1) 
            separable1k(newurks(ms, m1 = m1, m2 = m2, m = m), gamma = gamma)
        else if (m == 1) 
            separable1k(ms, gamma = gamma)
    }

