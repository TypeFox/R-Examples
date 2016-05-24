bwfw1.hmm <-
function (x, pii, A, f0, f1) 
{
    NUM <- length(x)
    f0x <- dnorm(x, f0[1], f0[2])
    f1x <- dnorm(x, f1[1], f1[2])
    alpha <- matrix(rep(0, NUM * 2), NUM, 2, byrow = TRUE)
    c0 <- rep(0, NUM)
    alpha[1, 1] <- pii[1] * f0x[1]
    alpha[1, 2] <- pii[2] * f1x[1]
    c0[1] <- 1/sum(alpha[1, ])
    alpha[1, ] <- c0[1] * alpha[1, ]
    for (k in 1:(NUM - 1)) {
        alpha[k + 1, 1] <- (alpha[k, 1] * A[1, 1] + alpha[k, 
            2] * A[2, 1]) * f0x[k + 1]
        alpha[k + 1, 2] <- (alpha[k, 1] * A[1, 2] + alpha[k, 
            2] * A[2, 2]) * f1x[k + 1]
        c0[k + 1] <- 1/sum(alpha[k + 1, ])
        alpha[k + 1, ] <- c0[k + 1] * alpha[k + 1, ]
    }
    beta <- matrix(rep(0, NUM * 2), NUM, 2, byrow = TRUE)
    beta[NUM, 1] <- c0[NUM]
    beta[NUM, 2] <- c0[NUM]
    for (k in (NUM - 1):1) {
        beta[k, 1] <- A[1, 1] * f0x[k + 1] * beta[k + 1, 1] + 
            A[1, 2] * f1x[k + 1] * beta[k + 1, 2]
        beta[k, 2] <- A[2, 1] * f0x[k + 1] * beta[k + 1, 1] + 
            A[2, 2] * f1x[k + 1] * beta[k + 1, 2]
        beta[k, ] <- c0[k] * beta[k, ]
    }
    lfdr <- rep(0, NUM)
    for (k in 1:NUM) {
        q1 <- alpha[k, 1] * beta[k, 1]
        q2 <- alpha[k, 2] * beta[k, 2]
        lfdr[k] <- q1/(q1 + q2)
    }
    gamma <- matrix(1:(NUM * 2), NUM, 2, byrow = TRUE)
    gamma[NUM, ] <- c(lfdr[NUM], 1 - lfdr[NUM])
    dgamma <- array(rep(0, (NUM - 1) * 4), c(2, 2, (NUM - 1)))
    for (k in 1:(NUM - 1)) {
        denom <- 0
        for (i in 0:1) {
            for (j in 0:1) {
                fx <- (1 - j) * f0x[k + 1] + j * f1x[k + 1]
                denom <- denom + alpha[k, i + 1] * A[i + 1, j + 
                  1] * fx * beta[k + 1, j + 1]
            }
        }
        for (i in 0:1) {
            gamma[k, i + 1] <- 0
            for (j in 0:1) {
                fx <- (1 - j) * f0x[k + 1] + j * f1x[k + 1]
                dgamma[i + 1, j + 1, k] <- alpha[k, i + 1] * 
                  A[i + 1, j + 1] * fx * beta[k + 1, j + 1]/denom
                gamma[k, i + 1] <- gamma[k, i + 1] + dgamma[i + 
                  1, j + 1, k]
            }
        }
    }
    bwfw.var <- list(bw = alpha, fw = beta, lf = lfdr, pr = gamma, 
        ts = dgamma, rescale = c0)
    return(bwfw.var)
}
