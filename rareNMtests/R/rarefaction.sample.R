rarefaction.sample <-
function (x, method = "sample-size", q = 0)  {

    x <- x[,apply(x, 2, sum)>0]
    T <- dim(x)[1]
    sp <- dim(x)[2]
    U <- sum(x>0)
    inds <- 1:T

    if (q == 0) {
        div <- specaccum(x)$richness
    } else 
    if (q == 1) {
        H1 <- function(x, t) {
            Y <- colSums(x>0)
            Y <- Y[Y > 0]
            k <- 1:t
            Qk_t <- sapply(k, function(k) sum(exp((lchoose(Y, k) + lchoose(T - Y, t - k)) - lchoose(T, t))))
            out <- exp(sum((-k/(t*U/T))*log(k/(t*U/T))*Qk_t))
            out
        }
        div <- sapply(inds, function(t) H1(x, t = t))
    } else
    if (q == 2) {
        H2 <- function(x, t) {
            Y <- colSums(x>0)
            Y <- Y[Y > 0]
            k <- 1:t
            Qk_t <- sapply(k, function(k) sum(exp((lchoose(Y, k) + lchoose(T - Y, t - k)) - lchoose(T, t))))
            out <- 1/(sum((k/(t*U/T))^2*Qk_t))
            out
        }
        div <- sapply(inds, function(t) H2(x, t = t))
    }

    if (method == "coverage") {
        coverFUN <- function(x, sample) {
            Y <- colSums(x>0)
            Y <- Y[Y > 0]
            if (sample < T) {
                num <- lchoose(T - Y, sample)
                den <- lchoose(T - 1, sample)
                out <- 1 - sum(Y*exp(num-den)/U)
                } else {
                f1 <- sum(Y==1)
                f2 <- sum(Y==2)
                ifelse(f1 == 0 && f2 == 0, 1, 1 - ((f1/U) * (((U-1)*f1/((U-1)*f1 + 2*f2)))))
            }
        }
    inds <- sapply(inds, function(i) coverFUN(x, sample=i))
    }

    out <- data.frame(inds, div)
    colnames(out) <- c(method, paste("Hill (q=", q, ")", sep=""))
    out
}
