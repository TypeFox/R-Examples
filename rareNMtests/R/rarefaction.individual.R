rarefaction.individual <-
function (x, method = "sample-size", q = 0, powerfun = 1, log.scale = FALSE, inds = NULL)  {

    x <- x[x > 0]
    n <- sum(x)

    if(is.null(inds)) {
        length.out <- round(n^powerfun)
        if (log.scale == TRUE) {
            inds <- round(10^(seq(log10(1), log10(n), length.out = length.out)))
        } else {
            inds <- round(seq(1, n, length.out = length.out))
        }
        inds <- inds[!duplicated(inds)]
    }

    if (q == 0) {
        div <- sapply(inds, function(t) rarefy(x, sample = t))
    } else 
    if (q == 1) {
        H1 <- function(x, sample) {
            k <- 1:sample
            fk_m <- sapply(k, function(t) sum(exp((lchoose(x, t) + lchoose(n - x, sample - t)) - lchoose(n, sample))))
            out <- exp(sum((-k/sample)*log(k/sample)*fk_m))
            out
        }
        div <- sapply(inds, function(t) H1(x, sample=t))
    } else
    if (q == 2) {
        H2 <- function(x, sample) {
            k <- 1:sample
            fk_m <- sapply(k, function(t) sum(exp((lchoose(x, t) + lchoose(n - x, sample - t)) - lchoose(n, sample))))
            out <- 1/(sum((k/sample)^2*fk_m))
            out
        }
        div <- sapply(inds, function(t) H2(x, t))
    }

    if (method == "coverage") {
        coverFUN <- function(x, sample) {
            x <- x[x > 0]
            n <- sum(x)
            if (sample < n) {
                num <- lchoose(n - x, sample)
                den <- lchoose(n - 1, sample)
                out <- 1 - sum(x*exp(num-den)/n)
                } else {
                f1 <- sum(x==1)
                f2 <- sum(x==2)
                ifelse(f1 == 0 && f2 == 0, 1, 1 - ((f1/n) * (((n-1)*f1/((n-1)*f1 + 2*f2)))))
            }
        }
    inds <- sapply(inds, function(i) coverFUN(x, sample=i))
    }

    out <- data.frame(inds, div)
    colnames(out) <- c(method, paste("Hill (q=", q, ")", sep=""))
    out
}
