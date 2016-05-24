`GetFitFGN` <-
function(z, MeanZeroQ=FALSE, algorithm=c("emle", "wmle"), ciQ=FALSE){
    y<-z
    if (!MeanZeroQ)
        y<-z-mean(z)
    algorithm <- match.arg(algorithm)
    if (identical(algorithm, "emle"))
        EntropyH<-function(H) -LLFGN(H,y) else {
        Ip <- (spec.pgram(z, fast=FALSE, detrend=FALSE, plot=FALSE, taper=0)$spec)/(2*pi)
        EntropyH<-function(H) {
            fp <- sdfFGN(H, length(z))
            sigHat <- mean(Ip/fp)
            2*sum(log(sigHat*fp))
            }
        }
    out<-optimize(f=EntropyH, interval=c(0.01, 0.99))
    H<-out[[1]] #sometimes optimize incorrectly chooses the endpoint (very rare)
    if (H == 0.01 || H==0.99) {
        warning("end point reached. Switching to L-BFGS-B ...")
        H0<-HurstK(y)
        out<-optim(H0, fn=EntropyH, method="L-BFGS-B", lower=0.01, upper=0.99)
        err<-out$convergence
        if (err != 0) {
            warning(" err = ", err , ". Trying Nelder-Mead...")
            out<-optim(H0, fn=EntropyH, method="Nelder-Mead")
        }
    }
    HHat <- out[[1]]
    ans<-list(H=HHat, loglikelihood=-out[[2]], alpha=2-2*HHat, algorithm=algorithm)
    if (ciQ) {
        f <- function(x) exp(EntropyH(HHat)-EntropyH(x))-0.1464
        out1<-try(ciU <- uniroot(f, lower=HHat, upper=0.97)$root, silent=TRUE)
        out2<-try(ciL <- uniroot(f, upper=HHat, lower=0.03)$root, silent=TRUE)
        if(is.numeric(out1)&&is.numeric(out2)){
            ci <- c(ciL, ciU)
            f <- function(x) exp(EntropyH(HHat)-EntropyH(x))-0.0044
            out1 <- try(U <- uniroot(f, lower=HHat, upper=0.97)$root, silent=TRUE)
            out2 <- try(L <- uniroot(f, upper=HHat, lower=0.03)$root, silent=TRUE)
            if (is.numeric(out1)&&is.numeric(out2)) {
                x <- seq(from=L, to=U, length.out=50)
                fx <- numeric(length(x))
                for (i in 1:length(x))
                    fx[i] <- f(x[i])
                if (identical(algorithm,"emle"))
                    ylab <- "Relative likelihood (exact)" else
                    ylab <- "Relative likelihood (Whittle)"
                plot(x, fx, type="l", xlab="H", ylab=ylab)
                ylo <- par("usr")[3]
                text(x=ciU, y=0.8, labels=bquote(hat(H) == .(round(HHat,3))))
                polygon(x=c(ciL, ciL, ciU, ciU), y=c(ylo, 0.1464, 0.1464, ylo), col="lightblue")
                text(x=ciL+0.3*(ciU-ciL), y=ylo+(0.1464-ylo)/2, "95% c.i.")
                lines(x=c(HHat, HHat), y=c(ylo,1), lty=2, col="orange")
                }
            } else ci <- NA
            ans <- c(ans, list(ci=ci))
        }
    ans
}
