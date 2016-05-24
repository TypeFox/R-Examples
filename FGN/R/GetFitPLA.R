`GetFitPLA` <-
function(z, MeanZeroQ=FALSE, algorithm=c("emle", "wmle"), ciQ=FALSE){
    y<-z
    if (!MeanZeroQ)
        y<-z-mean(z)
    algorithm <- match.arg(algorithm)
    if (identical(algorithm, "emle"))
        Entropya<-function(a) -LLPLA(a,y) else {
        Ip <- (spec.pgram(z, fast=FALSE, detrend=FALSE, plot=FALSE, taper=0)$spec)/(2*pi)
        Entropya<-function(a) {
            fp <- sdfPLA(a, length(z))
            sigHat <- mean(Ip/fp)
            2*sum(log(sigHat*fp))
            }
        }
    out<-optimize(f=Entropya, interval=c(0.001, 1.999))
    a<-out[[1]] #sometimes optimize incorrectly chooses the endpoint (very rare)
    if (a == 0.001 || a==1.999) {
        warning("end point reached. Switching to L-BFGS-B ...")
        H0 <- HurstK(y)
        a0 <- 2-2*H0
        out<-optim(a0, fn=Entropya, method="L-BFGS-B", lower=0.01, upper=2.0)
        err<-out$convergence
        if (err != 0) {
            warning(" err = ", err , ". Trying Nelder-Mead...")
            out<-optim(a0, fn=Entropya, method="Nelder-Mead")
        }
    }
    aHat <- out[[1]]
    ans<-list(a=aHat, loglikelihood=-out[[2]], alpha=aHat, algorithm=algorithm)
    if (ciQ) {
        f <- function(x) exp(Entropya(aHat)-Entropya(x))-0.1464
        out1<-try(ciU <- uniroot(f, lower=aHat, upper=1.99)$root, silent=TRUE)
        out2<-try(ciL <- uniroot(f, upper=aHat, lower=0.01)$root, silent=TRUE)
        if(is.numeric(out1)&&is.numeric(out2)){
            ci <- c(ciL, ciU)
            f <- function(x) exp(Entropya(aHat)-Entropya(x))-0.0044
            out1 <- try(U <- uniroot(f, lower=aHat, upper=1.99)$root, silent=TRUE)
            out2 <- try(L <- uniroot(f, upper=aHat, lower=0.01)$root, silent=TRUE)
            if (is.numeric(out1)&&is.numeric(out2)) {
                x <- seq(from=L, to=U, length.out=50)
                fx <- numeric(length(x))
                for (i in 1:length(x))
                    fx[i] <- f(x[i])
                if (identical(algorithm,"emle"))
                    ylab <- "Relative likelihood (exact)" else
                    ylab <- "Relative likelihood (Whittle)"
                plot(x, fx, type="l", xlab="a", ylab=ylab)
                ylo <- par("usr")[3]
                text(x=ciU, y=0.8, labels=bquote(hat(a) == .(round(aHat,3))))
                polygon(x=c(ciL, ciL, ciU, ciU), y=c(ylo, 0.1464, 0.1464, ylo), col="lightblue")
                text(x=ciL+0.3*(ciU-ciL), y=ylo+(0.1464-ylo)/2, "95% c.i.")
                lines(x=c(aHat, aHat), y=c(ylo,1), lty=2, col="orange")
                }
            } else ci <- NA
         ans <- c(ans, list(ci=ci))
        }
    ans
}
