`GetFitFD` <-
function(z, MeanZeroQ=FALSE, algorithm=c("emle", "wmle"), ciQ=FALSE){
    dlim <- 0.4999
    y<-z
    if (!MeanZeroQ)
        y<-z-mean(z)
    algorithm <- match.arg(algorithm)
    if (identical(algorithm, "emle"))
        Entropyd<-function(d) -LLFD(d,y) else {
        Ip <- (spec.pgram(z, fast=FALSE, detrend=FALSE, plot=FALSE, taper=0)$spec)/(2*pi)
        n <- length(y)
        Entropyd <- function(d)  (n/2)*log(mean((Ip/sdfFD(d, n))))
        }
    out<-optimize(f=Entropyd, interval=c(-dlim, dlim))
    d<-out[[1]] #sometimes optimize incorrectly chooses the endpoint (very rare)
    if (d == -dlim || d==dlim) {
        warning("end point reached. Switching to L-BFGS-B ...")
        H0 <- HurstK(y)
        d0 <- 0.5-H0
        out<-optim(d0, fn=Entropyd, method="L-BFGS-B", lower=-0.49, upper=0.49)
        err<-out$convergence
        if (err != 0) {
            warning(" err = ", err , ". Trying Nelder-Mead...")
            out<-optim(d0, fn=Entropyd, method="Nelder-Mead")
        }
    }
    dHat <- out[[1]]
    ans<-list(d=dHat, loglikelihood=-out[[2]], alpha=1-2*dHat, algorithm=algorithm)
    if (ciQ) {
        f <- function(x) exp(Entropyd(dHat)-Entropyd(x))-0.1464
        out1<-try(ciU <- uniroot(f, lower=dHat, upper=dlim)$root, silent=TRUE)
        out2<-try(ciL <- uniroot(f, upper=dHat, lower=-dlim)$root, silent=TRUE)
        if(is.numeric(out1)&&is.numeric(out2)){
            ci <- c(ciL, ciU)
            f <- function(x) exp(Entropyd(dHat)-Entropyd(x))-0.07
            out1 <- try(U <- uniroot(f, lower=dHat, upper=dlim)$root, silent=TRUE)
            out2 <- try(L <- uniroot(f, upper=dHat, lower=-dlim)$root, silent=TRUE)
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
                text(x=ciU, y=0.8, labels=bquote(hat(d) == .(round(dHat,3))))
                polygon(x=c(ciL, ciL, ciU, ciU), y=c(ylo, 0.1464, 0.1464, ylo), col="lightblue")
                text(x=ciL+0.3*(ciU-ciL), y=ylo+(0.1464-ylo)/2, "95% c.i.")
                lines(x=c(dHat, dHat), y=c(ylo,1), lty=2, col="orange")
                }
            } else ci <- NA
        ans <- c(ans, list(ci=ci))
        }
    ans
}
