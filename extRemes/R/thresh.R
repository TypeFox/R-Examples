mrlplot <- function(x, nint=100, alpha=0.05, na.action=na.fail, xlab="Threshold", ...) {

    x <- na.action(x)
    n <- length(x)

    r <- range(x, finite=TRUE)
    u.i <- matrix(seq(r[1], r[2] - 1,, nint), ncol=1)

    z.alpha <- abs(qnorm(1 - alpha/2))

    mrlfun <- function(u, x, z) {
	eid <- x > u
	xu <- x[eid]
	n <- sum(xu, na.rm=TRUE)
	m <- mean(xu - u)
	s <- sqrt(var(xu)/n)
	res <- c(m - s * z, m, m + s * z)
	return(res)
    } # end of internal 'mrlfun' function.

    out <- t(apply(u.i, 1, mrlfun, x=x, z=z.alpha))
    conf <- paste((1 - alpha)*100, "%", sep="") 
    colnames(out) <- c(paste(conf, " lower", sep=""), "Mean Excess", paste(conf, " upper", sep=""))

    yl <- range(c(out), finite=TRUE)
    plot(u.i, out[,2], type="l", xlab=xlab, ylab="Mean Excess", ylim=yl, ...)
    lines(u.i, out[,1], lty=2, col="gray", lwd=1.5)
    lines(u.i, out[,3], lty=2, col="gray", lwd=1.5)

    invisible(out)
} # end of 'mrlplot' function.

threshrange.plot <- function(x, r, type=c("GP","PP","Exponential"), nint=10, alpha=0.05, na.action=na.fail, set.panels=TRUE, verbose=FALSE, ...) {

    type <- match.arg(type)

    x <- na.action(x)
    n <- length(x)

    if(missing(r)) r <- quantile(x, probs=c(0.75,0.99))
    u.i <- matrix(seq(r[1],r[2],,nint), ncol=1)

    thfun <- function(u, x, type, a, verbose, ...) {
	fit <- try(fevd(x=x, threshold=u, type=type, verbose=verbose, ...), silent=verbose)
	if(verbose) print(fit)
	if(class(fit) != "try-error") {
	    if(!is.element(type,c("PP","Exponential"))) res <- try(ci(fit,type="parameter", alpha=a, R=100, tscale=TRUE), silent=verbose)
	    else res <- try(ci(fit,type="parameter", alpha=a, R=100), silent=verbose)
	    if(verbose) print(res)
# 	    if((class(res) != "try-error") && !is.element(type, c("PP","Exponential")) && fit$method != "Lmoments") {
# 		if(is.element(fit$method, c("MLE","GMLE"))) {
# 		    pcov <- try(parcov.fevd(fit), silent=TRUE)
# 		    if(class(pcov) != "try-error") {
# 		        res["scale",2] <- res["scale",2] - res["shape",2] * u
# 		        d <- rbind(1, -u)
# 		        v <- sqrt(t(d) %*% pcov %*% d)
# 		        res["scale",c(1,3)] <- res["scale",2] + qnorm(c(a/2,1 - a/2))
# 			}
# 		    } else if(fit$method == "Bayesian") {
# 			ptmp <- fit$results
# 			pnames <- colnames(ptmp)
# 			if(is.element("log.scale",pnames)) {
# 			    ptmp[,"log.scale"] <- exp(ptmp[,"log.scale"])
# 			    pnames[pnames=="log.scale"] <- "scale"
# 			    names(ptmp) <- pnames
# 			}
# 			ptmp <- ptmp[,"scale"] <- ptmp[,"scale"] - ptmp[,"shape"] * u
# 			qp <- quantile(c(ptmp[,"scale"]), probs=c(a/2,1-a/2))
# 			res["scale",] <- c(qp[1], ptmp[,"scale"], qp[2])
# 		    } else {
# 		    res["scale",c(1,3)] <- numeric(2)+NA
# 		}
# 	    }
	} else res <- fit
	if(class(res) == "try-error") {
	    if(type=="PP") res <- matrix(NA, 3, 3)
	    else if(type != "Exponential") res <- matrix(NA, 2, 3)
	    else res <- rep(NA, 3)
	}
	return(res)
    } # end of internal 'thfun' function.

    out <- apply(u.i, 1, thfun, x=x, type=type, a=alpha, verbose=verbose, ...)
    if(type == "PP") rownames(out) <- c("low.loc","low.scale","low.shape", "location", "scale", "shape", "up.loc", "up.scale", "up.shape")
    else if(type != "Exponential") rownames(out) <- c("low.t.scale", "low.shape", "t.scale", "shape", "up.t.scale", "up.shape")
    else rownames(out) <- c("low.scale", "scale", "up.scale")

    if(set.panels) {
	if(type=="PP") par(mfrow=c(3,1))
	else if(type != "Exponential") par(mfrow=c(2,1))
	xlb <- ""
    } else xlb <- "Threshold"

    m1 <- deparse(match.call())
    if(type=="PP") {
	yl <- range(c(out[c("low.loc","location","up.loc"),]), finite=TRUE)
	plot(u.i, out["location",], ylim=yl, xlab=xlb, ylab="location", type="b", main=m1)
	for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[c("low.loc","up.loc"),j])

	yl <- range(c(out[c("low.scale","scale","up.scale"),]), finite=TRUE)
	plot(u.i, out["scale",], ylim=yl, xlab=xlb, ylab="scale", type="b")
        for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[c("low.scale","up.scale"),j])

	yl <- range(c(out[c("low.shape","shape","up.shape"),]), finite=TRUE)
        plot(u.i, out["shape",], ylim=yl, xlab="Threshold", ylab="shape", type="b")
        for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[c("low.shape","up.shape"),j])
    } else if(type != "Exponential") {
	yl <- range(c(out[c("low.t.scale","t.scale","up.t.scale"),]), finite=TRUE)
        plot(u.i, out["t.scale",], ylim=yl, xlab=xlb, ylab="reparameterized scale", type="b", main=m1)
        for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[c("low.t.scale","up.t.scale"),j])

        yl <- range(c(out[c("low.shape","shape","up.shape"),]), finite=TRUE)
        plot(u.i, out["shape",], ylim=yl, xlab="Threshold", ylab="shape", type="b")
        for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[c("low.shape","up.shape"),j])
    } else {
	yl <- range(c(out[c("low.scale","scale","up.scale"),]), finite=TRUE)
        plot(u.i, out["scale",], ylim=yl, xlab="Threshold", ylab="scale", type="b", main=m1)
        for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[c("low.scale","up.scale"),j])
    }

    invisible(t(out))
} # end of 'threshrange.plot' function.
