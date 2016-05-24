compresid2way <-
    function(aov, data=NULL, fac=1:2,
	     label = TRUE, numlabel = FALSE, xlab=NULL, ylab=NULL, main=NULL,
	     col=c(2,3,4,4),lty=c(1,1,2,4), pch=c(1,2))
{
    ## Zweck: forget-it-plot   Autor: Stahel  Datum: Dez 89
    ## Arguments:
    ##	 aov	    either a aov object with a formula of the form
    ##		    y ~ a + b , where  a  and  b  are factors
    ##		    or such a formula
    ##	 data	    data frame containing  a  and  b
    ##	 fac	    the two factors used for plotting
    ##	 label	    show levels of factors in the plot
    ##	 numlabel   show effects of factors in the plot
    ##	 col,lty,pch  colors, line types, plotting characters to be used
    ##	   [1]	    positive residuals
    ##	   [2]	    negative residuals
    ##	   [3]	    grid
    ##	   [4]	    labels

    if (inherits(aov,"aov")) {
	lform <- formula(aov)
	if (is.null(data)) {
	    datanm <- as.character(aov$call)[3]
	    if (is.na(datanm))
		stop("no data found")
	    data <- eval(parse(text=datanm))
	}
    } else {
	if (!is.data.frame(data))
	    stop("unsuitable argument  data")
	lform <- aov
	aov <- aov(lform,data)
    }
    lmm <- model.frame(aov)
    fac <- if (is.numeric(fac)) fac+1 else match(fac,names(lmm))
    if (any(is.na(fac)))
	stop("factor(s) not found")
    if (!all(vapply(lmm[,fac], is.factor, NA)))
	stop("variables are not both factors")
    ## coefficients, components of the fit
    lcf <- dummy.coef(aov)
    lic <- lcf[["(Intercept)"]]
    if (is.na(lic)) lic <- 0
    lia <- fac[1]
    lib <- fac[2]
    lfa <- lmm[,lia]
    lfb <- lmm[,lib]
    lcfa <- lcf[[lia]]
    lcfb <- lcf[[lib]]
    lmna <- min(lcfa)
    lmnb <- min(lcfb)
    lcfa <- lcfa-lmna
    lcfb <- lcfb-lmnb
    lic <- lic+lmna+lmnb
    lefa <- lcfa[lfa]
    lefb <- lcfb[lfb]
    lfit <- lic+lefa+lefb
    lfnames <- names(lmm)[c(lia,lib)]
    lyname <- names(lmm)[1]
    ly <- lfit+resid(aov)
    ## prepare plot
    lx <- lefb-lefa
    if (is.null(main))
	main <- format(lform)
    if (is.null(ylab))
	ylab <- lyname
    if (is.null(xlab))
	xlab <- paste("-",paste(lfnames,collapse = "    + "))
    lty <- rep(lty,length = 4)
    if (length(pch) <= 1) pch <- rep(c(pch,pch,1),length = 2)
    lrgy <- range(c(lfit, ly))
    lrgx <- range(lx)
    lht <- 0.05 * diff(lrgy)
    lwd <- 0.05 * diff(lrgx)
    plot(lrgx+lwd*c(-1,1), lrgy+lwd*c(-1,1), type = "n", xlab = "", ylab = ylab)
    mtext(main, 3, 1,
	  cex = par("cex.main"), col = par("col.main"), font = par("font.main"))
    mtext(xlab,1, par("mgp")[1], at = 0)
    ## residuals
    li <- ly > lfit
    if (any(li)) {
	lpch <- if (length(pch) >= length(li)) pch[li] else pch[1]
	segments(lx[li], lfit[li], lx[li], ly[li], lty = lty[1], col = col[1])
	points(lx[li], ly[li], col = col[1], pch = lpch)
    }
    li <- !li
    if (any(li)) {
	lpch <- if (length(pch) >= length(li)) pch[li] else pch[2]
	segments(lx[li], lfit[li], lx[li], ly[li], lty = lty[2], col = col[2])
	points(lx[li], ly[li], col = col[2], pch = lpch)
    }
    ## grid
    lmxa <- max(lcfa)
    segments(lcfb, lic + lcfb, lcfb - lmxa, lic + lmxa + lcfb,
	     lty = lty[3], col = col[3])
    lmxb <- max(lcfb)
    segments( - lcfa, lic + lcfa, lmxb - lcfa, lic + lmxb + lcfa,
	     lty = lty[3], col = col[3])
    ## labels
    if(label)
	text(c(lcfb - lmxa - lwd, lmxb - lcfa + lwd),
	     c(lmxa + lcfb, lmxb + lcfa) + lic + lht,
	     c(levels(lfb), levels(lfa)), col = col[4])
    if(numlabel) {
	ldg <-	- min(0, floor(log10(max(abs(lrgy)))) - 3)
	text(c(lcfb + lwd, - lcfa - lwd), lic + c(lcfb, lcfa) - lht,
	     round(c(lcfb, lcfa), ldg), col = col[4])
    }
    lcf <- list(lic,lcfa,lcfb)
    names(lcf) <- c("(Intercept)",lfnames)
    lcompy <- data.frame(ly,lefa,lefb)
    names(lcompy) <- c(paste("part",lyname,sep = "."),
		       paste("eff",lfnames,sep = "."))
    invisible(list(compy = lcompy,coef = lcf))
}
