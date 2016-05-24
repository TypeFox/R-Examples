Cline.plot <-
function(cfit){
	llogit.p <-function(x,coef){
		u <- as.numeric(coef)[1]
		v <- as.numeric(coef)[2]
		x^v/(x^v + (1-x)^v*exp(u))
		}
	Barton.p <- function(x,coef){
		coef <- as.numeric(coef)
		x+2*x*(1-x)*(coef[1]+coef[2]*(2*x-1))
		}
	Beta.p <- function(x,coef){
		coef <- as.numeric(coef)
		pbeta(x,coef[1]*coef[2],(1-coef[1])*coef[2])
		}
	L.p <- function(x,coef){
		coef <- as.numeric(coef)
		exp(coef[1]+coef[2]*x)/(1+exp(coef[1]+coef[2]*x))
		}
	Mn.p <- function(x,coef){
		coef <- as.numeric(coef)
		Hom <- exp(coef[3]+x*coef[4])/(1+exp(coef[1]+x*coef[2])+exp(coef[3]+x*coef[4]))
		Het <- exp(coef[1]+x*coef[2])/(1+exp(coef[1]+x*coef[2])+exp(coef[3]+x*coef[4]))
		data.frame(Hom,Het)
		}
	Rich.p <- function(x,coef){
		coef <- as.numeric(coef)
		U <- coef[1];L <- coef[2]; b <- coef[3]; m <- coef[4]
		U+(L-U)/(1+exp(b*(x-m)))
		}
	nmod <- 0
	Models <- character()
	if(!is.null(cfit$Barton)){nmod <- nmod+1; Models[1] <- "Barton"}
	if(!is.null(cfit$beta)){nmod <- nmod+1; Models[2] <- "beta"}
	if(!is.null(cfit$logit.logistic)){nmod <- nmod+1; Models[3] <- "logit-logistic"}
	if(!is.null(cfit$binom)){nmod <- nmod+1; Models[4] <- "binomial"}
	if(!is.null(cfit$multinom)){nmod <- nmod+1; Models[5] <- "multinomial"}
	if(!is.null(cfit$Richards)){nmod <- nmod+1; Models[6] <- "Richards"}
	
	par(mfrow=c(nmod,3),mar=c(2,2,.5,.5),mgp=c(1,.25,0),tcl=-.2)
	
	Sx <- seq(from=0,to=1,length.out=200)
	cols <- c("black","red")
	
	if(!is.null(cfit$Barton)){
		plot(0:1,0:1,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," Barton")))
		for(i in 1:dim(cfit$Barton)[1]){
			lines(Sx,Barton.p(Sx,cfit$Barton[i,1:2]),col=cols[1+1*cfit$Barton[i,7]])
			}
		plot(cfit$Barton[,1:2],xlab=expression(paste(italic(a))),ylab=expression(paste(italic(b))),col=cols[1+1*cfit$Barton[,7]])
		o <- order(cfit$Barton[,6])
		XQ <- qchisq(rank(cfit$Barton[,6])/(1+length(cfit$Barton[,6])),df=2,lower.tail=FALSE)
		plot(XQ,cfit$Barton[,5],xlab=expression(paste(chi^2," quantiles")),ylab=expression(paste(italic(D)^2)),col=cols[1+1*cfit$Barton[,7]])
	lines(XQ[o],XQ[o])
		}
	if(!is.null(cfit$beta)){
			plot(0:1,0:1,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," beta")))
		for(i in 1:dim(cfit$beta)[1]){
			lines(Sx,Beta.p(Sx,cfit$beta[i,1:2]),col=cols[1+1*cfit$beta[i,7]])
			}
		plot(cfit$beta[,1:2],xlab=expression(paste(italic(mu))),ylab=expression(paste(italic(nu))),col=cols[1+1*cfit$beta[,7]])
		o <- order(cfit$beta[,6])
		XQ <- qchisq(rank(cfit$beta[,6])/(1+length(cfit$beta[,6])),df=2,lower.tail=FALSE)
		plot(XQ,cfit$beta[,5],xlab=expression(paste(chi^2," quantiles")),ylab=expression(paste(italic(D)^2)),col=cols[1+1*cfit$beta[,7]])
	lines(XQ[o],XQ[o])
}
	if(!is.null(cfit$logit.logistic)){
		plot(0:1,0:1,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," logit-logistic")))
		for(i in 1:dim(cfit$logit)[1]){
			lines(Sx,llogit.p(Sx,cfit$logit[i,1:2]),col=cols[1+1*cfit$logit[i,7]])
			}
		plot(cfit$logit[,1:2],xlab=expression(paste(italic(u))),ylab=expression(paste(italic(v))),col=cols[1+1*cfit$logit[,7]])
		o <- order(cfit$logit[,6])
		XQ <- qchisq(rank(cfit$logit[,6])/(1+length(cfit$logit[,6])),df=2,lower.tail=FALSE)
		plot(XQ,cfit$logit[,5],xlab=expression(paste(chi^2," quantiles")),ylab=expression(paste(italic(D)^2)),col=cols[1+1*cfit$logit[,7]])
	lines(XQ[o],XQ[o])
	}
	if(!is.null(cfit$binom)){
		plot(0:1,0:1,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," binomial")))
		for(i in 1:dim(cfit$binom)[1]){
			lines(Sx,L.p(Sx,cfit$binom[i,1:2]),col=cols[1+1*cfit$binom[i,7]])
			}
		plot(cfit$binom[,1:2],xlab=expression(paste(italic(alpha))),ylab=expression(paste(italic(beta))),col=cols[1+1*cfit$binom[,7]])
		o <- order(cfit$binom[,6])
		XQ <- qchisq(rank(cfit$binom[,6])/(1+length(cfit$binom[,6])),df=2,lower.tail=FALSE)
		plot(XQ,cfit$binom[,5],xlab=expression(paste(chi^2," quantiles")),ylab=expression(paste(italic(D)^2)),col=cols[1+1*cfit$binom[,7]])
	lines(XQ[o],XQ[o])
	}
	if(!is.null(cfit$Richards)){
		plot(0:1,0:1,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," Richards")))
		for(i in 1:dim(cfit$Richards)[1]){
			lines(Sx,Rich.p(Sx,cfit$Richards[i,1:4]),col=cols[1+1*cfit$Richards[i,9]])
			}
		plot(1,1,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
		o <- order(cfit$Richards[,8])
		XQ <- qchisq(rank(cfit$Richards[,8])/(1+length(cfit$Richards[,6])),df=2,lower.tail=FALSE)
		plot(XQ,cfit$Richards[,7],xlab=expression(paste(chi^2," quantiles")),ylab=expression(paste(italic(D)^2)),col=cols[1+1*cfit$Richards[,9]])
	lines(XQ[o],XQ[o])
}
	if(!is.null(cfit$multinom)){
		plot(0:1,0:1,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," multinomial")))
		for(i in 1:dim(cfit$multinom)[1]){
			lines(Sx,Mn.p(Sx,cfit$multinom[i,1:4])$Hom,col=cols[1+1*cfit$multinom[i,9]])
			}
		for(i in 1:dim(cfit$multinom)[1]){
			lines(Sx,Mn.p(Sx,cfit$multinom[i,1:4])$Het,col="grey")
			if(cfit$multinom[i,9]){lines(Sx,Mn.p(Sx,cfit$multinom[i,1:4])$Het,col="red",lty=2)}
			}
		plot(1,1,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
		o <- order(cfit$multinom[,8])
		XQ <- qchisq(rank(cfit$multinom[,8])/(1+length(cfit$multinom[,6])),df=2,lower.tail=FALSE)
		plot(XQ,cfit$multinom[,7],xlab=expression(paste(chi^2," quantiles")),ylab=expression(paste(italic(D)^2)),col=cols[1+1*cfit$multinom[,9]])
	lines(XQ[o],XQ[o])
	}
}
