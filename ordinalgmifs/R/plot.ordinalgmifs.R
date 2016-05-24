plot.ordinalgmifs <-
function (x, type = "trace", xlab=NULL, ylab=NULL, main=NULL, ...) 
{
	if (is.null(xlab)) xlab="Step"
	if (is.null(ylab)) {
		if (type == "AIC") {
			ylab="AIC"
		} else if (type == "BIC") {
			ylab="BIC"
		} else if (type == "trace") {
			ylab=expression(hat(beta))
		} else if (type == "logLik") {
			ylab = "logLikelihood"
		}
	}
    if (is.null(x$x)) {
        stop("Penalized model not requested\n")
    }
    else if (type == "trace") {
        matplot(x$beta, ylab = ylab, xlab = xlab, 
            type = "l")
    }
    else if (type == "AIC") {
        plot(x$AIC, xlab = xlab, ylab = ylab)
    }
    else if (type == "BIC") {
        plot(x$BIC, xlab = xlab, ylab = ylab)
    }
    else if (type == "logLik") {
        plot(x$logLik, xlab = xlab, ylab = ylab)
    }
    if (is.null(main)) {
    	if (x$probability.model == "Cumulative" | x$probability.model == 
        "ForwardCR" | x$probability.model == "BackwardCR") {
        	title(paste(x$probability.model, "model using a", x$link, 
            "link", sep = " "))
    	}
    	else if (x$probability.model == "AdjCategory") {
        	title(paste(x$probability.model, "model using a loge link", 
            sep = " "))
    	}
    	else if (x$probability.model == "Stereotype") {
        	title(paste(x$probability.model, "model using a logit link", 
            sep = " "))
    	}
    } else title(main)
}
