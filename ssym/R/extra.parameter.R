extra.parameter <-
function(object, lower, upper, grid){
    new.c <- as.list(object$call)
	new.c$std.out <- FALSE
	if(missingArg(grid)) grid <- 1
	grid <- max(floor(abs(grid)),10)
    if(object$family!="Sinh-t" & object$family!="Contnormal"){
		xis <- seq(lower,upper,length=grid)
		conver <- matrix(0,grid,1)
		result <- matrix(0,grid,1)
		result2 <- matrix(0,grid,1)
		i <- 1
		bar <- txtProgressBar(min=1, max=grid, initial=0, width=50, char="+", style=3)
		while(i <= grid){
			new.c$xi <- xis[i]
			temp <- try(eval(as.call(new.c),envir = parent.frame()), silent=TRUE)
			if(is.list(temp)){
			  if(object$censored==FALSE) temp2 <- qqnorm(qnorm(temp$cdfz),plot.it=FALSE)
			  else{
			      surv0 <- survfit(Surv(temp$z_es,1-temp$event)~1)
				  ids <- ifelse(surv0$n.event>0,TRUE,FALSE)
				  survs <- ifelse(1-surv0$surv[ids] < 1e-30,1 - 1e-30, 1- surv0$surv[ids])
				  survs <- ifelse(survs > 1 - 1e-15,1 - 1e-15, survs)
				  probs <- temp$cdfz(surv0$time[ids])
			  	  probs <- ifelse(probs < 1e-30,1e-30, probs)
			  	  probs <- ifelse(probs > 1 - 1e-15,1 - 1e-15, probs)
			  	  temp2 <- list(x=qnorm(survs),y=qnorm(probs))
			  }
			  result[i] <- mean(abs(sort(temp2$x)-sort(temp2$y)))
			  result2[i] <- -2*sum(temp$lpdf)
			  conver[i] <- 1
			}
			i <- i + 1
            setTxtProgressBar(bar,i)
		}
		close(bar)
        cat("\n")
		xis <- xis[conver==1]
		result <- as.matrix(result[conver==1])
		par(mfrow=c(1,2))
		plot(xis,result,type="b",xlim=range(xis),ylim=range(result),xlab=expression(eta),ylab=expression(Upsilon(eta)))
		title("Behaviour of the overall goodness-of-fit statistic")
		result2 <- as.matrix(result2[conver==1])
		plot(xis,result2,type="b",xlim=range(xis),ylim=range(result2),xlab=expression(eta),ylab="-2*log-Likelihood")
		title("Behaviour of -2*log-Likelihood")

	}else{
		xis1 <- seq(lower[1],upper[1],length=grid)
		xis2 <- seq(lower[2],upper[2],length=3)
		conver <- matrix(0,grid,3)
		result <- matrix(0,grid,3)
		result2 <- matrix(0,grid,3)
		bar <- txtProgressBar(min=1, max=3*grid, initial=0, width=50, char="+", style=3)
		for(i in 1:3){
			j <- 1
			while(j <= grid){
				new.c$xi <- c(xis1[j],xis2[i])
				temp <- try(eval(as.call(new.c),envir = parent.frame()), silent=TRUE)
				if(is.list(temp)){
				  if(object$censored==FALSE) temp2 <- qqnorm(qnorm(temp$cdfz),plot.it=FALSE)
				  else{surv0 <- survfit(Surv(temp$z_es,1-temp$event)~1)
	  			  ids <- ifelse(surv0$n.event>0,TRUE,FALSE)
				  survs <- ifelse(1-surv0$surv[ids] < 1e-30,1 - 1e-30, 1- surv0$surv[ids])
				  survs <- ifelse(survs > 1 - 1e-15,1 - 1e-15, survs)
				  probs <- temp$cdfz(surv0$time[ids])
				  probs <- ifelse(probs < 1e-30,1e-30, probs)
				  probs <- ifelse(probs > 1 - 1e-15,1 - 1e-15, probs)
				  temp2 <- list(x=qnorm(survs),y=qnorm(probs))}
				  result[j,i] <- mean(abs(sort(temp2$x)-sort(temp2$y)))
				  result2[j,i] <- -2*sum(temp$lpdf)
				  conver[j,i] <- 1
				}
				j <- j + 1
                setTxtProgressBar(bar,(i-1)*grid + j)
			}
		}
		close(bar)
		cat("\n")
		par(mfrow=c(1,2))
		ylim <- range(result[conver[,1]==1,1],result[conver[,2]==1,2],result[conver[,3]==1,3])
		xlim <- range(xis1[conver[,1]==1],xis1[conver[,2]==1],xis1[conver[,3]==1])
		plot(xis1[conver[,1]==1],result[conver[,1]==1,1],type="b",xlim=xlim,ylim=ylim,xlab="",ylab="")
		par(new=TRUE)
		plot(xis1[conver[,2]==1],result[conver[,2]==1,2],type="b",xlim=xlim,ylim=ylim,xlab="",ylab="",col="red")
		par(new=TRUE)
		plot(xis1[conver[,3]==1],result[conver[,3]==1,3],type="b",xlim=xlim,ylim=ylim,xlab=expression(eta[1]),ylab=expression(Upsilon(eta)),col="blue")
		legend(xlim[1],ylim[2],lty=1,col=c("black","red","blue"),title=expression(eta[2]),legend=c(xis2[1],xis2[2],xis2[3]))
		title("Behaviour of the overall goodness-of-fit statistic")

		ylim <- range(result2[conver[,1]==1,1],result2[conver[,2]==1,2],result2[conver[,3]==1,3])
		xlim <- range(xis1[conver[,1]==1],xis1[conver[,2]==1],xis1[conver[,3]==1])
		plot(xis1[conver[,1]==1],result2[conver[,1]==1,1],type="b",xlim=xlim,ylim=ylim,xlab="",ylab="")
		par(new=TRUE)
		plot(xis1[conver[,2]==1],result2[conver[,2]==1,2],type="b",xlim=xlim,ylim=ylim,xlab="",ylab="",col="red")
		par(new=TRUE)
		plot(xis1[conver[,3]==1],result2[conver[,3]==1,3],type="b",xlim=xlim,ylim=ylim,xlab=expression(eta[1]),ylab="AIC",col="blue")
		legend(xlim[1],ylim[2],lty=1,col=c("black","red","blue"),title=expression(eta[2]),legend=c(xis2[1],xis2[2],xis2[3]))
		title("Behaviour of -2*log-Likelihood")
	}
}
