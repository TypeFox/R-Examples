#
# WRITE ALPHA-NOMINATE FUNCTION
# "rcObject" is a roll call matrix written as a pscl object
#
#

summary.anominate <- function(object, ...){
        x <- object 
        if(!class(x)=="anominate") stop("Input is not of class 'anominate'.")
	
	cat("\n\nSUMMARY OF ANOMINATE OBJECT")
	cat("\n---------------------------\n")
	cat("\nNumber of Legislators:\t  ", dim(na.omit(x$wnom.result$legislators))[1],
	" (", dim(x$wnom.result$legislators)[1]-dim(na.omit(x$wnom.result$legislators))[1],
	" legislators deleted)", sep="")
	cat("\nNumber of Votes:\t  ", dim(na.omit(x$wnom.result$rollcalls))[1],
	" (", dim(x$wnom.result$rollcalls)[1]-dim(na.omit(x$wnom.result$rollcalls))[1],
	" votes deleted)", sep="")
	cat("\nNumber of Dimensions:\t ", x$wnom.result$dimensions)
	cat("\nalpha Mean:", round(mean(x$alpha),3))
		alpha.distribution <- round(quantile(x$alpha, probs=c(0.025,0.05,0.5,0.95,0.975)),3)
		alpha.distribution <- matrix(alpha.distribution,ncol=5)
		rownames(alpha.distribution) <- "     "
		colnames(alpha.distribution) <- c("2.5%","5%","50%","95%","97.5%")
	cat("\nalpha Percentiles:\n")
	print(alpha.distribution,...)
	cat("\n\n")

}



densplot.anominate <- function(x, ...){
	if(!class(x)=="anominate") stop("Input is not of class 'anominate'.")
        densplot(x$alpha, main=paste("Density plot of alpha\n Mean = ", round(mean(x$alpha),3),
		" [",round(quantile(x$alpha, probs=0.025),3),", ", round(quantile(x$alpha, probs=0.975),3),"]",sep=""), ... )
}


traceplot.anominate <- function(x, ...){
	
	if(!class(x)=="anominate") stop("Input is not of class 'anominate'.")
	
	traceplot(x$alpha, main=paste("Trace plot of alpha\n Mean = ", round(mean(x$alpha),3),
	" [",round(quantile(x$alpha, probs=0.025),3),", ",round(quantile(x$alpha, probs=0.975),3),"]",sep=""), ...)

}

plot.anominate <- function(x, ...){

	if(!class(x)=="anominate") stop("Input is not of class 'anominate'.")

	legislators <- x$legislators
	wnom.result <- x$wnom.result

	dims <- length(x$legislators)
	nlegis <- ncol(x$legislators[[1]])

	anom.means <- matrix(NA,nrow=nlegis,ncol=dims)
	for (i in 1:dims){
	anom.means[,i] <- summary(as.mcmc(legislators[[i]]))[[1]][,"Mean"]
	}

	wnom.coords <- na.omit(wnom.result$legislators[,grepl("coord",colnames(wnom.result$legislators))])
	
	if (dims==1){
	wnom.coords <- matrix(wnom.coords,nrow=nlegis,ncol=1)
	}
	
	wnom.coords <- na.omit(wnom.coords)

	dimension.name <- c("First","Second","Third","Fourth","Fifth","Sixth","Seventh","Eighth","Ninth","Tenth")

	if (dims == 1){
	dev.new(width=4, height=4)
	plot(wnom.coords[,i], anom.means[,i],
             main="W-NOMINATE vs. a-NOMINATE\nIdeal Points (w/ 95% CIs)", 
             xlab=paste("W-NOM (",dimension.name[i]," Dimension)",sep=""),
             ylab=paste("a-NOM (",dimension.name[i]," Dimension)",sep=""),
             type="n", bty="n", xlim=c(-1,1), ylim=c(floor(min(anom.means[,i],na.rm=TRUE)), ceiling(max(anom.means[,i],na.rm=TRUE))))
       	segments(wnom.coords[,i], summary(legislators[[i]])[[2]][,1],
		wnom.coords[,i], summary(legislators[[i]])[[2]][,5], col="grey", lwd=2)
        points(wnom.coords[,i], anom.means[,i], pch=20, cex=1.1)
        }

	if (dims == 2){
	dev.new(width=7, height=4)
	par(mfrow=c(1,2))
	for (i in 1:dims){
	plot(wnom.coords[,i], anom.means[,i],
	main="W-NOMINATE vs. a-NOMINATE\nIdeal Points (w/ 95% CIs)", 
	xlab=paste("W-NOM (",dimension.name[i]," Dimension)",sep=""),
	ylab=paste("a-NOM (",dimension.name[i]," Dimension)",sep=""),
	segments(wnom.coords[,i],summary(legislators[[i]])[[2]][,1],
		wnom.coords[,i],summary(legislators[[i]])[[2]][,5]))
	}
	}

	if (dims >= 3 & dims <= 4){
	dev.new(width=8, height=8)
	par(mfrow=c(2,2))
	for (i in 1:dims){
	plot(wnom.coords[,i], anom.means[,i],
	main="W-NOMINATE vs. a-NOMINATE\nIdeal Points (w/ 95% CIs)", 
	xlab=paste("W-NOM (",dimension.name[i]," Dimension)",sep=""),
	ylab=paste("a-NOM (",dimension.name[i]," Dimension)",sep=""),
	segments(wnom.coords[,i],summary(legislators[[i]])[[2]][,1],
		wnom.coords[,i],summary(legislators[[i]])[[2]][,5]))
	}
	}

	if (dims >= 5 & dims <= 9){
	dev.new(width=9, height=9)
	par(mfrow=c(3,3))
	for (i in 1:dims){
	plot(wnom.coords[,i], anom.means[,i],
	main="W-NOMINATE vs. a-NOMINATE\nIdeal Points (w/ 95% CIs)", 
	xlab=paste("W-NOM (",dimension.name[i]," Dimension)",sep=""),
	ylab=paste("a-NOM (",dimension.name[i]," Dimension)",sep=""),
	segments(wnom.coords[,i],summary(legislators[[i]])[[2]][,1],
		wnom.coords[,i],summary(legislators[[i]])[[2]][,5]))
	}
	}

}



anominate <- function(rcObject, dims=1, nsamp=1000, thin=1, burnin=500, 
	minvotes=20, lop=0.025, polarity=1,
	random.starts=TRUE, verbose=FALSE, constrain=FALSE) 
{
    cat("\nPreparing to run alpha-NOMINATE...\n\n")
    start <- proc.time()

    wnom.result <- wnominate(rcObject, ubeta = 15, uweights = 0.5,
                             dims = dims, minvotes = minvotes, 
                             lop = lop, trials = 3, polarity = polarity,
                             verbose = FALSE) 

    use.votes <- !is.na(wnom.result$rollcalls$spread1D)
    use.legis <- !is.na(wnom.result$legislators$coord1D)
    nvotes <- sum(use.votes)
    nlegis <- sum(use.legis)


    if(random.starts==TRUE){
      legis.starts <- runif(dims*nlegis,min=-1,max=1)
	for (i in 1:dims){
		legis.starts[((i-1)*nlegis)+polarity[i]] <- abs(legis.starts[((i-1)*nlegis)+polarity[i]])
		}
      bill.starts <- runif(2*dims*nvotes,min=-1,max=1)
    }

    if(random.starts==FALSE){
      legis.starts <- as.vector(t(t(wnom.result$legislators[,grepl("coord",colnames(wnom.result$legislators))])))
      legis.starts <- legis.starts[!is.na(legis.starts)]
      midpoints <- wnom.result$rollcalls[,grepl("midpoint",colnames(wnom.result$rollcalls))]
      spreads <- wnom.result$rollcalls[,grepl("spread",colnames(wnom.result$rollcalls))]
      yea.starts <- as.vector(t(t(midpoints - spreads)))
      yea.starts <- yea.starts[!is.na(yea.starts)]
      nay.starts <- as.vector(t(t(midpoints + spreads)))
      nay.starts <- nay.starts[!is.na(nay.starts)]
      bill.starts <- c(yea.starts,nay.starts)
    }

#    if(random.starts==FALSE){
#      legis.starts <- as.vector(t(wnom.result$legislators[,grepl("coord",colnames(wnom.result$legislators))]))
#      legis.starts <- legis.starts[!is.na(legis.starts)]
#      bill.midpoints <- as.matrix(wnom.result$rollcalls[,grepl("midpoint",colnames(wnom.result$rollcalls))])
#      bill.spreads <- as.matrix(wnom.result$rollcalls[,grepl("spread",colnames(wnom.result$rollcalls))])
#      D <- vector("list",2)
#      for(i in 1:dims){
#        D[[i]] <- cbind(bill.midpoints[,i] - bill.spreads[,i], bill.midpoints[,i] + bill.spreads[,i])
#      }
#      bill.starts <- as.vector(t(do.call(cbind,D)))
#      bill.starts <- bill.starts[!is.na(bill.starts)]
#    }



    if (!is.null(rcObject$codes)) {
        rcObject$votes[rcObject$votes %in% c(rcObject$codes$missing, 
            rcObject$codes$notInLegis)] <- -1
        rcObject$votes[rcObject$votes %in% rcObject$codes$yea] <- 1
        rcObject$votes[rcObject$votes %in% rcObject$codes$nay] <- 0
    }

  
    rcmatrix <- as.vector(rcObject$votes[use.legis,use.votes])
    nparam <- 2 + (dims*nlegis) + (dims*2*nvotes) + (dims*dims) + (dims*dims)
    output.len <- (floor(nsamp/thin)+1)*nparam

#   print(table(rcmatrix))
#   cat(sprintf("Number of votes = %5i\n", sum(rcmatrix != -1)))

#    print(summary(legis.starts))

if(verbose==TRUE){
	verbose <- 100
}

if (verbose==FALSE){	
	verbose <- 100000000
}

    anom.res <- .C("Canominate",
               rcdata = as.integer(rcmatrix),
               legis.starts = legis.starts,
               bill.starts = bill.starts,
               output=numeric(output.len),
               thin=as.integer(thin),
               ncol=as.integer(nvotes),
               nrow=as.integer(nlegis),
               nsamples=as.integer(nsamp),
               dim=as.integer(dims),
               verbose=as.integer(verbose),
	       constrain=as.integer(constrain))
    
     anom.mat <- matrix(anom.res$output,ncol=nparam,byrow=TRUE)
     anom.mat <- anom.mat[(burnin+1):nrow(anom.mat),]
     anom.mcmc <- as.mcmc(anom.mat[-nrow(anom.mat),])
     beta <- anom.mcmc[,1]
     alpha <- anom.mcmc[,2]
     #legislators <- anom.mcmc[,3:(2+nlegis)]
     #rollcall.yeas <- anom.mcmc[,(3+nlegis):(2+nlegis+nvotes)]
     #rollcall.nays <- anom.mcmc[,(3+nlegis+nvotes):(2+nlegis+nvotes+nvotes)]

legislators <- vector("list",dims)
	for(i in 1:dims){
	legislators[[i]] <- anom.mcmc[,(((i-1)*nlegis)+3):(2+(nlegis*i))]
	}

rollcalls <- anom.mcmc[,(3+(nlegis*dims)):(2+(nlegis*dims)+(2*nvotes*dims))]

all.yeas <- rollcalls[,1:(nvotes*dims)]
all.nays <- rollcalls[,(1+(nvotes*dims)):(2*nvotes*dims)]

rollcall.yeas <- vector("list",dims)
	for(i in 1:dims){
	rollcall.yeas[[i]] <- all.yeas[,(nvotes*(i-1)+1):(nvotes*i)]
	}

rollcall.nays <- vector("list",dims)
	for(i in 1:dims){
	rollcall.nays[[i]] <- all.nays[,(nvotes*(i-1)+1):(nvotes*i)]
	}


     legis.names <- rownames(wnom.result$legislators[use.legis,])
     vote.names <- paste("Vote",rownames(wnom.result$rollcalls[use.votes,]))

for (i in 1:dims){
colnames(legislators[[i]]) <- legis.names
colnames(rollcall.yeas[[i]]) <- vote.names
colnames(rollcall.nays[[i]]) <- vote.names
}

### ROTATION 

# Center each iteration at 0 by subtracting legislator means
for (i in 1:dims){
rollcall.yeas[[i]] <- sweep(rollcall.yeas[[i]], 1, rowMeans(legislators[[i]]), "-")
rollcall.nays[[i]] <- sweep(rollcall.nays[[i]], 1, rowMeans(legislators[[i]]), "-")
legislators[[i]] <- sweep(legislators[[i]], 1, rowMeans(legislators[[i]]), "-")
}


### THIS IS JUST A  CHECK TO MAKE SURE THE ROTATION IS WORKING
#legislators.old <- legislators
#rollcall.yeas.old <- rollcall.yeas
#rollcall.nays.old <- rollcall.nays
###


# Rotate it to the target W-NOMINATE
wnom.legis.coords <- wnom.result$legislators[,grepl("coord",colnames(wnom.result$legislators))]
wnom.legis.coords <- na.omit(wnom.legis.coords)
wnom.legis.coords <- as.matrix(wnom.legis.coords)

midpoints <- wnom.result$rollcalls[,grepl("midpoint",colnames(wnom.result$rollcalls))]
spreads <- wnom.result$rollcalls[,grepl("spread",colnames(wnom.result$rollcalls))]
wnom.yeas <- midpoints - spreads
wnom.yeas <- na.omit(wnom.yeas)
wnom.yeas <- as.matrix(wnom.yeas)
wnom.nays <- midpoints + spreads
wnom.nays <- na.omit(wnom.nays)
wnom.nays <- as.matrix(wnom.nays)

nchains <- nrow(legislators[[1]])

for (j in 1:nchains){
for (i in 1:dims){
	legislators[[i]][j,] <- procrustes(t(lapply(j, function(i) do.call(rbind, lapply(legislators, "[", i, TRUE)))[[1]]), as.matrix(wnom.legis.coords))$X.new[,i]
	rollcall.yeas[[i]][j,] <- procrustes(t(lapply(j, function(i) do.call(rbind, lapply(rollcall.yeas, "[", i, TRUE)))[[1]]), as.matrix(wnom.yeas))$X.new[,i]
	rollcall.nays[[i]][j,] <- procrustes(t(lapply(j, function(i) do.call(rbind, lapply(rollcall.nays, "[", i, TRUE)))[[1]]), as.matrix(wnom.nays))$X.new[,i]
}
}



    anomObject <- list(alpha = as.mcmc(alpha),
	beta = as.mcmc(beta),
	legislators = as.mcmc(legislators),
	yea.locations = as.mcmc(rollcall.yeas),
	nay.locations = as.mcmc(rollcall.nays),
	wnom.result = wnom.result)

    class(anomObject) <- c("anominate")

    cat("alpha-NOMINATE estimation completed successfully.")
    cat("\nalpha-NOMINATE took", (proc.time() - start)[3], "seconds to execute.\n\n")

    anomObject

}


anominate.sim <- function(nvotes=500, nlegis=101, seed=123345, utility="normal") 
{

  if(utility !="normal" | utility!="quadratic") cat("\n Utility must be 'normal' or 'quadratic'")
    
  set.seed(seed)
  yealocs <- matrix(runif(nvotes,min=-1.0,max=0.5),nrow=nvotes,ncol=1)
  naylocs <- matrix(runif(nvotes,min=-0.5,max=1.0),nrow=nvotes,ncol=1)
  ideallocs <- matrix(runif(nlegis,min=-0.9,max=0.9),nrow=nlegis,ncol=1)

  if(utility=="normal"){
  nom <- 1*(matrix(runif(nvotes*nlegis),nlegis,nvotes) < nomprob(yealocs, naylocs, ideallocs, 10, 0.5, normal=1))
  rcObject <- rollcall(nom)
  }

  if(utility=="quadratic"){
  qn <- 1*(matrix(runif(nvotes*nlegis),nlegis,nvotes) < qnprob(yealocs, naylocs, ideallocs, 3.25, 0.5, normal=1))
  rcObject <- rollcall(qn)
  }

  return(rcObject)
}
