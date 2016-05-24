"check.trait" <-
function (trait,data,fdrate=0.05,graph=TRUE,binshow=FALSE,qoption="bh95") {
  if (!is.data.frame(data)) {
    if (is(data,"gwaa.data")) data <- data@phdata else stop("data should be a data frame or gwaa.data")
  } 
#  attach(data,warn.conflicts=FALSE)
	with(data, {
  rmvec <- rep(TRUE,length(trait))
  for (i in 1:length(trait)) 
    if (!is.numeric(get(trait[i]))) {
    		cat("--------------------------------\n")
		cat("trait",trait[i],"is not numeric, skipping\n");
		rmvec[i] = FALSE
	}
  trait <- trait[rmvec]
  for (i in 1:length(trait)) {
    tra <- get(trait[i])
    mu <- mean(tra,na.rm=TRUE)
    rvar <- 1./var(tra,na.rm=TRUE)
    Pv <- rep(1,length(tra))
    if (length(levels(as.factor(tra)))>2) 
     for (j in 1:length(tra)) {
      if (is.na(tra[j]) || is.na(rvar)) {Pv[j] = 1.0} else
      { if (j!=length(tra) && j!=1)
		tt <- tra[c(seq(1:(j-1)),seq((j+1),length(tra)))]
	else if (j==1)
		tt <- tra[seq(2:length(tra))]
	else if (j==length(tra))
		tt <- tra[seq(1:(length(tra)-1))]
	muj <- mean(tt,na.rm=TRUE)
	rvarj <- 1./var(tt,na.rm=TRUE)
	Pv[j] = pchisq((tra[j]-muj)*(tra[j]-muj)*rvarj,1,lower.tail = FALSE)
	}
     }
    cat("--------------------------------\n")
    cat("Trait",trait[i],"has",sum(!is.na(get(trait[i]))),"measurements\n");
    cat("Missing:",sum(is.na(get(trait[i]))),"(",100*sum(is.na(get(trait[i])))/length(get(trait[i])),"%)\n");
    cat("Mean =",mean(get(trait[i]),na.rm=TRUE),"; s.d. =",sqrt(var(get(trait[i]),na.rm=TRUE)),"\n");
    if (min(Pv,na.rm=TRUE)<fdrate) {
      qobj<-qvaluebh95(Pv,fdrate = fdrate)
      if (sum(qobj$signif)>=1) {
        cat("Outliers discovered for trait",trait[i],":\n");
        cat(data$id[qobj$signif])
        cat("\n")
        cat(tra[qobj$signif])
        cat("\n")
      } else {
        cat("NO outliers discovered for trait",trait[i],"\n");
      }
    } else {
      cat("NO outliers discovered for trait",trait[i],"\n");
    }
  }
  if (graph) {
  	rmvec <- rep(TRUE,length(trait))
  	if (!binshow) for (i in 1:length(trait)) {
		if (length(levels(as.factor(get(trait[i]))))<=2) {
			rmvec[i] = FALSE
		}
	}
	trait <- trait[rmvec]
  	nt <- length(trait)
	spar <- par(no.readonly=TRUE)
	if (nt>1) {
		par(mfrow=c(nt,nt))
	  	for (i in 1:nt) for (j in 1:nt) {
			plot(get(trait[i]),get(trait[j]),xlab=trait[i],ylab=trait[j])
		}
	} else if (nt==1) {
		hist(get(trait[1]),xlab=trait[1],main=paste("Histogram of",trait[1]))
	}
	par(spar)
  }
  })
#  detach(data)
}

