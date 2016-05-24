predict.mixor <-
function(object, newdata=NULL, na.action=na.fail, ...) {
    y <- object$Y
    link <- object$link
    IADD<-object$IADD
    weights<-object$weights
    KG <- object$KG
    KS <-object$KS
    k <- object$MAXJ
    ICEN <- object$ICEN
    which.random.slope<-object$which.random.slope
    if (IADD==1) {
	   	if (object$random.effect.mean==TRUE & ICEN==0) {
    		alpha <- c(as.numeric(object$MU[1]), as.numeric(object$MU[1])+as.numeric(object$GAM)[1:(k-2)])
    	} else if (object$random.effect.mean==TRUE & ICEN==1) {
      		alpha <- c(as.numeric(object$MU[1]), as.numeric(object$MU[1])+as.numeric(object$GAM)[1:(k-1)]) 	
    	} else {
    		alpha <- as.numeric(object$GAM)   	
    	}
    } else {
	   	if (object$random.effect.mean==TRUE & ICEN==0) {
    		alpha <- c(as.numeric(object$MU[1]), as.numeric(object$MU[1])-as.numeric(object$GAM)[1:(k-2)])
    	} else if (object$random.effect.mean==TRUE & ICEN==1) {
      		alpha <- c(as.numeric(object$MU[1]), as.numeric(object$MU[1])-as.numeric(object$GAM)[1:(k-1)]) 	
    	} else {
    		alpha <- -as.numeric(object$GAM)   	
    	}   
    }
    tt <- terms(object)
    if (missing(newdata) || is.null(newdata)) {
        neww <- object$W
        newx <- object$X
    } else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action, #na.omit
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, m)
		mf <- object$call
        mfmatch <- match(c("formula", "data", "id", "subset", "weights"), 
	        names(mf), 0L)
        IDS <- newdata[,as.character(mf[["id"]])]
        mmDone <- FALSE    
        neww <- as.matrix(m)
#       newx <- as.matrix(newx) # No random effects estimates for newdata!
    }
    n <- dim(neww)[1]
    levels <- sort(unique(y))
	zeta<-object$ALPHA
	reps<-table(object$id)
	if (missing(newdata) || is.null(newdata)) {
		if (dim(object$EBmean)[2]==1) {
			if (dim(newx)[2]==1  ) {                 #### does ( | is.numeric(newx) ) need to be added?
				theta<-as.matrix(rep(object$EBmean,reps), ncol=1)
				Xb <- neww %*% zeta + newx * theta
			} else {
				Xb<-matrix(0, ncol=1, nrow=dim(neww)[1])
				dimnames(Xb)[[1]]<-object$id
				unique.id<-unique(object$id)
				for (i in 1:length(unique.id)) {	  
					dimnames(newx)[[1]]<-object$id
					dimnames(neww)[[1]]<-object$id
					x<-newx[dimnames(newx)[[1]]==unique.id[i],]
					w<-neww[dimnames(neww)[[1]]==unique.id[i],]
					if (!is.null(which.random.slope) ) {
						theta<-rep(object$EBmean[i,], length(which.random.slope))			
					} else {
						theta<-rep(object$EBmean[i,], sum(dimnames(newx)[[1]]==unique.id[i]))
					}
					Xb[dimnames(neww)[[1]]==unique.id[i],] <- w %*% zeta + x %*% theta 
				}
	  		} 
		} else {
			Xb<-matrix(0, ncol=1, nrow=dim(neww)[1])
			dimnames(Xb)[[1]]<-object$id
			for (i in 1:length(unique(object$id))) {
				unique.id<-unique(object$id)
				dimnames(newx)[[1]]<-object$id
				dimnames(neww)[[1]]<-object$id
				x<-newx[dimnames(newx)[[1]]==unique.id[i],]
				w<-neww[dimnames(neww)[[1]]==unique.id[i],]
				theta<-object$EBmean[i,]
				Xb[dimnames(neww)[[1]]==unique.id[i],] <- w %*% zeta + x %*% theta 
			}
		}
	} else {
		Xb <- neww %*% zeta
	}
	if (ICEN==1) k=k+1
    z <- matrix(ncol = k - 1, nrow = n)
    for (i in 1:(k - 1)) {
        if (IADD==1) {
        	if (length(alpha)==0L & KG==0 & KS==0) {
        		z[, i] <- Xb
        	} else if (KG==0 & KS==0) {
            	z[, i] <- alpha[i] + Xb
            } else if (KG==0 & KS>0) {
				scale<-object$TAU
				if (KS>=2) { ###KS scaling parameter MONDAY! ###
					z[, i] <- (alpha[i] + Xb)/exp(neww[,1:KS]%*%scale)
				} else {
					z[, i] <- (alpha[i] + Xb)/exp(neww[,1:KS]*scale)
				}
			} else if (KG>0 & KS==0) {
           		new.thresh<-c(0,as.numeric(object$GAM)[-(1:(k-2))])
            	if (KG>=2) {
            		interact <- apply(new.thresh[seq(i,length(new.thresh),by=2)]*neww[,1:KG], 1, sum)
            	} else {
            		interact <- new.thresh[i]*neww[,1:KG]
            	}
            	z[, i] <- alpha[i] - interact + Xb      
            }
        } else {
	        if (length(alpha)==0L & KG==0 & KS==0) {
        		z[, i] <- -Xb
        	} else if (KG==0 & KS==0) {
            	z[, i] <- -alpha[i] - Xb
            } else  if (KG==0 & KS>0) {
				scale<-object$TAU
				if (KS>=2) { 
					z[, i] <- -(alpha[i] + Xb)/exp(neww[,1:KS]%*%scale)
				} else {
					z[, i] <- -(alpha[i] + Xb)/exp(neww[,1:KS]*scale)
				}
			} else if (KG>0 & KS==0) {
           		new.thresh<-c(0, as.numeric(object$GAM)[-(1:(k-2))])
            	if (KG>=2) {
            		interact<- apply(new.thresh[seq(i,length(new.thresh),by=2)]*neww[,1:KG], 1, sum)
            	} else {
            		interact<-new.thresh[i]*neww[,1:KG]
            	}
            	z[, i] <- -alpha[i] + interact - Xb               
            }
        }
    }
    pi <- matrix(ncol = k, nrow = n)
    if (missing(newdata) || is.null(newdata)) {
    	dimnames(pi)[[1]]<-object$id
    } else {
    	dimnames(pi)[[1]]<-IDS
    }
    if (ICEN==0) {
    	dimnames(pi)[[2]]<-sort(unique(y))
    } else {
    	dimnames(pi)[[2]]<-c(sort(unique(y)),paste(sort(unique(y))[k-1],"+",sep=""))   
    }
        if (link == "logit") {
            for (i in 1:k) {
                if (i == 1) {
                  pi[, i] <- exp(z[, i])/(1 + exp(z[, i]))
                }
                else if (i <= k - 1) {
                  pi[, i] <- exp(z[, i])/(1 + exp(z[, i])) - 
                    exp(z[, i - 1])/(1 + exp(z[, i - 1]))
                }
                else if (i == k) {
                  pi[, i] <- 1 - exp(z[, i - 1])/(1 + exp(z[, 
                    i - 1]))
                }
            }
        } else if (link == "probit") {
            for (i in 1:k) {
                if (i == 1) {
                  pi[, i] <- pnorm(z[, i])
                }
                else if (i <= k - 1) {
                  pi[, i] <- pnorm(z[, i]) - pnorm(z[, i - 1])
                }
                else if (i == k) {
                  pi[, i] <- 1 - pnorm(z[, i - 1])
                }
            }
        } else if (link == "cloglog") {
            for (i in 1:k) {
                if (i == 1) {
                  pi[, i] <- 1 - exp(-exp(z[, i]))
                }
                else if (i <= k - 1) {
                  pi[, i] <- exp(-exp(z[, i - 1])) - exp(-exp(z[, 
                    i]))
                }
                else if (i == k) {
                  pi[, i] <- exp(-exp(z[, i - 1]))
                }
            }
        }  else if (link == "loglog") {
            for (i in 1:k) {
                if (i == 1) {
                  pi[, i] <- exp(-exp(z[, i]))
                }
                else if (i <= k - 1) {
                  pi[, i] <- exp(-exp(z[, i ])) - exp(-exp(z[, 
                    i-1]))
                }
                else if (i == k) {
                  pi[, i] <- 1 - exp(-exp(z[, i - 1]))
                }
            }        
        }
    levels<-dimnames(pi)[[2]]
    class <- levels[apply(pi, 1, which.max)]
    list(predicted = pi, class = class)
}
