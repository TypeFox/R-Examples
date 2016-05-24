ordinal.gmifs <-
function(formula, data, x=NULL, subset, epsilon=0.001, tol=1e-5, scale=TRUE, 
                                  probability.model="Cumulative",
								  link="logit", verbose=FALSE, ...){
    mf <- match.call(expand.dots = FALSE)
    cl <- match.call()
    m <- match(c("formula", "data", "subset"), names(mf), 0L)    
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf)
    w <- model.matrix(mt, mf)
    if (!is.null(x)) {
    	if (missing(subset)) 
        	r <- TRUE 
    	else {  
        	e <- substitute( subset) 
         	r <- eval( e, data) 
         	if (!is.logical(r))  
             	stop("'subset' must evaluate to logical" ) 
         	r <- r & !is.na(r)  
     	} 
        if (class( x)=="character") {
            nl <- as.list( 1:ncol(data))  
            names(nl) <- names( data) 
            vars <- eval(substitute( x), nl, parent.frame())  
         	x <- data [r , vars, drop=FALSE ]
            x <- as.matrix(x )
        } else if (class(x)== "matrix" || class(x)== "data.frame") {
                x <- x[r,, drop =FALSE]
                x <- as.matrix(x)
        }
    }
    if (is.na(match(probability.model,c("Cumulative", "AdjCategory", "ForwardCR", "BackwardCR", "Stereotype")))) {
    	stop("Error:")
    	cat("Only Cumulative, AdjCategory, ForwardCR, BackwardCR, Stereotype available for probability.model parameter.\n")
    }
    if (probability.model=="Stereotype" & link!="logit") {
    	warning("Warning: Stereotype model only uses a logit link.")
    }
    if (probability.model=="AdjCategory" & link!="loge") {
    	warning("Warning: AdjCategory model only uses a loge link.")
    }
    if ( (probability.model=="Cumulative" | probability.model=="ForwardCR" | probability.model=="BackwardCR") &
    	is.na(match(link,c("logit","probit","cloglog")))) {
    	stop("Error: ")
    	cat("Only logit, probit, and cloglog links available for ", probability.model, " model.\n")
    }
    levels <- sort(unique(y))
	k <- length(unique(y))
    data <- data.matrix(data)
    is.intercept<-grep("Intercept",dimnames(w)[[2]])
    if (length(is.intercept)==1) {
    	w<-w[,-is.intercept, drop=FALSE]
    }
    zeta <- rep(0, dim(w)[2])
    Ymat <- matrix(0, nrow = length(y), ncol = k)
    for (i in levels) {
        Ymat[which(y == i), which(levels == i)] <- 1
    }
	alpha <- numeric()
    pi.0 <- table(y) / length(y)
    tab <- table(y)
	if (probability.model=="Cumulative") {
		if (link=="logit") {
			alpha <- log(cumsum(pi.0) / (1 - cumsum(pi.0)))[1:(k - 1)]
		} else if (link=="probit") {
		    alpha <- qnorm(cumsum(pi.0))[1:(k - 1)]
		} else if (link=="cloglog") {
			alpha <- log(-log(1-cumsum(pi.0)))[1:(k - 1)]
		}
	} else if (probability.model=="AdjCategory") {
		for(i in 1:(k-1)) {
			alpha[i] <- log(tab[i+1] / sum(tab[i:(i+1)])/(1-tab[i+1] / sum(tab[i:(i+1)])))
		}	
	} else if (probability.model=="ForwardCR") {
		Cum.Ymat <- matrix(0, nrow=nrow(Ymat), ncol=k-1)
		for(i in 1:(k-1)) {
			if (link=="logit") {
				alpha[i] <- log(tab[i] / sum(tab[i:k]) / (1-tab[i]/sum(tab[i:k])))
				Cum.Ymat[, i] <- rowSums(Ymat[, i:k])
			} else if (link=="probit") {
 		   		alpha[i]<- qnorm(tab[i]/sum(tab[i:k]))
				Cum.Ymat[,i]<-rowSums(Ymat[,i:k])  #		
			} else if (link=="cloglog") {
				alpha[i]<- log(-log(1-(tab[i]/sum(tab[i:k]))))
				Cum.Ymat[,i]<-rowSums(Ymat[,i:k])
			}
		}
	} else if (probability.model=="BackwardCR") {
   		Cum.Ymat <- matrix(0, nrow=nrow(Ymat), ncol=k-1)
		for(i in 1:(k-1)) {
			if (link=="logit") {
				alpha[i] <- log(tab[i+1] / 
						sum(tab[1:(i+1)]) / 
						(1 - tab[i+1] / sum(tab[1:(i+1)])))
				#  swap in some C code?
				Cum.Ymat[,i] <- apply(matrix(Ymat[,(i+1):k], 
					nrow=dim(Ymat)[1]), 1, sum)
			} else if (link=="probit") {
			    alpha[i]<- qnorm(tab[i+1]/sum(tab[1:(i+1)]))
				Cum.Ymat[,i] <- apply(matrix(Ymat[,(i+1):k], 
					nrow=dim(Ymat)[1]), 1, sum)
			} else if (link=="cloglog") {
				alpha[i]<- log(-log(1-(tab[i+1]/sum(tab[1:(i+1)]))))
				Cum.Ymat[,i] <- apply(matrix(Ymat[,(i+1):k], 
					nrow=dim(Ymat)[1]), 1, sum)
			}
		}
	} else if (probability.model=="Stereotype") {
	  for(i in 1:(k-1)) {
      	alpha[i]<- log(pi.0[i]/pi.0[k]) 
    	}
      phi <- rep(0.1,k-2)          
      phi.update <- matrix(phi, ncol = length(phi))
	}
	names(alpha) <- paste("(Intercept)", 1:(k-1), sep=":")
	ui<-matrix(0,nrow=length(alpha)-1,ncol=length(alpha)+length(zeta))  #constraint matrix
	for (i in 1:dim(ui)[1]) {
		ui[i,i]<- -1
		ui[i,i+1]<- 1
	}
    ci<-rep(0,dim(ui)[1])
	if (!is.null(x)) {
		vars <- dim(x)[2]
		oldx <- x
    	if (scale) {
        	sd <- apply(x, 2, sd)
        	for (i in 1:dim(x)[2]) {
        		if (sd[i]==0) {
            		x[,i] <- scale(x[,i], center = TRUE, scale = FALSE)
            	} else {
            		x[,i] <- scale(x[,i], center = TRUE, scale = TRUE)
            	}
        	}
        }
    	x <- cbind(x, -1 * x)
    	beta <- rep(0, dim(x)[2])
    	names(beta) <- dimnames(x)[[2]]
    	step <- 0
    	Estimates <- matrix(0,ncol=vars)
    	alpha.update <- matrix(alpha, ncol = k-1)
    	logLikelihood <- numeric()
    	#initialize zeta and re-update alpha if dim(w)[2]!=0
    	if (probability.model=="Cumulative") {
        	if (dim(w)[2]!=0) {
	  	        initialize <- optim(c(alpha,zeta), fn.cum, w=w, x=x, beta=beta, y=y, k=k, levels=levels, Ymat=Ymat, link=link, method="BFGS")
    			alpha <- initialize$par[1:(k-1)]
		    	alpha.update <- matrix(alpha, ncol = k-1)
    			zeta <- initialize$par[k:length(initialize$par)]
 				zeta.update <- matrix(zeta, ncol = length(zeta))         
			} 		
		} else if (probability.model=="AdjCategory") {
        	if (dim(w)[2]!=0) {
    	        initialize <- optim(c(alpha,zeta), fn.acat, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, method="BFGS")
    			alpha <- initialize$par[1:(k-1)]
		    	alpha.update <- matrix(alpha, ncol = k-1)
    			zeta <- initialize$par[k:length(initialize$par)]
 				zeta.update <- matrix(zeta, ncol = length(zeta))         
			} 		
		} else if (probability.model=="ForwardCR") {
        	if (dim(w)[2]!=0) {
    	       	initialize <- optim(c(alpha,zeta), fn.fcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    			alpha <- initialize$par[1:(k-1)]
			    alpha.update <- matrix(alpha, ncol = k-1)
    			zeta <- initialize$par[k:length(initialize$par)]
				zeta.update <- matrix(zeta, ncol = length(zeta))
 			}
		} else if (probability.model=="BackwardCR") {
        	if (dim(w)[2]!=0) {
    	       	initialize <- optim(c(alpha,zeta), fn.bcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    			alpha <- initialize$par[1:(k-1)]
			    alpha.update <- matrix(alpha, ncol = k-1)
    			zeta <- initialize$par[k:length(initialize$par)]
				zeta.update <- matrix(zeta, ncol = length(zeta))
 			}
		} else if (probability.model=="Stereotype") {
			if (dim(w)[2]!=0) {           
			initialize <- optim(c(alpha, phi, zeta), fn.stereo, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, 
				method="L-BFGS-B",upper=c(rep(Inf,k-1),rep(1,k-2)),lower=c(rep(-Inf,k-1),rep(0,k-2)))
			alpha <- initialize$par[1:(k-1)]
			alpha.update <- matrix(alpha, ncol=k-1)
			phi <- initialize$par[k:(2*k-3)]
			phi.update <- matrix(phi, ncol = length(phi))
			zeta <- initialize$par[(2*k-2):length(initialize$par)]
			zeta.update <- matrix(zeta, ncol = length(zeta)) 
			}
		}	
    	repeat {
			# Z is a k-1 column matrix filled with 
			# alpha_i + x*beta where alpha_i is the ith threshold
        	z <- matrix(ncol = k - 1, nrow = length(y))
        	if (dim(w)[2]!=0) {
          		xbeta <- w%*%zeta + x %*% beta
			} else {
				xbeta <- x%*%beta
			}
        	for (i in 1:(k - 1)) {
            	z[, i] <- alpha[i] + xbeta
        	}
			# u is the matrix that stores elements of the 
			# first derivative for each of the variables
			if (probability.model=="Cumulative") {
				update.j <- du.cum(z=z, Ymat=Ymat, k=k, x=x, link=link)
			} else if (probability.model=="AdjCategory") {
				update.j <- du.adjcat(eta=z, x=x, k=k, Ymat=Ymat)
			} else if (probability.model=="ForwardCR") {
				update.j <- du.fcr(eta=z, Ymat=Ymat, k=k, Cum.Ymat=Cum.Ymat, x=x, link=link)
			} else if (probability.model=="BackwardCR") {
				update.j <- du.bcr(eta=z, Ymat=Ymat, k=k, Cum.Ymat=Cum.Ymat, x=x, link=link)
			} else if (probability.model=="Stereotype") {
				update.j <- du.stereo(Xb=xbeta, x=x, alpha=alpha, phi=phi, Ymat=Ymat, k=k)
			}
    	    if (update.j$update.value < 0) {
				# update the identified variables beta by adding epsilon
            	beta[update.j$update.j] <- beta[update.j$update.j] + epsilon
        	}
        	Estimates <- rbind(Estimates, beta[1:vars] - beta[(vars+1):length(beta)])
	    	if (probability.model=="Cumulative") {
    			if (dim(w)[2]!=0) {
     				out <- constrOptim(theta=c(alpha,zeta), f=fn.cum, grad=gradient, ui=ui, ci=ci, w=w, x=x[,1:vars], beta=beta[1:vars]-beta[(vars+1):length(beta)], y=y, k=k, levels=levels, Ymat=Ymat, link=link, method="BFGS")
    				alpha<-out$par[1:(k-1)]
					alpha.update <- rbind(alpha.update, alpha)
    				zeta<-out$par[k:length(out$par)]
					zeta.update <- rbind(zeta.update, zeta)
				} else {
     				out <- constrOptim(theta=alpha, f=fn.cum, grad=gradient, ui=ui, ci=ci, w=w, x=x[,1:vars], beta=beta[1:vars]-beta[(vars+1):length(beta)], y=y, k=k, levels=levels, Ymat=Ymat, link=link, method="BFGS")
    				alpha<-out$par[1:(k-1)]
					alpha.update <- rbind(alpha.update, alpha)		
				}
			} else if (probability.model=="AdjCategory") {
        		if (dim(w)[2]!=0) {
    	        	out <- optim(par=c(alpha,zeta), fn=fn.acat, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, method="BFGS")
    				alpha <- out$par[1:(k-1)]
 					alpha.update <- rbind(alpha.update, alpha)
    				zeta <- out$par[k:length(out$par)]
					zeta.update <- rbind(zeta.update, zeta)
 				} else {
    				out <- optim(par=alpha, fn=fn.acat, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, method="BFGS")
    				alpha<-out$par[1:(k-1)]
					alpha.update <- rbind(alpha.update, alpha)		
 				}      
			} else if (probability.model=="ForwardCR") {
        		if (dim(w)[2]!=0) {
     	        	out <- optim(par=c(alpha,zeta), fn=fn.fcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    				alpha <- out$par[1:(k-1)]
 					alpha.update <- rbind(alpha.update, alpha)
    				zeta <- out$par[k:length(out$par)]
					zeta.update <- rbind(zeta.update, zeta)
 				} else {
    				out <- optim(par=alpha, fn=fn.fcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    				alpha<-out$par[1:(k-1)]
					alpha.update <- rbind(alpha.update, alpha)		
 				}      
			} else if (probability.model=="BackwardCR") {
        		if (dim(w)[2]!=0) {
    	        	out <- optim(par=c(alpha,zeta), fn=fn.bcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    				alpha <- out$par[1:(k-1)]
 					alpha.update <- rbind(alpha.update, alpha)
    				zeta <- out$par[k:length(out$par)]
					zeta.update <- rbind(zeta.update, zeta)
 				} else {
    				out <- optim(par=alpha, fn=fn.bcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    				alpha<-out$par[1:(k-1)]
					alpha.update <- rbind(alpha.update, alpha)		
 				}      
			} else if (probability.model=="Stereotype") {
				if (dim(w)[2]!=0) {
    				out<- optim(par=c(alpha,phi,zeta), fn=fn.stereo, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, 
    					method="L-BFGS-B",upper=c(rep(Inf,k-1),rep(1,k-2)),lower=c(rep(-Inf,k-1),rep(0,k-2)))
    				alpha <- out$par[1:(k-1)]
 					alpha.update <- rbind(alpha.update, alpha)
			    	phi <- out$par[k:(2*k-3)]
    				phi.update <- rbind(phi.update, phi)
        			zeta<-out$par[(2*k-2):length(out$par)] ### May need to differ for stereotype logit!
        			zeta.update <- rbind(zeta.update, zeta)
	        	} else {
    				out<- optim(c(alpha,phi), fn.stereo, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, 
    					method="L-BFGS-B",upper=c(rep(Inf,k-1),rep(1,k-2)),lower=c(rep(-Inf,k-1),rep(0,k-2)))
    				alpha <- out$par[1:(k-1)]
 					alpha.update <- rbind(alpha.update, alpha)
			    	phi <- out$par[k:(2*k-3)]
    				phi.update <- rbind(phi.update, phi)
	        	}
			}
			# extract the logLikelihood from the
			# optimization step
        	logLikelihood[step+1]<- LL1<- -out$value
			if (verbose)  cat("step = ", step, "\n")
			# if the difference between two successive 
			#  log-logLikelihoods is less than tol, stop updating
        	if (step >= 1 && LL1 - LL0 < tol) {
            	break
        	}
			# assign the current log-logLikelihood to LL0 so 
			# in the next step the difference can be taken
        	LL0 <- LL1
			# increment the step number
        	step <- 1 + step
    	}
	    beta  <- Estimates[-1,]
		# remove the initial alpha terms
		alpha <- alpha.update[-1,]
    	if (probability.model=="Stereotype") {
    		if (dim(w)[2]!=0) {
 			    p <- apply(beta, 1, function(x) length(x[x!=0]) ) + length(alpha) + +length(phi) + length(zeta)
    		} else {
    			p <- apply(beta, 1, function(x) length(x[x!=0]) ) + length(alpha) + length(phi)
    		}
		} else if (dim(w)[2]!=0) {
			p <- apply(beta, 1, function(x) length(x[x!=0]) ) + length(alpha) + length(zeta)
		} else {
			p <- apply(beta, 1, function(x) length(x[x!=0]) ) + length(alpha) 
		}
     	AIC<-2*p-2*logLikelihood
       	BIC<-p*(log(length(y))-log(2*pi))-2*logLikelihood
		# returns which row has minimum AIC for easy
		# model extraction
		model.select <- which.min(AIC)
 		if (dim(w)[2]!=0) {
			zeta <- zeta.update[-1,]
			if (class(zeta)!="matrix") {
				zeta<-matrix(zeta, ncol=dim(w)[2])
			}
			colnames(zeta)<-colnames(w)
		}
	    output <- list(beta=beta, alpha=alpha, zeta=zeta, x=oldx, y=y, w=w, scale=scale,
			logLik=logLikelihood, AIC=AIC, BIC=BIC, model.select=model.select, probability.model=probability.model, link=link)
		if (probability.model=="Stereotype") {
			phi.update <- cbind(1, phi.update)
			phi <- phi.update[-1,]
			colnames(phi)<- paste("phi",1:(k-1),sep=".")
			output <- list(beta=beta, alpha=alpha, phi=phi, zeta=zeta, x=oldx, y=y, w=w, scale=scale,
				logLik=logLikelihood, AIC=AIC, BIC=BIC, model.select=model.select, probability.model=probability.model, link=link)
			} 
    } else {
    	beta <- NULL
	    if (probability.model=="Cumulative") {   	
    		if (dim(w)[2]!=0) {
	  	        out <- optim(c(alpha,zeta), fn.cum, w=w, x=x, beta=beta, y=y, k=k, levels=levels, Ymat=Ymat, link=link, method="BFGS")
    			alpha<-out$par[1:(k-1)]
    			zeta<-out$par[k:length(out$par)]
    			names(zeta)=colnames(w)
    			logLik <- -out$value
    			output<-list(alpha=alpha, zeta=zeta, y=y, w=w, x=x, probability.model=probability.model, link=link, logLik=logLik)
    		} else {
	  	        out <- optim(alpha, fn.cum, w=w, x=x, beta=beta, y=y, k=k, levels=levels, Ymat=Ymat, link=link, method="BFGS")
    			alpha<-out$par[1:(k-1)]
    			logLik <- -out$value
    	    	output<-list(alpha=alpha, y=y, w=w, x=x, probability.model=probability.model, link=link, logLik=logLik)
    		}
		} else if (probability.model=="AdjCategory") {
        	if (dim(w)[2]!=0) {
    	       	out <- optim(par=c(alpha,zeta), fn=fn.acat, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, method="BFGS")
    			alpha <- out$par[1:(k-1)]
    			zeta <- out$par[k:length(out$par)]
    			names(zeta)=colnames(w)
    			logLik <- -out$value
    			output<-list(alpha=alpha, zeta=zeta, y=y, w=w, x=x, probability.model=probability.model, link=link, logLik=logLik)
 			} else {
    			out <- optim(par=alpha, fn=fn.acat, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, method="BFGS")
    			alpha<-out$par[1:(k-1)]
    			logLik <- -out$value
    			output<-list(alpha=alpha, y=y, w=w, x=x, probability.model=probability.model, link=link, logLik=logLik)
 			}      
		} else if (probability.model=="ForwardCR") {
        	if (dim(w)[2]!=0) {
    	       	out <- optim(par=c(alpha,zeta), fn=fn.fcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    			alpha <- out$par[1:(k-1)]
    			zeta <- out$par[k:length(out$par)]
    			names(zeta)=colnames(w)
    			logLik <- -out$value
    			output<-list(alpha=alpha, zeta=zeta, y=y, w=w, x=x, probability.model=probability.model, link=link, logLik=logLik)
 			} else {
    			out <- optim(par=alpha, fn=fn.fcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    			alpha<-out$par[1:(k-1)]
    			logLik <- -out$value
    			output<-list(alpha=alpha, y=y, w=w, x=x, probability.model=probability.model, link=link, logLik=logLik)
 			}      
		} else if (probability.model=="BackwardCR") {
        	if (dim(w)[2]!=0) {
    	       	out <- optim(par=c(alpha,zeta), fn=fn.bcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    			alpha <- out$par[1:(k-1)]
    			zeta <- out$par[k:length(out$par)]
    			names(zeta)=colnames(w)
    			logLik <- -out$value
    			output<-list(alpha=alpha, zeta=zeta, y=y, w=w, x=x, probability.model=probability.model, link=link, logLik=logLik)
 			} else {
    			out <- optim(par=alpha, fn=fn.bcr, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, Cum.Ymat=Cum.Ymat, link=link, method="BFGS")
    			alpha<-out$par[1:(k-1)]
    			logLik <- -out$value
    			output<-list(alpha=alpha, y=y, w=w, x=x, probability.model=probability.model, link=link, logLik=logLik)
 			}      
		} else if (probability.model=="Stereotype") {
        	if (dim(w)[2]!=0) {
    	       	out <- optim(par=c(alpha, phi, zeta), fn=fn.stereo, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, 
    	       		method="L-BFGS-B",upper=c(rep(Inf,k-1),rep(1,k-2)),lower=c(rep(-Inf,k-1),rep(0,k-2)))
    			alpha <- out$par[1:(k-1)]
    			phi <- c(1,out$par[k:(2*k-3)])
    			names(phi)<- paste("phi",1:(k-1),sep=".")
    			zeta <- out$par[(2*k-2):length(out$par)]
    			names(zeta)=colnames(w)
    			logLik <- -out$value
    			output<-list(alpha=alpha, phi=phi, zeta=zeta, y=y, w=w, x=x, scale=scale, probability.model=probability.model, link=link, logLik=logLik)
 			} else {
    			out <- optim(par=c(alpha, phi), fn=fn.stereo, w=w, x=x, beta=beta, y=y, k=k, Ymat=Ymat, 
    				method="L-BFGS-B",upper=c(rep(Inf,k-1),rep(1,k-2)),lower=c(rep(-Inf,k-1),rep(0,k-2)))
    			alpha<-out$par[1:(k-1)]
    			phi<-c(1,out$par[k:length(out$par)])
    			names(phi)<- paste("phi",1:(k-1),sep=".")
    			logLik <- -out$value
    			output<-list(alpha=alpha, phi=phi, y=y, w=w, x=x, scale=scale, probability.model=probability.model, link=link, logLik=logLik)
 			}      
		}
    }
class(output)<-"ordinalgmifs"
output 
}
