plotsingle <-
function(fit, d = "best", pl = -Inf, pu = Inf, ql = NA, qu = NA, sf = 3, ex = 1){
	
	if(d == "best"){
		ssq <- fit$ssq[ex, is.na(fit$ssq[ex,])==F]
		best.index <- which(ssq == min(ssq))[1]
	}
	index<-switch(which(d==c("normal", "t", "gamma", "lognormal", "logt","beta", "hist", "best")), 1, 2, 3, 4, 5, 6, 7, best.index)
	
	par(ps=15)
	par(mar = c(5.1, 5.1, 4.1, 2.1))
	
	if(index==1){
		
		if(pl == -Inf){pl <- qnorm(0.001, fit$Normal[ex,1], fit$Normal[ex,2])}
		if(pu == Inf){pu <- qnorm(0.999, fit$Normal[ex,1], fit$Normal[ex,2])}
		x <- seq(from = pl, to = pu, length = 200)
		fx <- dnorm(x, fit$Normal[ex,1], fit$Normal[ex,2]) 
		plot(x, fx , type = "l", lwd = 2, xlab = "x", ylab = expression(f[X](x)), main=paste("Normal (mean = ",signif(fit$Normal[ex,1], sf), ", sd = ", signif(fit$Normal[ex,2], sf), ")", sep=""))
	
		if(is.na(ql) == F){
			x.q1 <- qnorm(ql, fit$Normal[ex,1], fit$Normal[ex,2])
			if(x.q1 > pl){
			  x1<-seq(from = pl, to = x.q1 , length = ceiling((x.q1-pl)/(pu - pl)*200))
			  lines(x1, dnorm(x1, fit$Normal[ex,1], fit$Normal[ex,2] ), type="h", col = "red" )
			}
		}
	
		if(is.na(qu) == F){
			x.q2 <- qnorm(qu, fit$Normal[ex,1], fit$Normal[ex,2])
			if(x.q2 < pu){
			  x2<-seq(from = x.q2, to = pu, length = ceiling((pu - x.q2)/(pu - pl)*200))
			  lines(x2, dnorm(x2, fit$Normal[ex,1], fit$Normal[ex,2]) , type="h", col = "red" )
			}
		}
	}
	
	if(index==2){
		
		if(pl == -Inf){pl <- fit$Student.t[ex,1] + fit$Student.t[ex,2] * qt(0.001, fit$Student.t[ex,3])}
		if(pu == Inf){pu <- fit$Student.t[ex,1] + fit$Student.t[ex,2] * qt(0.999, fit$Student.t[ex,3])}
		
		x <- seq(from = pl, to = pu, length = 200)
		fx <- dt((x - fit$Student.t[ex,1])/fit$Student.t[ex,2], fit$Student.t[ex,3])/fit$Student.t[ex,2]
		plot(x, fx, type = "l", lwd = 2, xlab = "x", ylab = expression(f[X](x)), main=paste("Student-t(",signif(fit$Student.t[ex,1], sf), ", ", signif(fit$Student.t[ex,2], sf), ")", sep=""))
	
		if(is.na(ql) == F){
			x.q1 <- fit$Student.t[ex,1] + fit$Student.t[ex,2] * qt(ql, fit$Student.t[ex,3])
			if(x.q1 > pl){
			  x1<-seq(from = pl, to = x.q1 , length = ceiling((x.q1-pl)/(pu - pl)*200))
			  lines(x1, dt((x1 - fit$Student.t[ex,1])/fit$Student.t[ex,2], fit$Student.t[ex,3])/fit$Student.t[ex,2], type="h", col = "red" )
			}
		}
	
		if(is.na(qu) == F){
			x.q2 <- fit$Student.t[ex,1] + fit$Student.t[ex,2] * qt(qu, fit$Student.t[ex,3])
			
			if(x.q2 < pu){
			  x2<-seq(from = x.q2, to = pu, length = ceiling((pu - x.q2)/(pu - pl)*200))
			  lines(x2, dt((x2 - fit$Student.t[ex,1])/fit$Student.t[ex,2], fit$Student.t[ex,3])/fit$Student.t[ex,2] , type="h", col = "red" )
			}
		}
	}
	
	if(index==3){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		
		if(pl == -Inf){pl <- xl + qgamma(0.001, fit$Gamma[ex,1], fit$Gamma[ex,2])}
		if(pu == Inf){pu <- xl + qgamma(0.999, fit$Gamma[ex,1], fit$Gamma[ex,2])}
		x <- seq(from = pl, to = pu, length = 200)
		fx <- dgamma(x - xl, fit$Gamma[ex,1], fit$Gamma[ex,2])  
		
		plot(x, fx, type = "l", lwd = 2, xlab = "x", ylab = expression(f[X](x)), main=paste("Gamma(",signif(fit$Gamma[ex,1], sf), ", ", signif(fit$Gamma[ex,2], sf), ")", sep=""))

	
		if(is.na(ql) == F){
			x.q1 <- xl + qgamma(ql, fit$Gamma[ex,1], fit$Gamma[ex,2])
			if(x.q1 > pl){
			  x1<-seq(from = pl, to = x.q1 , length = ceiling((x.q1-pl)/(pu - pl)*200))
			  lines(x1, dgamma(x1 - xl, fit$Gamma[ex,1], fit$Gamma[ex,2] ), type="h", col = "red" )
			}
		}
	
		if(is.na(qu) == F){
			x.q2 <- xl + qgamma(qu, fit$Gamma[ex,1], fit$Gamma[ex,2])
			if(x.q2 < pu){
			  x2<-seq(from = x.q2, to = pu, length = ceiling((pu - x.q2)/(pu - pl)*200))
			  lines(x2, dgamma(x2 - xl, fit$Gamma[ex,1], fit$Gamma[ex,2]) , type="h", col = "red" )
			}
		}
	}
	
	if(index==4){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		
		if(pl == -Inf){pl <- xl + qlnorm(0.001, fit$Log.normal[ex,1], fit$Log.normal[ex,2])}
		if(pu == Inf){pu <- xl + qlnorm(0.999, fit$Log.normal[ex,1], fit$Log.normal[ex,2])}
		x <- seq(from = pl, to = pu, length = 200)
		fx <- dlnorm(x - xl, fit$Log.normal[ex,1], fit$Log.normal[ex,2]) 
		plot(x, fx, type = "l", lwd = 2, xlab = "x", ylab = expression(f[X](x)), main=paste("Log normal(",signif(fit$Log.normal[ex,1], sf), ", ", signif(fit$Log.normal[ex,2], sf), ")", sep=""))

	
		if(is.na(ql) == F){
			x.q1 <- xl + qlnorm(ql, fit$Log.normal[ex,1], fit$Log.normal[ex,2])
			if(x.q1 > pl){
			  x1<-seq(from = pl, to = x.q1 , length = ceiling((x.q1-pl)/(pu - pl)*200))
			  lines(x1, dlnorm(x1 - xl, fit$Log.normal[ex,1], fit$Log.normal[ex,2] ), type="h", col = "red" )
			}
		}
	
		if(is.na(qu) == F){
			x.q2 <- xl + qlnorm(qu, fit$Log.normal[ex,1], fit$Log.normal[ex,2])
			if(x.q2 < pu){
			  x2<-seq(from = x.q2, to = pu, length = ceiling((pu - x.q2)/(pu - pl)*200))
			  lines(x2, dlnorm(x2 - xl, fit$Log.normal[ex,1], fit$Log.normal[ex,2]) , type="h", col = "red" )
			}
		}
	}	
	
	if(index==5){ # log student t
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		
    # Calculate axes limits using the lognormal; log-t limits may be too extreme
		if(pl == -Inf){pl <- xl + qlnorm(0.001, fit$Log.normal[ex,1], fit$Log.normal[ex,2])}
		if(pu == Inf){pu <- xl + qlnorm(0.999, fit$Log.normal[ex,1], fit$Log.normal[ex,2])}
    
	#	if(pl == -Inf){pl <- xl + exp(fit$Log.Student.t[ex,1] + fit$Log.Student.t[ex,2] * qt(0.001, fit$Log.Student.t[ex,3]))}
	#	if(pu == Inf){pu <- xl + exp(fit$Log.Student.t[ex,1] + fit$Log.Student.t[ex,2] * qt(0.999, fit$Log.Student.t[ex,3]))}
		x <- seq(from = pl, to = pu, length = 200)
		fx <- dt( (log(x - xl) - fit$Log.Student.t[ex,1]) / fit$Log.Student.t[ex,2], fit$Log.Student.t[ex,3]) / ((x - xl) * fit$Log.Student.t[ex,2])
		plot(x, fx, type = "l", lwd = 2, xlab = "x", ylab = expression(f[X](x)), main=paste("Log T(",signif(fit$Log.Student.t[ex,1], sf), ", ", signif(fit$Log.Student.t[ex,2], sf), ")", sep=""))

	
		if(is.na(ql) == F){
			x.q1 <- xl + exp(fit$Log.Student.t[ex,1] + fit$Log.Student.t[ex,2] * qt(ql, fit$Log.Student.t[ex,3]))
			if(x.q1 > pl){
			  x1<-seq(from = pl, to = x.q1 , length = ceiling((x.q1-pl)/(pu - pl)*200))
			  lines(x1, dt( (log(x1 - xl) - fit$Log.Student.t[ex,1])/fit$Log.Student.t[ex,2], fit$Log.Student.t[ex,3]) / ((x1 - xl) * fit$Log.Student.t[ex,2]), type="h", col = "red" )
			}
		}
	
		if(is.na(qu) == F){
			x.q2 <- xl + exp(fit$Log.Student.t[ex,1] + fit$Log.Student.t[ex,2] * qt(qu, fit$Log.Student.t[ex,3]))
			if(x.q2 < pu){
			  x2<-seq(from = x.q2, to = pu, length = ceiling((pu - x.q2)/(pu - pl)*200))
			  lines(x2, dt( (log(x2 - xl) - fit$Log.Student.t[ex,1])/fit$Log.Student.t[ex,2], fit$Log.Student.t[ex,3]) / ((x2 - xl) * fit$Log.Student.t[ex,2]) , type="h", col = "red" )
			}
		}
	}	

	
	
	if(index==6){
		xl <- fit$limits[ex,1]
		xu <- fit$limits[ex,2]
		if(xl == -Inf){xl <- 0}
		if(xu == Inf){xu <- 1}
		
		if(pl == -Inf){pl <- xl + (xu - xl) * qbeta(0.001, fit$Beta[ex,1], fit$Beta[ex,2])}
		if(pu == Inf){pu <- xl + (xu - xl) * qbeta(0.999, fit$Beta[ex,1], fit$Beta[ex,2])}
		x <-  seq(from = pl, to = pu, length = 200)
		fx <-  1/(xu - xl) * dbeta( (x - xl) / (xu - xl), fit$Beta[ex,1], fit$Beta[ex,2])
		plot(x, fx, type = "l", lwd = 2, xlab = "x", ylab = expression(f[X](x)), main=paste("Beta(",signif(fit$Beta[ex,1], sf), ", ", signif(fit$Beta[ex,2], sf), ")", sep=""))
	
		if(is.na(ql) == F){
			x.q1 <- xl + (xu - xl) * qbeta(ql, fit$Beta[ex,1], fit$Beta[ex,2])
			if(x.q1 > pl){
			  x1<-seq(from = pl, to = x.q1 , length = ceiling((x.q1-pl)/(pu - pl)*200))
			  lines(x1, 1/(xu - xl) * dbeta( (x1 - xl) / (xu - xl), fit$Beta[ex,1], fit$Beta[ex,2] ), type="h", col = "red" )
			}
		}
	
		if(is.na(qu) == F){
			x.q2 <- xl + (xu - xl) * qbeta(qu, fit$Beta[ex,1], fit$Beta[ex,2])
			if(x.q2 < pu){
			  x2<-seq(from = x.q2, to = pu, length = ceiling((pu - x.q2)/(pu - pl)*200))
			  lines(x2, 1/(xu - xl) * dbeta( (x2 - xl) / (xu - xl), fit$Beta[ex,1], fit$Beta[ex,2]) , type="h", col = "red" )
			}
		}
	}
	
	if(index==7){
	  
    
	  if(pl == -Inf & fit$limits[ex,1] > -Inf){pl <- fit$limits[ex,1]}
	  if(pu == Inf & fit$limits[ex,2] < Inf){pu <- fit$limits[ex,2] }
	  if(pl == -Inf & fit$limits[ex,1] == -Inf){pl <- qnorm(0.001, fit$Normal[ex,1], fit$Normal[ex,2])}
	  if(pu == Inf & fit$limits[ex,2] == Inf){pu <- qnorm(0.999, fit$Normal[ex,1], fit$Normal[ex,2])}
    
    p <- c(0, fit$probs[ex,], 1)
    x <- c(pl, fit$vals[ex,], pu)
    
    h <- rep(0, length(x) -1)
    for(i in 1:length(h)){
      h[i]<-(p[i+1] - p[i]) / (x[i+1]-x[i])
    }
    fx <- c(0, rep(h, each = 2), 0)
    x2 <- rep(x, each = 2)
	  
	  
    if(min(fx)<0){
      plot(x2, fx , type = "n", lwd = 2, xlab = "x", ylab = expression(f[X](x)), main="histogram fit")}else{
      plot(x2, fx , type = "l", lwd = 2, xlab = "x", ylab = expression(f[X](x)), main="histogram fit")
	    if(is.na(ql) == F){
	      x.q1 <- qhist(ql, x, p)
	      if(x.q1 > pl){
	        x1<-seq(from = pl, to = x.q1 , length = ceiling((x.q1-pl)/(pu - pl)*200))
	        lines(x1, dhist(x1, x, p), type="h", col = "red" )
	      }
	    }
	  
	    if(is.na(qu) == F){
	      x.q2 <- qhist(qu, x, p)
	      if(x.q2 < pu){
	        x2<-seq(from = x.q2, to = pu, length = ceiling((pu - x.q2)/(pu - pl)*200))
	        lines(x2, dhist(x2, x, p) , type="h", col = "red" )
	      }
	    }
    }
    
	}
	
}
