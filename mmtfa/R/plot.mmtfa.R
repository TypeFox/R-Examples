plot.mmtfa <- function(x, xmarg=1, ymarg=2, res=200, levels=c(seq(.01,1,by=0.01), 0.001), what=c("contour","uncertainty"), main=NULL, xlab=NULL, legend=TRUE, ...){
	teigen <- x
	G <- ncol(teigen$bz)
	classcolours <- rainbow(length(unique(teigen$bestz)))
	numscr <- length(what)
	if(ncol(teigen$x)>1){
		if(numscr==2){
			par(mfrow=c(2,1))
		}
		if("contour"%in%what){
			plot(teigen$x[,c(xmarg,ymarg)], col=classcolours[teigen$bestz], pch=20, main="", ...)
      if(is.null(main)){title("Marginal Contour Plot")}else{title(main)}
			if(legend){
			  legendvec <- paste("Group",1:G)
			  legend("topleft", legendvec, col=c(classcolours), pch=rep(1,G))
			}
      lims <- par()$usr
			xseq <- seq(lims[1], lims[2], length.out=res)
			yseq <- seq(lims[3], lims[4], length.out=res)
			seqmat <- matrix(NA, res^2, 2)
			seqmat[,1] <- rep(xseq, each=res)
			seqmat[,2] <- rep(yseq, res)
			val <- matrix(NA,res,res)
			if(!teigen$info$univar){
				sigmainv <- array(NA, dim=c(2,2,teigen$bestg))
				for(g in 1:teigen$bestg){
					sigmainv[,,g] <- solve(teigen$par$sigma[c(xmarg,ymarg),c(xmarg,ymarg),g])
				}
			}
			delt <- deltaup(seqmat, matrix(teigen$par$mean[,c(xmarg,ymarg)],nrow=teigen$bestg,ncol=2), teigen$par$sigma[c(xmarg,ymarg),c(xmarg,ymarg),],sigmainv, teigen$bestg, res^2, teigen$info$univar)
			dens <- rowSums(exp(tft(seqmat,teigen$bestg,colSums(teigen$bz)/nrow(teigen$x),teigen$par$df,2,matrix(teigen$par$mean[,c(xmarg,ymarg)],nrow=teigen$bestg,ncol=2),sigmainv,res^2,array(teigen$par$sigma[c(xmarg,ymarg),c(xmarg,ymarg),], dim=c(2,2,teigen$bestg)),teigen$info$univar,delt,teigen$info$gauss)))
			val <- matrix(dens, res, res, byrow=TRUE)
			contour(x=xseq, y=yseq, z=val, add=TRUE, levels=levels, col=rgb(0.5,0.5,0.5,alpha=0.7))
		}
		if("uncertainty"%in%what){
			plot(teigen$x[,c(xmarg,ymarg)], col=classcolours[teigen$bestz],cex=(2*(1-apply(teigen$bz,1,max))), pch=20, ...)
			if(is.null(main)){title("Uncertainty Plot")}else{title(main)}
			if(legend){
			  legendvec <- paste("Group",1:G)
			  legend("topleft", legendvec, col=c(classcolours), pch=rep(1,G))
			}
			if(numscr==2){
  			par(mfrow=c(1,1))
			}
		}
	}
	else{
		pigs <- colSums(teigen$fuzzy)/nrow(teigen$x)
		dunivt <- function(xdum,df,sig,mean,pig,gauss){ 
      if(!gauss){
  			exp(log(pig)+lgamma((df+1)/2)-(1/2)*log(sig)-((1/2)*(log(pi)+log(df))+lgamma(df/2)+((df+1)/2)*(log(1+ mahalanobis(matrix(xdum,ncol=1), mean, 1/sig, inverted=TRUE)/df))))
      }
      else{
        dnorm(xdum, mean=mean, sd=sqrt(sig))
      }
		}
    mixuniv <- function(x){
      summat <- rep(0, length(x))
      for(g in 1:G){
         summat <- summat + dunivt(x, teigen$par$df[g], teigen$par$sigma[,,g], teigen$par$mean[g,], pigs[g], teigen$info$gauss)
      }
      summat
    }
		bigdens <- NA
		for(g in 1:G){
			bigdens[g] <- dunivt(teigen$par$mean[g,],teigen$par$df[g],teigen$par$sigma[,,g],teigen$par$mean[g,],pigs[g], teigen$info$gauss)
		}
		plot(density(teigen$x), ylim=c(0,max(bigdens)+0.04*max(bigdens)), lty=4, main="", xlab="", ...)
		if(is.null(main)){title("Univariate Density Plot")}else{title(main)}
		if(is.null(xlab)){title(xlab=colnames(teigen$x)[xmarg])}else{title(xlab=xlab)}
		#univt <- function(xdum){ log(pig[g])+lgamma((teigen$par$df[g]+1)/2)-(1/2)*log(teigen$par$sigma[,,g])-((p/2)*(log(pi)+log(teigen$par$df[g]))+lgamma(teigen$par$df[g]/2)+((teigen$par$df[g]+p)/2)*(log(1+ mahalanobis(matrix(xdum,nrow(teigen$fuzzy),1), teigen$par$mug[g,], 1/teigen$par$sigma[,,g], inverted=TRUE)/teigen$par$df[g])))}
		curve(mixuniv(x), add=TRUE, col="black", n=res, lwd=2)
    for(g in 1:G){
			curve(dunivt(x,teigen$par$df[g],teigen$par$sigma[,,g],teigen$par$mean[g,],pigs[g], teigen$info$gauss),add=TRUE, col=classcolours[g], n=res)
		}
    if(legend){
      legendvec <- c("density()", "mixture", paste("Group",1:G)) 
      legend("topleft", legendvec, col=c("black", "black", classcolours), lty=c(4,1,rep(1,G)), lwd=c(1,2,rep(1,G)))
    } 
	}
}
