dplot.cavgb2 <- function(group,x,shape1, scale, shape2, shape3, pl0, pl, w=rep(1,length(x)),
		xmax = max(x)*(2/3), ymax=2e-05, decomp="r", choicecol=1:length(levels(group)),xlab=""){
  par(mfrow=c(2,1))
  K <- length(levels(group))
  L <- length(pl0)
  error <- "FALSE"
  for (k in 1:K){
    dPk <- length(unique(pl[group==levels(group)[k]]))  
    if (dPk>L){
      error <- "TRUE"
      warning("the estimated probabilities are not uniquely defined for group ", levels(group)[k])
      }     
  	}
  if (error) return()
  for (k in 1:K){
    pk <- as.vector( unique(pl[group==levels(group)[k],]))
    fk <- function(x) dcgb2(x,shape1, scale, shape2, shape3,pl0,pk,decomp=decomp)
    sub=paste("pl0 = (",round(pl0[1],3))
    pl1 <- length(pl0)-1
    if (pl1 >= 2){   
      for (i in 2:pl1) {
        sub <- paste(sub,",", round(pl0[i],3))
      }
    }
    sub <- paste(sub,",",round(pl0[pl1+1],3),")")
#    xmax <- max(x)*2/3                                # change 28.04.2014

    if (k==1){
      curve(fk,col=choicecol[k],from=0,to=xmax,lwd=2,ylab="Density",xlab=xlab,
        main="Compound densities per group",  ylim=c(0,ymax))
    }
    else {
      curve(fk,col=choicecol[k],lwd=2,lty=k,add=TRUE)
    }
  }

#  print("Please, place the cursor for the legend",quote = FALSE)           # change 2014-05-19
    legend("topright",levels(group), lwd=2,col=choicecol,lty=1:K)           # change 2014-05-19


# empirical counterparts
  for (k in 1:K){
    rdk <- x[group==levels(group)[k]]
		wk <- w[group==levels(group)[k]]
		wk <- wk/sum(wk)
		densk <- density(rdk,weights=wk,kernel="epanechnikov")
		if (k==1){
			plot(densk,col=choicecol[k],lwd=2,main="Kernel density estimate per group",
				xlab=xlab,xlim=c(0,xmax), ylim=c(0,ymax))
		}
		else { lines(densk,col=choicecol[k],lwd=2,lty=k)
	  
		}
	}
  par(mfrow=c(1,1))
}