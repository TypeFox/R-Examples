"plot.imds"<-
function(x, xylim="auto",clab=1:nrow(X),lab.cex=1,lab.col="black",...){
		if(class(x)[2]=="sph"){
			X <- x$X
			r <- x$r
			start=0
			end=360
			length=100
			func=lines
			theta <- c(seq(start/180*pi, end/180*pi, length=length), end/180*pi)
			n <- nrow(X)
			if(xylim=="auto"){
				xmax <- max(X[,1] + r)
				xmin <- min(X[,1] - r)
				ymax <- max(X[,2] + r)
				ymin <- min(X[,2] - r)
				xrange <- xmax - xmin
				yrange <- ymax - ymin
				range <- max(c(xrange,yrange))
				xmar <- (range - xrange)/2
				ymar <- (range - yrange)/2
				xylim <- matrix(0,2,2)
				xylim[1,] <- c(xmin-xmar,xmax+xmar)
				xylim[2,] <- c(ymin-ymar,ymax+ymar)
			}
			plot(X,xlim=xylim[1,],ylim=xylim[2,],type="n",xlab="",ylab="",main="")
			text(X,labels=clab,cex=lab.cex,xlim=xylim[1,],ylim=xylim[2,])
			for(i in 1:n){
				func(r[i]*cos(theta)+X[i,1], r[i]*sin(theta)+X[i,2], ...)
			}
		}else if(class(x)[2]=="box"){
			X <- x$X
			R <- x$R
			n <- nrow(X)
			if(xylim=="auto"){
				xmax <- max(X[,1] + R[,1])
				xmin <- min(X[,1] - R[,1])
				ymax <- max(X[,2] + R[,2])
				ymin <- min(X[,2] - R[,2])
				xrange <- xmax - xmin
				yrange <- ymax - ymin
				range <- max(c(xrange,yrange))
				xmar <- (range - xrange)/2
				ymar <- (range - yrange)/2
				xylim <- matrix(0,2,2)
				xylim[1,] <- c(xmin-xmar,xmax+xmar)
				xylim[2,] <- c(ymin-ymar,ymax+ymar)
			}
			plot(X,xlim=xylim[1,],ylim=xylim[2,],type="n",xlab="",ylab="",main="")
			text(X,labels=clab,cex=lab.cex,xlim=xylim[1,],ylim=xylim[2,],col=lab.col)
			for(i in 1:n){
				xleft <- X[i,1] - R[i,1]
				ybottom <-  X[i,2] - R[i,2]
				xright <- X[i,1] + R[i,1]
				ytop <- X[i,2] + R[i,2]
				rect(xleft,ybottom,xright,ytop,...)
			}
		}
		}#end