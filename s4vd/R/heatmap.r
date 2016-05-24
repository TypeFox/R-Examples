BCheatmap <- function(
		X,res,
		cexR=1.5,
		cexC=1.25,
		axisR=FALSE,
		axisC=TRUE,
		heatcols = maPalette(low="blue",mid="white",high="red", k=50),
		clustercols= c(1:5),
		allrows=FALSE,
		allcolumns=TRUE
)
{
	number <- res@Number
	layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE),widths=c(8.5,1.5),heights=c(1,1))
	if(number==1){
		rowmat <- res@RowxNumber
		colmat <- t(res@NumberxCol)
		roworder <- c(which(res@RowxNumber[,1]),which(!res@RowxNumber[,1]))
		colorder <- c(which(res@NumberxCol[1,]),which(!res@NumberxCol[1,]))
		X <- X[roworder,colorder]
		par(mar=c(5, 1, 1, 3))
		image(t(X)[,nrow(X):1],col=heatcols,x=c(1:ncol(X)),y=c(1:nrow(X)),axes=F,ylab="",xlab="")
		if(axisC)axis(1,1:ncol(X),labels=colnames(X),las=2,line = -0.5, tick = 0,cex.axis = cexC)
		if(axisR)axis(4,1:nrow(X),labels=rownames(X),las=2,line = -0.5, tick = 0,cex.axis = cexR)
		rin1 <- which(roworder %in% which(rowmat[,1]))
		cin1 <- which(colorder %in% which(colmat[,1]))
		nr <- length(roworder)
		nc <- length(colorder)
		if(allrows) nr <- nrow(X)
		if(allcolumns) nc <- ncol(X)
		xl <- 0.5
		yb <- nr-length(rin1)+0.5
		xr <- length(cin1)+0.5
		yt <- nr+0.5
		rect(xleft=xl,ybottom=yb,xright=xr,ytop=yt,density=0,angle=25,lwd=4,col=clustercols[1])
	}else{
	rowmat <- res@RowxNumber
	overlap <- rowSums(rowmat)
	roworder <- which(overlap==number) # in all
	if(number>2){
		for(i in 1:(number-2)){
			innext <- intersect(which(rowmat[,i]),which(rowmat[,i+1]))
			nooverlap <- which(rowmat[,i]&rowSums(rowmat)==1)
			for(l in 1:(number-1-i)){ 
				temp   <- intersect(which(rowmat[,i]),which(rowmat[,i+1+l])) 
				temp   <- temp[!intersect(which(rowmat[,i]),which(rowmat[,i+1+l])) %in% innext]
				roworder <- unique(c(roworder,temp))
			}
			roworder <- unique(c(roworder,nooverlap))
			roworder <- unique(c(roworder,innext))
		}
	}
	innext <- intersect(which(rowmat[,number-1]),which(rowmat[,number]))
	nooverlap <- which(rowmat[,number-1]&rowSums(rowmat)==1)
	roworder <- unique(c(roworder,nooverlap))
	roworder <- unique(c(roworder,innext))
	nooverlap <- which(rowmat[,number]&rowSums(rowmat)==1)
	roworder <- unique(c(roworder,nooverlap))
	if(allrows) roworder <- c(roworder,which(!1:nrow(rowmat)%in%roworder)) 
	colmat <- t(res@NumberxCol)
	overlap <- rowSums(colmat)
	colorder <- which(overlap==number) # in all
	if(number>2){
		for(i in 1:(number-2)){
			innext <- intersect(which(colmat[,i]),which(colmat[,i+1]))
			nooverlap <- which(colmat[,i]&rowSums(colmat)==1)
			for(l in 1:(number-1-i)){ 
				temp   <- intersect(which(colmat[,i]),which(colmat[,i+1+l])) 
				temp   <- temp[!intersect(which(colmat[,i]),which(colmat[,i+1+l])) %in% innext]
				colorder <- unique(c(colorder,temp))
			}
			colorder <- unique(c(colorder,nooverlap))
			colorder <- unique(c(colorder,innext))
		}
	}	
	innext <- intersect(which(colmat[,number-1]),which(colmat[,number]))
	nooverlap <- which(colmat[,number-1]&rowSums(colmat)==1)
	colorder <- unique(c(colorder,nooverlap))
	colorder <- unique(c(colorder,innext))
	nooverlap <- which(colmat[,number]&rowSums(colmat)==1)
	colorder <- unique(c(colorder,nooverlap))
	if(allcolumns) colorder <- c(colorder,which(!1:nrow(colmat)%in%colorder)) 
	X <- X[roworder,colorder]
	par(mar=c(5, 1, 1, 3))# c(bottom, left, top, right)
	image(t(X)[,nrow(X):1],col=heatcols,x=c(1:ncol(X)),y=c(1:nrow(X)),axes=F,ylab="",xlab="")
	if(axisC)axis(1,1:ncol(X),labels=colnames(X),las=2,line = -0.5, tick = 0,cex.axis = cexC)
	if(axisR)axis(4,1:nrow(X),labels=rownames(X),las=2,line = -0.5, tick = 0,cex.axis = cexR)
	rin1 <- which(roworder %in% which(rowmat[,1]))
	cin1 <- which(colorder %in% which(colmat[,1]))
	nr <- length(roworder)
	nc <- length(colorder)
	if(allrows) nr <- nrow(X)
	if(allcolumns) nc <- ncol(X)
	xl <- 0.5
	yb <- nr-length(rin1)+0.5
	xr <- length(cin1)+0.5
	yt <- nr+0.5
	rect(xleft=xl,ybottom=yb,xright=xr,ytop=yt,density=0,angle=25,lwd=4,col=clustercols[1])
	for(i in 2:number){
			rin <- which(roworder %in% which(rowmat[,i]))
			rstart <- numeric()
			rstop <- numeric()
			e <- 1
			rstart[e] <- rin[1]
			for(j in 2:length(rin)){
				if(rin[j-1]-rin[j]!=-1){
					rstop[e] <- rin[j-1]
					e <- e+1
					rstart[e] <- rin[j]
				}
			}
			rstop[e] <- rin[j]
			
			cin <- which(colorder %in% which(colmat[,i]))
			cstart <- numeric()
			cstop <- numeric()
			e <- 1
			cstart[e] <- cin[1]
			for(j in 2:length(cin)){
				if(cin[j-1]-cin[j]!=-1){
					cstop[e] <- cin[j-1]
					e <- e+1
					cstart[e] <- cin[j]
				}
			}
			cstop[e] <- cin[j]
			
			for(j in 1:length(rstart)){
				for(k in 1:length(cstart)){
					xl <- cstart[k] - 0.5
					yb <- nr - rstop[j] + .5
					xr <- cstop[k]+0.5
					yt <- nr - rstart[j] + 1.5
					rect(xleft=xl,ybottom=yb,xright=xr,ytop=yt,density=0,angle=45*i,lwd=4,col=clustercols[i])
				}
			} 
		}
	}
	min.raw <- min(X)
	max.raw <- max(X)
	z <- seq(min.raw, max.raw, length=length(heatcols))
	image(z=t(matrix(z, ncol=1)),col=heatcols, 
			xaxt="n", yaxt="n")
	axis(4,at=seq(0,1,by=.5),cex.axis=1.5,labels=c(round(min.raw,digits=1),0,round(max.raw,digits=1)))
}

maPalette <- 
function (low = "white", high = c("green", "red"), mid = NULL, 
          k = 50) 
{
  low <- col2rgb(low)/255
  high <- col2rgb(high)/255
  if (is.null(mid)) {
    r <- seq(low[1], high[1], len = k)
    g <- seq(low[2], high[2], len = k)
    b <- seq(low[3], high[3], len = k)
  }
  if (!is.null(mid)) {
    k2 <- round(k/2)
    mid <- col2rgb(mid)/255
    r <- c(seq(low[1], mid[1], len = k2), seq(mid[1], high[1], 
                                              len = k2))
    g <- c(seq(low[2], mid[2], len = k2), seq(mid[2], high[2], 
                                              len = k2))
    b <- c(seq(low[3], mid[3], len = k2), seq(mid[3], high[3], 
                                              len = k2))
  }
  rgb(r, g, b)
}