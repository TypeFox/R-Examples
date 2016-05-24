regimes.plot.MSAR <-
function(res,data,ex=1,col.l="red",nc=1,ylim=NULL,xlab="time",ylab="series",d=NULL,dt=1,lwd=1){    
	if (!inherits(res, "MSAR")) 
        stop("use only with \"MSAR\" objects")
    data = as.array(data)
    T = dim(data)[1]
    N.samples = dim(data)[2]
    if(is.null(N.samples)|is.na(N.samples)){N.samples <- 1}
    order = attributes(res$theta)$order
    M = attributes(res$theta)$NbRegimes
    if (N.samples==1){res$smoothedprob = array(res$smoothedprob,c(1,T-order-1,M))}
    d <- dim(data)[3]
	if(is.null(d)|is.na(d)){d <- 1}
	data = array(data,c(T,N.samples,d))
	if (is.null(ylim)) {ylim = c(min(data[,ex,nc]),max(data[,ex,nc]))}
	M = attributes(res$theta)$NbRegimes
	col = rev(gray(1:max(M,5)/max(M,5)))
	Reg = apply(res$smoothedprob[ex,,],1,which.max) 
	Reg = c(rep(Reg[1],length(data[,ex,nc])-length(Reg)),Reg) # A revoir
	plot((0:(T-1))*dt,data[,ex,nc],type="n",xlim=c(0,length(data[,ex,nc])*dt),ylim=ylim,xlab=xlab,ylab=ylab)
	if (length(unique(Reg))>1) {
		Reg = c(M+1,Reg,M+1)
		R = matrix(0,M,length(Reg))
		
		for (m in 1:M) { # loop could start at m=2...
			R[m,] = Reg==m
			if (sum(R[m,])>0) {
				xp = which(diff(as.numeric(R[m,]))>0)*dt
				xm = which(diff(as.numeric(R[m,]))<0)*dt
				if (xm[1]<xp[1]) {xm = xm[2:length(xm)]}
				for (k in 1:length(xp)) {
				polygon((c(xp[k], xm[k], xm[k], xp[k], NA)+c(xp[k]+dt, xm[k]+dt, xm[k]+dt, xp[k]+dt, NA))/2-dt, c(ylim[1], ylim[1], ylim[2], ylim[2], NA), col = col[m], border = NA, lwd = 3)
			}
}		}
		Reg = c(rep(NA,order),Reg[2:(length(Reg)-1)])
	} else {
		polygon(c(1,T,T,1,NA), c(ylim[1], ylim[1], ylim[2], ylim[2], NA), col = col[Reg[2]], border = NA, lwd = 3)
	}	
	d <- attributes(res$theta)$NbComp
	if (d==1) {lines((1:(T))*dt,data[,ex,],typ="l",col=col.l,lwd=lwd)}
	else {lines((1:(T))*dt,data[,ex,nc],typ="l",col=col.l,lwd=lwd)}
	return(Reg)
}
