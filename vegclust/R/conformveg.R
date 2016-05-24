conformveg<-function(x,y, fillvalue=0, verbose=FALSE) {
	nx = names(x)
	ny = names(y)
	if(is.null(nx) || is.null(ny)) {
		stop("x and y must have both column names")
	}
	nn = ny[!(ny %in% nx)]
	nm = c(nx,nn)
	if(verbose) {
		cat(paste(sum(ny %in% nx)," names in Y that are in X."))
		cat(paste(sum(!(ny %in% nx))," names in Y that are NOT in X: "))
		nynx = ny[!(ny %in% nx)]
		for(i in 1:length(nynx)) cat(paste(nynx[i],"-"))
		cat("\n")
		cat(paste(sum(nx %in% ny)," names in X that are in Y"))
		cat(paste(sum(!(nx %in% ny))," names in X that are NOT in Y: "))
		nxny = nx[!(nx %in% ny)]
		for(i in 1:length(nxny)) cat(paste(nxny[i],"-"))
		cat("\n")
	}
	
	xinf = data.frame(matrix(fillvalue,nrow=nrow(x),ncol=length(nm)))
	xinf[,1:length(nx)] = x	
	row.names(xinf) = row.names(x)
	names(xinf) = nm
	
	yinf = data.frame(matrix(fillvalue,nrow=nrow(y),ncol=length(nm)))
	#Shared columns
	sel = which(ny %in% nx)
	for(i in 1:length(sel)) yinf[,which(nx %in% ny[sel[i]])] = y[,sel[i]]
	
	if(length(nn)>0) {
		yinf[,((length(nx)+1):length(nm))] = y[,!(ny%in%nx)] #Remaining columns in y
	}
	row.names(yinf) = row.names(y)
	names(yinf) = nm
	
	return (list(x=xinf, y =yinf))
}