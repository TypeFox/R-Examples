depth.KFSD=function (fdataobj, fdataori = fdataobj, trim = 0.25,
                     h=NULL,scale = FALSE, draw = FALSE){
    if (is.fdata(fdataobj)) {
        fdat <- TRUE
        if (is.null(rownames(fdataobj$data))) 
            rownames(fdataobj$data) <- 1:nrow(fdataobj$data)
        nms <- rownames(fdataobj$data)
        m0 <- nrow(fdataobj)
        fdataobj <- na.omit.fdata(fdataobj)
        fdataori <- na.omit.fdata(fdataori)
        nas <- na.action(fdataobj)
        nullans <- !is.null(nas)
        data <- fdataobj[["data"]]
        data2 <- fdataori[["data"]]
        names1 <- names2 <- names <- fdataobj[["names"]]
        names1$main <- "depth.KFSD median"
        names2$main <- paste("depth.KFSD trim ", trim * 100, 
            "%", sep = "")
        tt = fdataobj[["argvals"]]
        rtt <- fdataobj[["rangeval"]]
    }
    else {
        stop("no fdata class object")
        data <- fdataobj
        data2 <- fdataori
        fdat <- FALSE
    }
    n <- nrow(data)
    m <- ncol(data)
    m2 <- ncol(data2)
    n2 <- nrow(data2)
    if (is.null(n) && is.null(m)) 
        stop("ERROR IN THE DATA DIMENSIONS")
    if (is.null(m) && is.null(m2)) 
        stop("ERROR IN THE DATA DIMENSIONS")
    mdist=matrix(NA,ncol=n2,nrow=n2)
    for (i in 1:(n2-1)){for (j in (i+1):n2){
      mdist[i,j]<-mdist[j,i]<-norm.fdata(fdataori[i]-fdataori[j])
    }}
    if (is.null(h))   {
	  h<-0.15
	  hq2=quantile(mdist,probs=h,na.rm=TRUE)
	  #print("es nulo h")  
	}
	else {
	  #cat("no es nulo h ",h,"\n")    
	  if (is.numeric(h))    hq2<-h  
	  else hq2=quantile(mdist,probs=as.numeric(substring(h,first=3)),na.rm=TRUE)
	}	
	kern=function(x,y,h=hq2){exp(-norm.fdata(x-y)^2/h^2)}
	K0=rep(1,nrow(fdataobj))
	M1=matrix(NA,nrow=n2,ncol=n2)
	M2=matrix(NA,nrow=n,ncol=n2)
	M=array(NA,dim=c(n,n2,n2))
	for (i in 1:n2){for (j in i:n2){
	if (i==j) M1[i,i]=K0[i] else M1[i,j]<-M1[j,i]<-kern(fdataori[i],fdataori[j],h=hq2)
	}}
	same.dim <- FALSE
	if (n==n2 & m==m2){ same.dim <- TRUE}
	if (same.dim)
	  if (all(fdataobj$data==fdataori$data)) {
	    M2=M1
	  } else {sam.dim <- FALSE}
	if (!same.dim){	
  	for (i in 1:n){for (j in 1:n2){
  	if (all(fdataobj[i]$data == fdataori[j]$data)) M2[i,j]=K0[i] else M2[i,j]=kern(fdataobj[i],fdataori[j],h=hq2)
  	}}
	}  
	
#	print(M1)
#	print(M2)
	for (i in 1:n){for (j in 1:n2){for (k in 1:n2){
#	M[i,j,k]<-(K0-kern(fdataori[j],fdataori[k],h=h)-kern(fdataobj[i],fdataori[j],h=h)-kern(fdataobj[i],fdataori[k],h=h))/
#	(sqrt(2*K0-2*kern(fdataobj[i],fdataori[j],h=h))*sqrt(2*K0-2*kern(fdataobj[i],fdataori[k],h=h)))
	M[i,j,k]<-(K0[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(2*K0[i]-2*M2[i,j])*sqrt(2*K0[i]-2*M2[i,k]))
	}}}
	l=which(!is.finite(M))
	M[l]=NA
#	print(M)
#	print(apply(M,1,sum,na.rm=TRUE))
	ans=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
    if (scale) {
        mn <- min(ans, na.rm = TRUE)
        mx <- max(ans, na.rm = TRUE)
        ans = as.vector(ans/mx)
    }
    k = which.max(ans)
    med = data[k, ]
    lista = which(ans >= quantile(ans, probs = trim, na.rm = T))
    if (nullans) {
        ans1 <- rep(NA, len = m0)
        ans1[-nas] <- ans
        ans <- ans1
    }
    names(ans) <- nms
    if (length(lista) == 1) {
        mtrim <- data[lista, ]
        if (draw) {
            draw = FALSE
            warning("The plot is not shown")
        }
    }
    else mtrim = apply(fdataobj[lista]$data, 2, mean, na.rm = TRUE)
    tr <- paste("KFSD.tr", trim * 100, "%", sep = "")
    if (fdat) {
        med <- fdata(med, tt, rtt, names1)
        mtrim <- fdata(mtrim, tt, rtt, names2)
        rownames(med$data) <- "KFSD.med"
        rownames(mtrim$data) <- tr
        if (draw) {
            if (!scale) {
                mn <- min(ans, na.rm = TRUE)
                mx <- max(ans, na.rm = TRUE)
                scl <- mx - mn
            }
            ind1 <- !is.nan(ans)
            ans[is.nan(ans)] = NA
            cgray = 1 - (ans - mn)/(scl)
            plot(fdataori, col = gray(0.9), lty = 1, main = "KFSD Depth")
            lines(fdataobj[ind1, ], col = gray(cgray[ind1]))
            lines(mtrim, lwd = 2, col = "yellow")
            lines(med, col = "red", lwd = 2)
            legend("topleft", legend = c(tr, "Median"), lwd = 2, 
                col = c("yellow", "red"), box.col = 0)
        }
    }
    out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = lista, 
        dep = ans, h = h, hq=hq2)
    if (scale) 
        out$dscale = mx
    return(invisible(out))	
}

depth.FSD=function (fdataobj, fdataori = fdataobj, 
                    trim = 0.25, scale = FALSE, draw = FALSE){
    if (is.fdata(fdataobj)) {
        fdat <- TRUE
        if (is.null(rownames(fdataobj$data))) 
            rownames(fdataobj$data) <- 1:nrow(fdataobj$data)
        nms <- rownames(fdataobj$data)
        m0 <- nrow(fdataobj)
        fdataobj <- na.omit.fdata(fdataobj)
        fdataori <- na.omit.fdata(fdataori)
        nas <- na.action(fdataobj)
        nullans <- !is.null(nas)
        data <- fdataobj[["data"]]
        data2 <- fdataori[["data"]]
        names1 <- names2 <- names <- fdataobj[["names"]]
        names1$main <- "depth.FSD median"
        names2$main <- paste("depth.FSD trim ", trim * 100, 
            "%", sep = "")
        tt = fdataobj[["argvals"]]
        rtt <- fdataobj[["rangeval"]]
    }
    else {
        stop("no fdata class object")
        data <- fdataobj
        data2 <- fdataori
        fdat <- FALSE
    }
    n <- nrow(data)
    m <- ncol(data)
    m2 <- ncol(data2)
    n2 <- nrow(data2)
    if (is.null(n) && is.null(m)) 
        stop("ERROR IN THE DATA DIMENSIONS")
    if (is.null(m) && is.null(m2)) 
        stop("ERROR IN THE DATA DIMENSIONS")
	kern=function(x,y){inprod.fdata(x,y)}
	
	K01=norm.fdata(fdataobj)^2
	K02=norm.fdata(fdataori)^2
	M1=matrix(NA,nrow=n2,ncol=n2)
	M2=matrix(NA,nrow=n,ncol=n2)
	M=array(NA,dim=c(n,n2,n2))
	for (i in 1:n2){for (j in i:n2){
	if (i==j) M1[i,i]=K02[i] else M1[i,j]<-M1[j,i]<-kern(fdataori[i],fdataori[j])
	}}
  same.dim <- FALSE
  if (n==n2 & m==m2){ same.dim <- TRUE}
    if (same.dim)
      if (all(fdataobj$data==fdataori$data)) {
    	 M2=M1
	    } else {sam.dim <- FALSE}
  if (!same.dim){
  	for (i in 1:n){for (j in 1:n2){
  	if (all(fdataobj[i]$data == fdataori[j]$data)) M2[i,j]=K01[i] else M2[i,j]=kern(fdataobj[i],fdataori[j])
  }}
	}
	for (i in 1:n){for (j in 1:n2){for (k in 1:n2){
	if (all(fdataobj[i]$data == fdataori[j]$data) | all(fdataobj[i]$data == fdataori[k]$data)) M[i,j,k]=NA else M[i,j,k]<-(K01[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(K01[i]+K02[j]-2*M2[i,j])*sqrt(K01[i]+K02[k]-2*M2[i,k]))
	}}}
	l=which(!is.finite(M))
	M[l]=NA
	ans=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
    if (scale) {
        mn <- min(ans, na.rm = TRUE)
        mx <- max(ans, na.rm = TRUE)
        ans = as.vector(ans/mx)
    }
    k = which.max(ans)
    med = data[k, ]
    lista = which(ans >= quantile(ans, probs = trim, na.rm = T))
    if (nullans) {
        ans1 <- rep(NA, len = m0)
        ans1[-nas] <- ans
        ans <- ans1
    }
    names(ans) <- nms
    if (length(lista) == 1) {
        mtrim <- data[lista, ]
        if (draw) {
            draw = FALSE
            warning("The plot is not shown")
        }
    }
    else mtrim = apply(fdataobj[lista]$data, 2, mean, na.rm = TRUE)
    tr <- paste("FSD.tr", trim * 100, "%", sep = "")
    if (fdat) {
        med <- fdata(med, tt, rtt, names1)
        mtrim <- fdata(mtrim, tt, rtt, names2)
        rownames(med$data) <- "FSD.med"
        rownames(mtrim$data) <- tr
        if (draw) {
            if (!scale) {
                mn <- min(ans, na.rm = TRUE)
                mx <- max(ans, na.rm = TRUE)
                scl <- mx - mn
            }
            ind1 <- !is.nan(ans)
            ans[is.nan(ans)] = NA
            cgray = 1 - (ans - mn)/(scl)
            plot(fdataori, col = gray(0.9), lty = 1, main = "FSD Depth")
            lines(fdataobj[ind1, ], col = gray(cgray[ind1]))
            lines(mtrim, lwd = 2, col = "yellow")
            lines(med, col = "red", lwd = 2)
            legend("topleft", legend = c(tr, "Median"), lwd = 2, 
                col = c("yellow", "red"), box.col = 0)
        }
    }
    out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = lista, 
        dep = ans)
    if (scale) 
        out$dscale = mx
    return(invisible(out))
	
	}

