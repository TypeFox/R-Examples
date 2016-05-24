weig <- function(x,data=NULL, optcomp=NULL, yweights=FALSE){
	## 14/10/09 IH
	
	## Internal function for prmsDACV and sprmsDACV
	
	if (is.null(optcomp)){
		optcomp <- x$inputs$a
	}
	n <- dim(data)[1]
	ynames <- names(x$YMeans)
	# used.vars <- x$used.vars[[2*optcomp]] # bei prm?
	
	if(missing(data)){ 
		yp <- predict(x)
		y0 <- x$inputs$y0
		Xn <- as.matrix(x$inputs$Xs[,])
		n <- dim(Xn)[1]
	} else { 
		yp <- predict(x,data)
		yindex <- which(colnames(data)==ynames)
		if (sum(unique(data[,yindex])%in%c(-1,1))<2){
			yfac <- as.factor(data[,1])
			levels(yfac) <- c("-1", "1")
			data[,yindex] <- as.numeric(as.character(yfac))
		}
		y0 <- data[,yindex]		
		
		Xnames <- names(x$XMeans) #
		Xindex <- which(colnames(data)%in%Xnames) #
		if (length(Xindex)!=length(Xnames)){
			stop("Column names of data don't match variable names in the model.")
		}
		Xn <- scale(data[,Xindex], center=x$XMeans, scale=x$Xscales)
		if (length(rownames(data))==0){
			rownames(Xn) <- 1:dim(Xn)[1]
		} else{
			rownames(Xn) <- rownames(data) 
		}
	}
	# if(is.null(used.vars)){used.vars=1:dim(Xn)[2]}
	
	ind1 <- which(y0==1)
	mT1 <- apply(as.matrix(x$scores[which(x$input$y0==1),]),2,median)
	sT1 <- apply(as.matrix(x$scores[which(x$input$y0==1),]),2,qn)
	Tn1 <- Xn[ind1,]%*%x$R

	ind2 <- which(y0==-1)  
	mT2 <- apply(as.matrix(x$scores[which(x$input$y0==-1),]),2,median)
	sT2 <- apply(as.matrix(x$scores[which(x$input$y0==-1),]),2,qn)
	Tn2 <- Xn[ind2,]%*%x$R
	
	cov1 <- covMcd(x$scores)
	cov2 <- covMcd(x$scores)
	
	wlist <- int_weight(Tn1, Tn2, ind1, ind2, y0, fun=x$inputs$fun, probp1=x$inputs$constants[1], hampelp2=x$inputs$constants[2], hampelp3=x$inputs$constants[3], probp4=yweights, yweights=yweights, center1=cov1$center, cov1=cov1$cov, center2=cov2$center, cov2=cov2$cov)
	w <- wlist$we
	return(w)
	
}