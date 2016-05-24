setMethod("lda",
	signature(x = "IdtClMANOVA"),
	function(x,prior="proportions",selmodel=BestModel(H1res(x)),egvtol=1.0e-9,...)
	{
        	p <- 2*x@NIVar
     		nk <- as.numeric(table(x@grouping))
		N <- sum(nk)
     		k <- length(nk) 
		if (k==1) stop("The data belongs to one single group. A partition into at least two different groups is required\n")
		if (is.character(selmodel))  selmodel <- sapply(selmodel,function(mod) which(mod==x@H0res@ModelNames))
		means <- x@H1res@mleNmuE
		if (prior[1]=="proportions") prior <- nk/N
		names(prior) <- rownames(means)
		W <- x@H1res@Configurations[[selmodel]]$mleSigE	
		Wi <- solve(W)					
		B <- x@H0res@Configurations[[selmodel]]$mleSigE - W	
		WiBdecp <- eigen(Wi%*%B)	
		if (selmodel != 5) {
			if (selmodel==1)  r <- min(p,k-1)
			else  {
				WiBegval <- Re(WiBdecp$values)	
				posWiBegval <- WiBegval[WiBegval>egvtol]
				r <- length(posWiBegval)
			}
			eigvct <- Re(WiBdecp$vectors[,1:r])
			if (r==1) dim(eigvct) <- c(p,1)
			sclvar <- apply(eigvct,2,function(v) v%*%W%*%v)	
			scaling <- scale(eigvct,center=FALSE,scale=sqrt(sclvar))
			dimnames(scaling) <- list(rownames(W),paste("LD",1:r,sep=""))	
			attr(scaling,"scaled:scale") <- NULL
        	}
		else {
			Wd <- diag(x@H1res@Configurations[[5]]$mleSigE)
			scaling <- diag(1/sqrt(Wd))
			dimnames(scaling) <- list(rownames(W),paste("LD",1:p,sep=""))
        	}
		new("Idtlda",prior=prior,means=means,scaling=scaling,N=N) 
	}
)

setMethod("lda",
	signature(x = "IData"),
	function(x, grouping, prior="proportions", tol=1.0e-4, subset=1:nrow(x),Config=1:5,SelCrit=c("AIC","BIC"))
	{
	   	SelCrit <- match.arg(SelCrit)
		if (length(subset) < nrow(x)) {
			x <- x[subset,]
			grouping <- grouping[subset]
		}
   		levels(grouping)[table(grouping)==0] <- NA
		if (length(levels(grouping))==1) 
			stop("The data belongs to one single group. A partition into at least two different groups is required\n")
		lda(MANOVA(x,grouping,Mxt="Hom",Model="Normal",Config=Config,SelCrit=SelCrit,tol=tol),prior=prior)
	}
)

setMethod("predict",
	signature(object = "Idtlda"),
	function(object,newdata,prior=object@prior,...)
	{
   		if (is(newdata,"IData")) newdata <- as.matrix(cbind(newdata@MidP,newdata@LogR))
   		if (is(newdata,"data.frame")) newdata <- as.matrix(newdata)
   		n <- nrow(newdata)
   		k <- length(prior) 
		if (k==1) stop("The data belongs to one single group. A partition into at least two different groups is required\n")

   		sphdata <- newdata %*% object@scaling 
   		sphmeans <- object@means %*% object@scaling 
   		Mahdistovertwo <- apply(sphdata, 1, function(x) apply(sphmeans, 1, function(mu) (sum((mu-x)^2)/2)))
   		wghtdensities <- sweep(exp(-Mahdistovertwo),1,STATS=prior,FUN="*")
   		ncnst <- apply(wghtdensities,2,sum)  			# normalizing constants
   		posterior <- sweep(wghtdensities,2,STATS=ncnst,FUN="/")
   		clres <- apply(posterior, 2, function(pst) return(dimnames(sphmeans)[[1]][which.max(pst)]))
   		list(class=factor(clres,levels=dimnames(object@means)[[1]]),posterior=t(posterior))
	}
)

setMethod("show",					
	signature(object = "Idtlda"),
	function(object)
	{
		cat("Prior probabilities of groups:\n") ; print(object@prior) ; cat("\n")
		cat("Group means:\n") ; print(object@means) ; cat("\n")
		cat("Coefficients of linear discriminants:\n") ; print(object@scaling)
	}
)

setMethod("qda",
	signature(x = "IdtHetNMANOVA"),
	function(x,prior="proportions",selmodel=BestModel(H1res(x)),...)
	{
        	p <- 2*x@NIVar
   		nk <- as.numeric(table(x@grouping))
		N <- sum(nk)
     		k <- length(nk) 
		vnames <- colnames(x@H1res@mleNmuE)
		lev <- levels(x@grouping)
 		scaling <- array(dim=c(p,p,k),dimnames=list(vnames,paste("LD",1:p,sep=""),lev))
        	ldet <- numeric(k)
		if (is.character(selmodel))  selmodel <- sapply(selmodel,function(mod) which(mod==x@H1res@ModelNames))

		means <- x@H1res@mleNmuE
		if (prior[1]=="proportions") prior <- nk/N
		names(prior) <- lev
		for (g in 1:k)  {
	   		if (selmodel != 5) {
				W <- x@H1res@Configurations[[selmodel]]$mleSigE[,,g]
				scaling[,,g] <- backsolve(chol(W),diag(p))
				ldet[g] <- as.numeric(determinant(W,logarithm=TRUE)$modulus)/2
          		}
	   		else {
				Wd <- diag(x@H1res@Configurations[[5]]$mleSigE[,,g])
				scaling[,,g] <- diag(1/sqrt(Wd))
				ldet[g] <- sum(log(Wd))/2
           		}
        	}  

		new("Idtqda",prior=prior,means=means,scaling=scaling,ldet=ldet,lev=lev) 
	}
)

setMethod("qda",
	signature(x = "IData"),
	function(x, grouping, prior="proportions", tol=1.0e-4, subset=1:nrow(x), Config=1:5,SelCrit=c("AIC","BIC"))
	{
	   	SelCrit <- match.arg(SelCrit)
		if (length(subset) < nrow(x)) {
			x <- x[subset,]
			grouping <- grouping[subset]
		}
    		levels(grouping)[table(grouping)==0] <- NA
		if (length(levels(grouping))==1) 
			stop("The data belongs to one single group. A partition into at least two different groups is required\n")
		qda(MANOVA(x,grouping,Mxt="Het",Model="Normal",Config=Config,SelCrit=SelCrit,tol=tol),prior=prior)
	}
)

setMethod("predict",
	signature(object = "Idtqda"),
	function(object,newdata,prior=object@prior,...)
	{
		if (is(newdata,"IData")) newdata <- as.matrix(cbind(newdata@MidP,newdata@LogR))
   		if (is(newdata,"data.frame")) newdata <- as.matrix(newdata)
   		n <- nrow(newdata)
   		p <- ncol(newdata)
   		k <- length(prior) 
		grpnames <- dimnames(object@means)[[1]]
   		Mahdistovertwo <- matrix(nrow=k,ncol=n,dimnames=list(grpnames,rownames(newdata)))
   		for (g in 1:k)  {
   			sphdata <- newdata %*% object@scaling[,,g] 
   			sphmeans <- object@means[g,] %*% object@scaling[,,g]
   			Mahdistovertwo[g,] <- apply(sphdata, 1, function(x) apply(sphmeans, 1, function(mu) (sum((mu-x)^2)/2)))
   		} 
   		wghtdensities <- sweep(exp(sweep(-Mahdistovertwo,1,STATS=object@ldet,FUN="-")),1,STATS=prior,FUN="*")
   		ncnst <- apply(wghtdensities,2,sum)  			# normalizing constants
   		posterior <- sweep(wghtdensities,2,STATS=ncnst,FUN="/")
   		clres <- apply(posterior, 2, function(pst) return(grpnames[which.max(pst)]))
   		
		list(class=factor(clres,levels=grpnames),posterior=t(posterior))
      }
) 

setMethod("show",					
	signature(object = "Idtqda"),
	function(object)
	{
		cat("Prior probabilities of groups:\n") ; print(object@prior) ; cat("\n")
		cat("Group means:\n") ; print(object@means) ; cat("\n")
	}
)
