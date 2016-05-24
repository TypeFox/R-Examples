demo_hddc  <- function(){
	myEnv <- new.env()
	data(list="Crabs",envir = myEnv)
	myData <- get("Crabs",envir=myEnv)
	devAskNewPage(ask = FALSE)
	X <- as.matrix(myData[,-1])
	clx <- myData[,1]
	algorithm <- c("EM","CEM","SEM")
	initialization <- c("kmeans","random","param")
	model <- c("AKBKQKDK","ABQKD","AJBQD")
	while (1){
		while(!any((algo <- readline("Choose the algorithm:\n1: EM; 2: CEM; 3: SEM; q or quit: exit \n"))==as.character(1:3))){
			if(any(tolower(algo)%in%c("q","quit"))) return(invisible())
		}
		
		while(!any((init <- readline("Choose initialization:\n1: kmeans; 2: random; 3: param; q or quit: exit \n"))==as.character(1:3))){
			if(any(tolower(init)%in%c("q","quit"))) return(invisible())
		}
		
		while(!any((mod <- readline("Choose the model:\n1: AkBkQkDk; 2: ABQkD; 3: AjBQD; q or quit: exit \n"))==as.character(1:3))){
			if(any(tolower(mod)%in%c("q","quit"))) return(invisible())
		}
		cat("hddc(data,classes,model=\"",model[as.numeric(mod)],"\",init=\"",initialization[as.numeric(init)],"\",algo=\"",algorithm[as.numeric(algo)],"\")\n",sep="")
		demo_hddc_crabs(X,4,init=initialization[as.numeric(init)],algo=algorithm[as.numeric(algo)],model=model[as.numeric(mod)],ZZ=clx)
	}
}

demo_hddc_acp <- function(X,z,hd=NULL,...){
	Q <- hd$Q
	MU <- hd$mu
	if (is.matrix(Q)) {
		svg <- vector(mode='list',length=4)
		for (i in 1:4) svg[[i]] <- as.matrix(Q)
		Q <- svg
	}

	X <- as.matrix(X)
	p <- ncol(X)
	k <- max(z)
	
	CO <- cov(X)
	Z <- -eigen(CO,symmetric=T)$vectors
	coul <- 1:4
	
	V <- X%*%Z
	patch <- c(3,4,8,20)
	plot(V[,1],V[,2],col=z,pch=patch[z],...)
	
	#drawing the projective space (a line); matrix(10,p,1) is used only to have two points with the mean
	proj <- matrix(NA,k,p)
	for (i in 1:k) proj[i,] <- tcrossprod(Q[[i]])%*%matrix(10,p,1)+MU[i,]
	
	x <- proj%*%Z
	y <- MU%*%Z
	points(y[,1],y[,2],col=coul[1:k],pch=19,lwd=7)
	
	for (i in 1:k) {
		pente <- (x[i,2]-y[i,2])/(x[i,1]-y[i,1])
		oo <- x[i,2]-pente*x[i,1]
		xb <- (2*y[i,1]-sqrt(50^2/(pente^2+1)))/2
		xa <- (2*y[i,1]+sqrt(50^2/(pente^2+1)))/2
		lines(c(xa,xb),oo+pente*c(xa,xb),col=coul[i],type='l')
	}

}

demo_hddc_crabs <- function(DATA,k=4,model='AKBKQKD',threshold=0.2,method='C',algo='EM',itermax=50,eps=1e-2,init='kmeans',ZZ=NULL,min.individuals=2,noise.ctrl=1e-8,...){ 
	com_dim <- 1
	Mod <- c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD")
	p <- ncol(DATA)
	N <- nrow(DATA)
	t <- matrix(0,N,k)
	if(model%in%Mod[7:14]) method <- 1
	if (init=='param'){
		MU <- colMeans(DATA)
		prop <- rep(1/k,k)
		S <- crossprod(DATA-matrix(MU,N,p,byrow=TRUE))/N
		donnees <- eigen(S,symmetric=TRUE)
		ev <- donnees$values
		d <- if(is.numeric(method)) method else hdclassif_dim_choice(ev,N,method,threshold,FALSE,noise.ctrl)
		a <- ev[1:d]
		b <- sum(ev[(d[1]+1):p])/(p-d[1])
		
		Q <- donnees$vectors[,1:d]
		mu <- mvrnorm(k,MU,S)
		
		K <- diag((mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(DATA))+1/b*(diag(tcrossprod(mu))-2*mu%*%t(DATA)+2*(mu%*%Q)%*%(t(Q)%*%t(DATA))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
		
		t <- matrix(0,N,k)
		for (i in 1:k) t[,i]=1/rowSums(exp((K[i,]-t(K))/2))
	}
	else if (init=='kmeans') {
		mc <- match.call(expand.dots = FALSE)$...
		if (is.null(mc$algorithm)) alg="Hartigan-Wong"
		else alg=mc$algorithm
		if (is.null(mc$iter.max)) im=50
		else im=mc$iter.max
		if (is.null(mc$nstart)) nst=4
		else nst=mc$nstart
		cluster <- kmeans(DATA,k,iter.max=im,nstart=nst,algorithm=alg)$cluster
		for (i in 1:k) t[which(cluster==i),i] <- 1
	}
	else {
		t <- t(rmultinom(N,1,rep(1/k,k)))
		compteur=1
		while(min(colSums(t))<1 && (compteur <- compteur+1)<5) t <- t(rmultinom(N,1,rep(1/k,k)))
		if(min(colSums(t))<1) stop("Random initialization failed because of too many classes and too few observations")
	}
	
	likely <- c()
	I <- 0
	test <- Inf
	while ((I <- I+1)<=itermax && test>=eps){
		if (algo!='EM' && I!=1) t <- t2
		if (k>1 && (any(is.na(t)) || any(colSums(t>1/k)<min.individuals))) return(1)
		m <- hddc_m_step(DATA,k,t,model,threshold,method,noise.ctrl,1)
		t <- hddc_e_step(DATA,m)
		L <- t$L
		t <- t$t
		if (algo=='CEM') {
			t2 <- matrix(0,N,k)
			t2[cbind(1:N,max.col(t))] <- 1
		}
		else if(algo=='SEM') { 
			t2 <- matrix(0,N,k)
			for (i in 1:N)	t2[i,] <- t(rmultinom(1,1,t[i,]))
		}
		
		classes<-c()
		for (i in 1:N) classes[i]=which.max(t[i,])
		demo_hddc_acp(DATA,classes,m,xlab=paste('Iteration',I),ylab='',main="Clustering process",...)
		Sys.sleep(0.12)
		
		likely[I] <- L
		if (I!=1) test <- abs(likely[I]-likely[I-1])
	}
	
	cls <- max.col(t)
	ari = hddc_ari(ZZ,cls)
	cat("Adjusted Rand Index:",ari,"\n")
}

hdda  <- function(data,cls,model='AkjBkQkDk',graph=FALSE,d="Cattell",threshold=0.2,com_dim=NULL,show=TRUE,scaling=FALSE,cv.dim=1:10,cv.threshold=c(.001,.005,.05,1:9*0.1),cv.vfold=10,LOO=FALSE,noise.ctrl=1e-8){
	Mod <- c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD","ALL")
	Mod2 <- c("AKJBKQKDK","AKBKQKDK ","ABKQKDK  ","AKJBQKDK ","AKBQKDK  ","ABQKDK   ","AKJBKQKD ","AKBKQKD  ","ABKQKD   ","AKJBQKD  ","AKBQKD   ","ABQKD    ","AJBQD    ","ABQD     ")
	
	if (is.numeric(model)) model <- na.omit(Mod[model])
	else model <- toupper(model)
	num <- which(model%in%as.character(1:14))
	model[num] <- Mod[as.numeric(model[num])]
	if (any(!model%in%Mod)) stop("Invalid model name\n",call.=FALSE)
	mod_num <- c()
	for(i in 1:length(model)) mod_num[i] <- which(model[i]==Mod)
	mod_num <- sort(unique(mod_num))
	model <- Mod[mod_num]
	
	if(length(d)>1) stop("d cannot be a vector.\n",call.=FALSE)
	if(!is.numeric(d)) d <- toupper(d)
	
	if (d%in%c("CATTELL","C")) d <- "C"
	else if (d%in%c("BIC","B")) d <- "B"
	
	if(is.integer(d)) stop("d must be equal to \"BIC\" or \"Cattell\"",call.=FALSE)
	if (!is.numeric(threshold)) stop("The parameter 'threshold' must a double within 0 and 1.\n",call.=FALSE)
	if (threshold<0 | threshold>1) stop("The parameter 'threshold' must be within 0 and 1.\n",call.=FALSE)
	if (any(is.na(data))) stop("NA values are not allowed\n",call.=FALSE)
	if (nrow(data)!=length(cls)) stop ("The class vector does not fit with the data\n",call.=FALSE)
	if((length(model)>1 || (length(model)==1 && model=="ALL")) && d == "CV") stop("A specific model must be chosen for choosing the threshold/dimension with cross-validation.\n",call.=FALSE)
	if(any(model%in%Mod[c(1:6)]) && !d%in%c("C","B","CV")) stop("To use this model d must have the value: 'c', 'b' or 'cv'; it also may be an integer but only for common dimension model. See help for information.\n",call.=FALSE)
	if( any(model%in%Mod[7:14]) && !d%in%c("C","B","CV") && is.null(com_dim) ) stop("d must have the value: 'c', 'b' or 'cv'; it also may be an integer but only for common dimension model. See help for information.\n",call.=FALSE)
	
	
	save_d <- if(d=="B") "B" else "C"
	cls <- as.factor(cls)
	names <- levels(cls)
	z <- unclass(cls)
	data <- as.matrix(data)
	K <- max(z)
	if (scaling) {
		data <- scale(data)
		scaling <- list(mu=attr(data,"scaled:center"),sd=attr(data,"scaled:scale"))
	}
	else scaling <- NULL
	N <- nrow(data)
	p <- ncol(data)
	n <- table(z)
	
	if( (any(model%in%Mod[7:12]) && !is.null(com_dim) && com_dim > min(n,p)-1) )  stop("com_dim has to be lower or equal to ",min(n,p)-1, if(p>min(n)) paste(" because of the class",names[which.min(n)],"which has",min(n),"observations.\n") else paste(" because there are only",p,"dimensions.\n"),call.=FALSE)
	if(any(model%in%Mod[13:15]) && !is.null(com_dim) && com_dim > min(N,p)-1) stop("com_dim has to be lower or equal to ",min(N,p)-1, if(p<=N) paste(" because there are only",p,"dimensions.\n") else paste(" because there are only",N,"observations.\n"),call.=FALSE)
	
	
	if(LOO){
		if(d=="CV" && is.null(com_dim)) stop("To do LOO, d must be equal to 'c' or 'b', or you may select manually the dimension for common dimension models using 'com_dim'.\n",call.=FALSE)
		if(length(model)>1 || model=="ALL") stop("A specific model must be chosen for LOO.\n",call.=FALSE)
		
		if(model%in%Mod[1:14] && !is.null(com_dim)) d <- com_dim
		
		class <- rep(NA,N)
		posterior <- matrix(NA,N,K)
		for(i in 1:N){
			prms <- NULL
			try(prms <- hdda_prms(data[-i,],z[-i],model,threshold,d,names,noise.ctrl),silent=TRUE)
			if(!is.null(prms)) {
				res <- NULL
				try(res <- predict.hdc(prms,data[i,]),silent=TRUE)
				if(!is.null(res)){
					class[i] <- res$class
					posterior[i,] <- res$posterior
				}
			}
		}
		class <- factor(class,labels=names,levels=seq_along(names))
		return(list(class=class,posterior=posterior))
	}
	
	
	if(d=='CV'){
		d.max <- if(model%in%Mod[7:12]) min(n,p)-1 else min(N,p)-1
		cv.dim <- sort(cv.dim,decreasing=TRUE)
		cv.dim <- cv.dim[cv.dim<=d.max]
		if(length(cv.dim)==0) stop("cv.dim must be an integer stricly inferior \nto the dimension.\n",call.=FALSE)
		cv.threshold <- sort(cv.threshold) 
		cv.threshold <- cv.threshold[cv.threshold>=0  & cv.threshold<=1]
		if(length(cv.threshold)==0) stop("cv.threshold must be a float within 0 and 1.\n",call.=FALSE)
		cv.vfold <- if(cv.vfold<N) cv.vfold else N
		
		u <- sample(1:N)
		ind <- c()
		for(i in 1:cv.vfold) ind[i] <- if(i==1) floor(N/cv.vfold) else floor((N-sum(ind))/(cv.vfold+1-i))
		fin <- cumsum(ind)
		debut <- c(1,fin[-cv.vfold]+1)
		
		if(model%in%Mod[7:14]) {
			n_cv <- length(cv.dim)
			cv.threshold <- rep(.5,n_cv)
		}
		else{
			n_cv <- length(cv.threshold)
			cv.dim <- rep("C",n_cv)
		}
		
		res <- fails <- rep(0,n_cv)
		N2 <- rep(N,n_cv)
		for(j in 1:cv.vfold){
			ind <- u[debut[j]:fin[j]]
			prms <- NULL
			i <- 0
			while((i <- i+1)<=n_cv && is.null(prms) ){
				try(prms <- hdda_prms(data[-ind,],z[-ind],model,cv.threshold[i],cv.dim[i],names,noise.ctrl),silent=TRUE)
				if(!is.null(prms)) try(res[i] <-  res[i] + sum(predict.hdc(prms,data[ind,])$class==cls[ind]), silent=TRUE)
				else {
					N2[i] <- N2[i]-length(ind)
					fails[i] <- fails[i] + 1
				}
			}
			
			if(i<=n_cv) for(i in i:n_cv){
				if(model%in%Mod[1:6]) d <- hdclassif_dim_choice(prms$ev,as.vector(table(z[-ind])),"C",cv.threshold[i],FALSE,noise.ctrl)
				else d <- rep(cv.dim[i],K)
				if(model%in%Mod[13:14]) prms$Q <- prms$Q[,1:d[1]]
				else for(ii in 1:K) if(prms$d[ii]>1) prms$Q[[ii]] <- prms$Q[[ii]][,1:d[ii]]
				prms$d <- d	
				prms_bis <- hdda_prms_bis(model,prms,p)
				try(res[i] <-  res[i] + sum(predict.hdc(prms_bis,data[ind,])$class==cls[ind]), silent=TRUE)
			}
		}
		
		if(show){
			if(model%in%Mod[7:14]) cat("\t  Model   \t dim\t CV\n")
			else cat("\t  Model   \tthreshold\t CV\n")
			for(i in n_cv:1){
				if(model%in%Mod[7:14]) cat('\t',Mod2[model==Mod],'\t',cv.dim[i],"\t",res[i]/N2[i]*100,if(fails[i]>0) paste("  Info: failed",fails,"times"),'\n')
				else cat('\t',Mod2[model==Mod],'\t',cv.threshold[i],"\t\t",res[i]/N2[i]*100,if(fails[i]>0) paste("  Info: failed",fails,"times"),'\n')
			}
		}
		
		res <- res/N2*100
		res <- res[n_cv:1]
		cv.dim <- cv.dim[n_cv:1]
		cv.threshold <- cv.threshold[n_cv:1]
		if(model%in%Mod[7:14]){
			d <- com_dim <- cv.dim[which.max(res)]
			if(show) cat("Best dimension with respect to the CV results: ",d,".\n",sep="")
			if(graph){
				barplot(res-100/K,names.arg=cv.dim,offset=100/K,col="blue", xlab="Dimensions",ylab="Correct classification rate",axes=FALSE, main=paste("Cross-Validation\n(chosen dim=",d,")",sep=""))
				axis(2,at=floor(100/K+(max(res)-100/K)/5*0:5))
			}
		}
		else{
			d <- "C"
			threshold <- cv.threshold[which.max(res)]
			if(show) cat("Best threshold with respect to the CV results: ",threshold,".\n",sep="")
			if(graph){
				barplot(res-100/K,names.arg=cv.threshold,offset=100/K,col="blue", xlab="Thresholds",ylab="Correct classification rate",axes=FALSE, main=paste("Cross-Validation\nthreshold=",threshold,sep=""))
				axis(2,at=floor(100/K+(max(res)-100/K)/5*0:5))
			}
		}
	}	
	
	if(length(model)>1){
		e <- vector(mode="list",length=max(mod_num))
		BIC <- c()
		
		for(i in mod_num){
			e[[i]] <- hdda_prms(data,z,Mod[i],threshold,d,names,noise.ctrl,com_dim)
			BIC[i] <- hdclassif_bic(e[[i]],p)
		}
		
		prms <- e[[which.max(BIC)]]
		prms$BIC <- max(BIC,na.rm=TRUE)
		prms$scaling <- scaling
		prms$threshold <- if(save_d=="C") threshold else NULL
		
		if(show){
  			cat(" # :\t  Model  \t     BIC\n")
   			for(i in mod_num){
   				if(i<10) cat(' ')
				wng <- if(any(e[[i]]$b<10e-6) | any(e[[i]]$a<10e-6,na.rm=TRUE)) "info: b < 10e-6" else ""
				cat(i,':\t',Mod2[i],'\t',BIC[i],wng,'\n')
   			}
			cat("\nSELECTED: Model ",prms$model,", BIC=",prms$BIC,".\n",sep="")
   		}
		
		if(graph){			
			BIC <- BIC[!is.na(BIC)]
			min_b=min(BIC[BIC!=-Inf])
			max_b=max(BIC,na.rm=TRUE)
			BIC[BIC==-Inf] <- min_b
			barplot(BIC-min_b,names.arg=mod_num,offset=min_b,col="blue", xlab="models",ylab="BIC",axes=FALSE, main=paste("BIC for all models\n(chosen model=",prms$model,")",sep=""))
			axis(2,at=floor(min_b+(max_b-min_b)/5*0:5))
		}
		class(prms) <- 'hdc'
		return(prms)
	}
	else if(model=="ALL"){
		e <- vector(mode="list",length=14)
		BIC <- c()
		
		#models with var dim
		e[[1]] <- hdda_prms(data,z,Mod[1],threshold,d,names,noise.ctrl)
		for (i in 2:6) e[[i]] <- hdda_prms_bis(Mod[i],e[[1]],p)

		#models with common dim	
		e[[7]] <- hdda_prms(data,z,Mod[7],threshold,d,names,noise.ctrl,com_dim)
		for (i in 8:12) e[[i]] <- hdda_prms_bis(Mod[i],e[[7]],p)
		
		#models 13 and 14: common var/covar matrix
		e[[13]] <- hdda_prms(data,z,Mod[13],threshold,d,names,noise.ctrl,com_dim)
		e[[14]] <- hdda_prms_bis(Mod[14],e[[13]],p)
		
		#BIC calculation
		for(i in 1:14) BIC[i] <- hdclassif_bic(e[[i]],p)
   
   		prms <- e[[which.max(BIC)]]
		prms$BIC <- max(BIC,na.rm=TRUE)
		prms$scaling <- scaling
		prms$threshold <- if(save_d=="C") threshold else NULL
		
		if(show){
  			cat(" # :\t  Model  \t     BIC\n")
   			for(i in 1:14){
   				if(i<10) cat(' ')
				wng <- if(any(e[[i]]$b<10e-6) | any(e[[i]]$a<10e-6,na.rm=TRUE)) "info: b < 10e-6" else ""
				cat(i,':\t',Mod2[i],'\t',BIC[i],wng,'\n')
   			}
			cat("\nSELECTED: Model ",prms$model,", BIC=",prms$BIC,".\n",sep="")
   		}
		if(graph){
			min_b <- min(BIC[BIC!=-Inf])
			max_b <- max(BIC)
			BIC[BIC==-Inf] <- min_b
			barplot(BIC-min_b,names.arg=1:14,offset=min_b,col="blue", xlab="models",ylab="BIC",axes=FALSE, main=paste("BIC for all models\n(chosen model=",prms$model,")",sep=""))
			axis(2,at=floor(min_b+(max_b-min_b)/5*0:5))
		}
		class(prms) <- 'hdc'
		return(prms)
	}
	else {
		prms <- hdda_prms(data,z,model,threshold,d,names,noise.ctrl,com_dim)
		prms$BIC <- hdclassif_bic(prms,p)
		prms$scaling <- scaling
		prms$threshold <- if(save_d=="C") threshold else NULL
		
		class(prms) <- 'hdc'
		return(prms)
	}
}

hddc  <- function(data,K=1:10,model=c("AkjBkQkDk"),threshold=0.2,com_dim=NULL,itermax=60,eps=1e-3,graph=FALSE,algo='EM',d="Cattell",init='kmeans',show=TRUE,mini.nb=c(5,10),scaling=FALSE,min.individuals=2,noise.ctrl=1e-8,...){
	Mod <- c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD")
	Mod2 <- c("AKJBKQKDK","AKBKQKDK ","ABKQKDK  ","AKJBQKDK ","AKBQKDK  ","ABQKDK   ","AKJBKQKD ","AKBKQKD  ","ABKQKD   ","AKJBQKD  ","AKBQKD   ","ABQKD    ","AJBQD    ","ABQD     ")
	Alg <- c('EM','CEM','SEM')
	Init <- c('random','kmeans','mini-em','param')
	algo=toupper(algo)
	if(length(init)>1){
		init <- unclass(init)
		if(any(K!=max(init))) stop("The number of class of K and of the initialization vector are different\n")
		if(length(init)!=nrow(data)) stop("The size of the initialization vector is different of the size of the data\n")
	}
	
	if(!is.numeric(d)) d <- toupper(d)
	if (d%in%c("CATTELL","C")) d <- "C"
	else if (d%in%c("BIC","B")) d <- "B"
	save_d <- d
	
	if(length(model)==1 && toupper(model)=="ALL") model <- 1:14
	if (is.numeric(model)) model <- na.omit(Mod[model])
	else model <- toupper(model)
	num <- which(model%in%as.character(1:14))
	model[num] <- Mod[as.numeric(model[num])]
	if (any(!model%in%Mod)) stop("Invalid model name\n")
	mod_num <- c()
	for(i in 1:length(model)) mod_num[i] <- which(model[i]==Mod)
	mod_num <- sort(unique(mod_num))
	model <- Mod[mod_num]
	
	if(is.integer(d)) stop("d must be equal to \"BIC\" or \"Cattell\"")
	if (any(model%in%Mod[7:14]) && is.numeric(d) && d>ncol(data)) stop("d must be strictly inferior to the dimension, \nwhich is in this case ",ncol(data),'\n')
	if (!is.numeric(min.individuals) || min.individuals<2) stop("The minimum population control variable must be superior or equal to 2.\n")
	if (length(init)==1 && !any(init==Init)) stop("Invalid initialization name\n")
	if (is.numeric(threshold)==0 || threshold<=0 || threshold>=1) stop("The parameter 'threshold' must be a double strictly within ]0,1[\n")
	if (!any(Alg==algo)) stop("Invalid algorithm name\n")
	if (length(init)==1 && init=='param' && nrow(data)<ncol(data)) stop("The 'param' initialization can't be done when N<p\n")
	if (any(is.na(data))) stop("NA values are not supported\n")
	#if (length(init)==1 && init=='param' && library(MASS,logical.return=TRUE)==FALSE) stop("You need the library MASS to use the 'param' initialization\n") 
	if (length(init)==1 && init=='mini-em' && (length(mini.nb)!=2 | is.numeric(mini.nb)!=1)) stop("The parameter mini.nb must be a vector of length 2 with integers\n")
	if (!is.numeric(K) || min(K)<1) stop("K must be a vector of positive integers\n")

	data <- as.matrix(data)
	if (scaling) {
		data <- scale(data)
		scaling <- list(mu=attr(data,"scaled:center"),sd=attr(data,"scaled:scale"))
	}
	else scaling <- NULL
	BIC <- c()
	p <- ncol(data)
	e <- vector(mode="list",length=length(K))
	if (show) cat('\t  Model  \t   K\t   BIC\n')
	nm <- length(model)
	ind <- 1
	for (i in (K <- floor(sort(K)))){
		if (i==1){
			e[[1]] <- hddc_main(data,1,"AKJBKQKDK",threshold,d,algo,itermax,eps,init,mini.nb,min.individuals,noise.ctrl,...)
			BIC[1:nm] <- hdclassif_bic(e[[1]],p)
			if (show) cat('\t',"ALL      ",'\t',1,'\t',BIC[1],'\n')
			ind <- nm+1
		}
		else {
			for (M in model){
				e[[ind]] <- hddc_main(data,i,M,threshold,d,algo,itermax,eps,init,mini.nb,min.individuals,noise.ctrl,com_dim,...)
				if (length(e[[ind]])==1){
					if (show) cat('\t',Mod2[which(Mod==M)],'\t',i,'\t',"STOPPED: pop<min.indiv.\n")
					BIC[ind] <- -Inf
				}
				else {
					if(M%in%Mod[13:14]) BIC[ind] <- hdclassif_bic(e[[ind]],p,data)
					else BIC[ind] <- hdclassif_bic(e[[ind]],p)
					if (show) cat('\t',Mod2[which(Mod==M)],'\t',i,'\t',BIC[ind],'\n')
				}
				ind <- ind+1
			}
		}
	}
	
	if (max(BIC)==-Inf) return(NULL)
	if(graph){
		g <- matrix(BIC,nm,length(K))
		g[g==-Inf] <- NA
		if (length(K)==1) plot(as.factor(model),g,ylab="BIC",main=paste("K =",K),xlab="model")
		else{
			plot(K,g[1,],type='o',ylim=c(min(g,na.rm=TRUE),max(BIC)),pch=1,ylab="BIC")
			if (nm>1) for (i in 2:nm) lines(K,g[i,],col=i,pch=i,type='o',lty=i)
			legend(min(K,na.rm=TRUE),max(BIC),model,col=1:nm,pch=1:nm,bty="n",lwd=1,cex=0.85,lty=1:nm)
		}
	}
	prms <- e[[which.max(BIC)]]
	if (show & (length(model)>1 | length(K)>1)) cat("\nSELECTED: model ",prms$model," with ",prms$K," clusters, BIC=",max(BIC),".\n",sep="")
	prms$BIC <- max(BIC)
	prms$scaling <- scaling
	prms$threshold <- if(save_d=="C") threshold else NULL
	class(prms) <- 'hdc'
	return(prms)
}

hdclassif_dim_choice <- function(ev,n,method,threshold,graph,noise.ctrl){
	N <- sum(n)
	prop <- n/N
	K <- if(is.matrix(ev)) nrow(ev) else 1
	if(is.matrix(ev) && K>1){
  		p <- ncol(ev)
		if(method=="C"){
			dev <- abs(apply(ev,1,diff))
			max_dev <- apply(dev,2,max,na.rm=TRUE)
			dev <- dev/rep(max_dev,each=p-1)
			d <- apply((dev>threshold)*(1:(p-1))*t(ev[,-1]>noise.ctrl),2,which.max)
			
			if(graph){
				par(mfrow=c(K*(K<=3)+2*(K==4)+3*(K>4 && K<=9)+4*(K>9),1+floor(K/4)-1*(K==12)+1*(K==7)))
				for(i in 1:K){
					sub1 <- paste("Class #",i,", d",i,"=",d[i],sep="")
					Nmax <- max(which(ev[i,]>noise.ctrl))-1
					plot(dev[1:(min(d[i]+5,Nmax)),i],type="l",col="blue",main=paste("Cattell's Scree-Test\n",sub1,sep=""),ylab=paste("threshold =",threshold),xlab="Dimension",ylim=c(0,1.05))
					abline(h=threshold,lty=3)  	
					points(d[i],dev[d[i],i],col='red')
				}
			}
		}
		else if(method=="B"){
			d <- rep(0,K)
			if(graph) par(mfrow=c(K*(K<=3)+2*(K==4)+3*(K>4 && K<=9)+4*(K>9),1*(1+floor(K/4)-1*(K==12)+1*(K==7))))
			
			for (i in 1:K) {
				B <- c()
				Nmax <- max(which(ev[i,]>noise.ctrl))-1
				p2 <- sum(!is.na(ev[i,]))
				Bmax <- -Inf
				for (kdim in 1:Nmax){
					if ((d[i]!=0 & kdim>d[i]+10)) break
					a <- sum(ev[i,1:kdim])/kdim
					b <- sum(ev[i,(kdim+1):p2])/(p2-kdim)
					if (b<0 | a<0) B[kdim] <- -Inf
					else{
						L2 <- -1/2*(kdim*log(a)+(p2-kdim)*log(b)-2*log(prop[i])+p2*(1+1/2*log(2*pi))) * n[i]
						B[kdim] <- 2*L2 - (p2+kdim*(p2-(kdim+1)/2)+1) * log(n[i])
					}
					if ( B[kdim]>Bmax ){
						Bmax <- B[kdim]
						d[i] <- kdim
					}
				}
				
				if(graph){
					plot(B,type='l',col=4,main=paste("class #",i,", d=",d[i],sep=''),ylab='BIC',xlab="Dimension")
					points(d[i],B[d[i]],col=2)
				}
			}
		}
	}
  	else{
		ev <- as.vector(ev)
		p <- length(ev)
		if(method=="C"){
			dvp <- abs(diff(ev))
			Nmax <- max(which(ev>noise.ctrl))-1
			if (p==2) d <- 1
			else d <- max(which(dvp[1:Nmax]>=threshold*max(dvp[1:Nmax])))
			diff_max <- max(dvp[1:Nmax])
			
			if(graph){
				plot(dvp[1:(min(d+5,p-1))]/diff_max,type="l",col="blue",main=paste("Cattell's Scree-Test\nd=",d,sep=''),ylab=paste("threshold =",threshold,sep=' '),xlab='Dimension',ylim=c(0,1.05))
				abline(h=threshold,lty=3)	
				points(d,dvp[d]/diff_max,col='red')
			}
		}
		else if(method=="B"){
			d <- 0
			Nmax <- max(which(ev>noise.ctrl))-1
			B <- c()
			Bmax <- -Inf
			for (kdim in 1:Nmax){
				if (d!=0 && kdim>d+10) break
				a <- sum(ev[1:kdim])/kdim
				b <- sum(ev[(kdim+1):p])/(p-kdim)
				if (b<=0 | a<=0) B[kdim] <- -Inf
				else{
					L2 <- -1/2*(kdim*log(a)+(p-kdim)*log(b)+p*(1+1/2*log(2*pi)))*N
					B[kdim] <- 2*L2 - (p+kdim*(p-(kdim+1)/2)+1)*log(N)
				}
				if ( B[kdim]>Bmax ){
					Bmax <- B[kdim]
					d <- kdim
				}
			}
			
			if(graph){
				plot(B,type='l',col=4,main=paste("BIC criterion\nd=",d,sep=''),ylab='BIC',xlab="Dimension")
				points(d,B[d],col=2)
			}
		}
  	}
	return(d)
}

hdclassif_bic  <- function(par,p,data=NULL){
	model <- par$model
	K <- par$K
	d <- par$d
	b <- par$b
	a <- par$a
	mu <- par$mu
	N <- par$N
	prop <- par$prop
	
	if(length(b)==1){
		#update of b to set it as variable dimension models
		eps <- sum(prop*d)
		n_max <- if(model%in%c("ABQD","AJBQD")) length(par$ev) else ncol(par$ev)
		b <- b*(n_max-eps)/(p-eps)
		b <- rep(b,length=K)
	}	
	if (length(a)==1) a <- matrix(a,K,max(d))
	else if (length(a)==K) a <- matrix(a,K,max(d))
	else if (model=='AJBQD') a <- matrix(a,K,d[1],byrow=TRUE)
	
	if(min(a,na.rm=TRUE)<=0 | any(b<0)) return(-Inf)
	
	if (is.null(par$loglik)){
		som_a <- c()
		for (i in 1:K) som_a[i] <- sum(log(a[i,1:d[i]]))
		L <-  -1/2*sum(prop * (som_a + (p-d)*log(b) - 2*log(prop) + p*(1+log(2*pi))))*N
	}
	else if (model%in%c("ABQD","AJBQD")){
		Q <- rep(list(par$Q),K)
		K_pen <- matrix(0,K,N)
		for (i in 1:K) {
			s <- sum(log(a[i,1:d[i]]))
			X <- data-matrix(mu[i,],N,p,byrow=TRUE)
			proj <- (X%*%Q[[i]])%*%t(Q[[i]])
			A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
			B <- X-proj
			K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
		}
		A <- -1/2*t(K_pen)
		L <- sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
	}
	else L <- par$loglik[length(par$loglik)]
	
	
	ro <- K*p+K-1
	tot <- sum(d*(p-(d+1)/2))
	D <- sum(d)
	d <- d[1]
	to <- d*(p-(d+1)/2)
	if (model=='AKJBKQKDK') m <- ro+tot+2*K+D
	else if (model=='AKBKQKDK') m <- ro+tot+3*K
	else if (model=='ABKQKDK') m <- ro+tot+2*K+1
	else if (model=='AKJBQKDK') m <- ro+tot+K+D+1
	else if (model=='AKBQKDK') m <- ro+tot+2*K+1
	else if (model=='ABQKDK') m <- ro+tot+K+2
	else if (model=='AKJBKQKD') m <- ro+K*(to+d+1)+1
	else if (model=='AKBKQKD') m <- ro+K*(to+2)+1
	else if (model=='ABKQKD') m <- ro+K*(to+1)+2
	else if (model=='AKJBQKD') m <- ro+K*(to+d)+2
	else if (model=='AKBQKD') m <- ro+K*(to+1)+2
	else if (model=='ABQKD') m <- ro+K*to+3
	else if (model=='AJBQD') m <- ro+to+d+2
	else if (model=='ABQD') m <- ro+to+3
	bic <- -2*L+m*log(N)
	return(-bic)
}

hdda_prms  <- function(data,cls,model,threshold,method,kname,noise.ctrl,com_dim=NULL){
	p <- ncol(data)
	N <- nrow(data)
	K <- max(cls)
	com_ev <- NULL
	info <- NULL
	n <- as.vector(table(cls))
	prop <- matrix(n/N,1,K,dimnames=list(c(''),"Prior probabilities of groups:"=kname))
	
	mu <- matrix(rowsum(data,cls)/n,K,p,dimnames=list("Class"=kname,"Group means:"=paste('V',1:p,sep='')))

	#Calculation of Var/covar matrices and of eigenvectors
	
	if( model%in%c("AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD") ){
		if (N<p) {
			Y <- matrix(0,N,p)
			for (i in 1:K) Y[which(cls==i),] <- (data[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/sqrt(N)
			if(model%in%c("AJBQD","ABQD")) donnees <- eigen(tcrossprod(Y),symmetric=TRUE)
			else donnees <- eigen(tcrossprod(Y),symmetric=TRUE,only.values=TRUE)
		}
		else{
			W <- matrix(0,p,p)
			for (i in 1:K) W <- W + prop[i]*crossprod(data[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/n[i]
			if(model%in%c("AJBQD","ABQD")) donnees <- eigen(W,symmetric=TRUE)
			else donnees <- eigen(W,symmetric=TRUE,only.values=TRUE)
		}	
		ev <- com_ev <- donnees$values
	}
	
	if(!model%in%c("AJBQD","ABQD")){
		if(any(n<p)) Y <- vector(mode='list',length=K)
		Q <- vector(mode='list',length=K)
		ev <- matrix(NA,K,min(max(n),p))
		for(i in which(n<p)){
			Y[[i]] <- (data[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/sqrt(n[i])
			donnees <- eigen(tcrossprod(Y[[i]]),symmetric=TRUE)
			ev[i,1:n[i]] <- donnees$values
			Q[[i]] <- donnees$vectors
		}
		for(i in which(n>=p)){
			donnees <- eigen(crossprod(data[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/n[i],symmetric=TRUE)
			ev[i,] <- donnees$values
			Q[[i]] <- donnees$vectors
		}
	}
	
	#Intrinsic dimension calculations + graphical display
	
	if(model%in%c("AJBQD","ABQD")){
		if(!is.null(com_dim)) method <- com_dim
		if(method%in%c("C","B")) method <- hdclassif_dim_choice(com_ev,n,method,threshold,FALSE,noise.ctrl)
		d <- rep(method,K)
	}
	else if (model%in%c('AKJBKQKD','AKBKQKD','ABKQKD','AKJBQKD','AKBQKD','ABQKD')){
		if(!is.null(com_dim)) method <- com_dim
		if(method%in%c("C","B")) method <- hdclassif_dim_choice(com_ev,n,method,threshold,FALSE,noise.ctrl)
		d <- rep(method,K)
		if( d[1]>min(n,p)-1 ) {
			d[] <- min(n,p)-1
			info <- paste("Information: d has been lowered to",d[1],"because of the class",kname[which.min(n)],"which has",min(n),"observations.")
		}
		dmax <- if(any(ev<noise.ctrl,na.rm=TRUE)) max(min(unlist(apply(ev<noise.ctrl,1,which)))-2,1) else Inf
		if(d[1] > dmax) d[] <- dmax
	}
	else d <- hdclassif_dim_choice(ev,n,method,threshold,FALSE,noise.ctrl)
	
	#Setup of Qi matrices
	
	if (model%in%c("AJBQD","ABQD")){
		if (N>=p) {
			Q <- matrix(donnees$vectors[,1:d[1]],p,d[1])
		}
		else {
			Q <- matrix(t(Y)%*%donnees$vectors[,1:d[1]],p,d[1])
			normalise <- c()
			for(i in 1:d[1]) normalise[i] <- as.double(crossprod(Q[,i]))
			Q <- Q/matrix(sqrt(normalise),p,d[1],byrow=TRUE)
		}
	}
	else{
		for(i in which(n>=p)){
			Q[[i]] <- matrix(Q[[i]][,1:d[i]],p,d[i])
		}
		for(i in which(n<p)){
			Q[[i]] <- t(Y[[i]])%*%(Q[[i]][,1:d[i]])
			normalise <- c()
			for (j in 1:d[i]) normalise[j] <- as.double(crossprod(as.matrix(Q[[i]][,j])))
			Q[[i]] <- Q[[i]]/matrix(sqrt(normalise),p,d[i],byrow=TRUE)
		}
	}
	
	#Calculation of the remaining parameters of the selected model
	
	if ( model%in%c('AKJBKQKDK','AKJBQKDK','AKJBKQKD','AKJBQKD') ){
		ai <- matrix(NA,K,max(d),dimnames=list("Class"=kname,"Akj:"=paste("a",1:max(d),sep='')))
		for (i in 1:K) ai[i,1:d[i]] <- ev[i,1:d[i]]
	}
	else if ( model%in%c('AKBKQKDK','AKBQKDK' ,'AKBKQKD','AKBQKD') ){
		ai <- matrix(NA,1,K,dimnames=list(c("Ak:"),kname))
		for (i in 1:K) ai[i] <- sum(ev[i,1:d[i]])/d[i]
	}
	else if (model=="AJBQD"){
		ai <- matrix(ev[1:d[1]],1,d[1],dimnames=list(c("Aj:"),paste('a',1:d[1],sep='')))
	}
	else if (model=="ABQD"){
		ai <- matrix(sum(ev[1:d[1]])/d[1],dimnames=list(c("A:"),c('')))
	}
	else {
		a <- 0
		eps <- sum(prop*d)
		for (i in 1:K) a <- a + sum(ev[i,1:d[i]])*prop[i]
		ai <- matrix(a/eps,dimnames=list(c("A:"),c('')))
	}
	
	if ( model%in%c('AKJBKQKDK','AKBKQKDK','ABKQKDK','AKJBKQKD','AKBKQKD','ABKQKD') ){
		bi <- matrix(NA,1,K,dimnames=list(c("Bk:"),kname))
		for(i in which(n>=p)) bi[i] <- sum(ev[i,(d[i]+1):p])/(p-d[i])
		for(i in which(n<p)) bi[i] <- sum(ev[i,(d[i]+1):n[i]])/(p-d[i])
	}
	else if ( model%in%c("ABQD","AJBQD") ){
		if (N>=p) bi <- matrix(sum(ev[(d[1]+1):p])/(p-d[1]),dimnames=list(c("B:"),c('')))
		else bi <- matrix(sum(ev[(d[1]+1):N])/(N-d[1]),dimnames=list(c("B:"),c('')))
	}
	else{
		b <- 0
		eps <- sum(prop*d)
		for(i in which(n>=p)) b <- b + sum(ev[i,(d[i]+1):p])*prop[i]
		for(i in which(n<p)) b <- b + sum(ev[i,(d[i]+1):n[i]])*prop[i]
		bi <- matrix(b/(min(max(n),p)-eps),dimnames=list(c("B:"),c('')))
	}
	d <- matrix(d,1,K,dimnames=list(c('dim:'),"Intrinsic dimensions of the classes:"=kname))
	class(prop) <- class(mu) <- class(ai) <- class(bi) <- class(d) <- class(ev) <- "hd"
	list(model=model,K=K,d=d,a=ai,b=bi,mu=mu,prop=prop,ev=ev,Q=Q,kname=kname,info=info,N=N,com_ev=com_ev)
}

hdda_prms_bis  <- function(model,par,p){
	N <- par$N
	K <- par$K
	ev <- par$ev
	d <- par$d
	kname <- par$kname
	prop <- par$prop
	n <- prop*N
	
	if ( model%in%c('AKJBKQKDK','AKJBQKDK','AKJBKQKD','AKJBQKD') ){
		ai <- matrix(NA,K,max(d),dimnames=list("Class"=kname,"Akj:"=paste("a",1:max(d),sep='')))
		for (i in 1:K) ai[i,1:d[i]] <- ev[i,1:d[i]]
	}
	else if ( model%in%c('AKBKQKDK','AKBQKDK' ,'AKBKQKD','AKBQKD') ){
		ai <- matrix(NA,1,K,dimnames=list(c("Ak:"),kname))
		for (i in 1:K) ai[i] <- sum(ev[i,1:d[i]])/d[i]
	}
	else if (model=="AJBQD"){
		ai <- matrix(ev[1:d[1]],1,d[1],dimnames=list(c("Aj:"),paste('a',1:d[1],sep='')))
	}
	else if (model=="ABQD"){
		ai <- matrix(sum(ev[1:d[1]])/d[1],dimnames=list(c("A:"),c('')))
	}
	else {
		a <- 0
		eps <- sum(prop*d)
		for (i in 1:K) a <- a + sum(ev[i,1:d[i]])*prop[i]
		ai <- matrix(a/eps,dimnames=list(c("A:"),c('')))
	}
	
	if ( model%in%c('AKJBKQKDK','AKBKQKDK','ABKQKDK','AKJBKQKD','AKBKQKD','ABKQKD') ){
		bi <- matrix(NA,1,K,dimnames=list(c("Bk:"),kname))
		for(i in which(n>=p)) bi[i] <- sum(ev[i,(d[i]+1):p])/(p-d[i])
		for(i in which(n<p)) bi[i] <- sum(ev[i,(d[i]+1):n[i]])/(p-d[i])
	}
	else if ( model%in%c("ABQD","AJBQD") ){
		if (N>=p) bi <- matrix(sum(ev[(d[1]+1):p])/(p-d[1]),dimnames=list(c("B:"),c('')))
		else bi <- matrix(sum(ev[(d[1]+1):N])/(N-d[1]),dimnames=list(c("B:"),c('')))
	}
	else{
		b <- 0
		eps <- sum(prop*d)
		for(i in which(n>=p)) b <- b + sum(ev[i,(d[i]+1):p])*prop[i]
		for(i in which(n<p)) b <- b + sum(ev[i,(d[i]+1):n[i]])*prop[i]
		bi <- matrix(b/(min(max(n),p)-eps),dimnames=list(c("B:"),c('')))
	}
	class(ai) <- class(bi) <- "hd"

	list(model=model,K=K,d=d,a=ai,b=bi,mu=par$mu,prop=par$prop,ev=ev,Q=par$Q,kname=par$kname,info=par$info,N=N,com_ev=par$com_ev)
}

hddc_main <- function(DATA,K,model,threshold,method,algo,itermax,eps,init,mini.nb,min.individuals,noise.ctrl,com_dim=NULL,...){ 
	Mod <- c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD")
	p <- ncol(DATA)
	N <- nrow(DATA)
	com_ev <- NULL
	if ( any(model==Mod[7:14]) ){
		MU <- colMeans(DATA)
		if (N<p) {
			Y <- (DATA-matrix(MU,N,p,byrow=TRUE))/sqrt(N)
			YYt <- tcrossprod(Y)
			com_ev <- eigen(YYt,symmetric=TRUE,only.values=TRUE)$values
		}
		else{
			S <- crossprod(DATA-matrix(MU,N,p,byrow=TRUE))/N
			com_ev <- eigen(S,symmetric=TRUE,only.values=TRUE)$values
		}
		if(is.null(com_dim)) com_dim <- hdclassif_dim_choice(com_ev,N,method,threshold,FALSE,noise.ctrl)
	}
	if (K>1){
		t <- matrix(0,N,K)
		if(is.numeric(init)==1 | length(init)>1) {
			name <- unique(init)
			for (i in 1:K) t[which(init==name[i]),i] <- 1
		}
		else if (init=='param'){
			MU <- colMeans(DATA)
			prop <- rep(1/K,K)
			S <- crossprod(DATA-matrix(MU,N,p,byrow=TRUE))/N
			donnees <- eigen(S,symmetric=TRUE)
			ev <- donnees$values
			d <- if(is.numeric(method)) method else hdclassif_dim_choice(ev,N,method,threshold,FALSE,noise.ctrl)
			a <- ev[1:d]
			b <- sum(ev[(d[1]+1):p])/(p-d[1])
			
			Q <- donnees$vectors[,1:d]
			mu <- mvrnorm(K,MU,S)
			
			K_pen <- diag((mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(DATA))+1/b*(diag(tcrossprod(mu))-2*mu%*%t(DATA)+2*(mu%*%Q)%*%(t(Q)%*%t(DATA))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
			
			t <- matrix(0,N,K)
			for (i in 1:K) t[,i]=1/rowSums(exp((K_pen[i,]-t(K_pen))/2))
		}
		else if (init=='kmeans') {
			mc <- match.call(expand.dots = FALSE)$...
			if (is.null(mc$algorithm)) alg="Hartigan-Wong"
			else alg=mc$algorithm
			if (is.null(mc$iter.max)) im=50
			else im=mc$iter.max
			if (is.null(mc$nstart)) nst=4
			else nst=mc$nstart
			cluster <- kmeans(DATA,K,iter.max=im,nstart=nst,algorithm=alg)$cluster
			for (i in 1:K) t[which(cluster==i),i] <- 1
		}
		else if (init=='mini-em'){
			prms_best <- 1
			for (i in 1:mini.nb[1]){
				prms <- hddc_main(DATA,K,model,threshold,method,algo,mini.nb[2],0,'random',mini.nb,min.individuals,noise.ctrl,com_dim)
				if(length(prms)!=1){
					if (length(prms_best)==1) prms_best <- prms
					else if (prms_best$loglik[length(prms_best$loglik)]<prms$loglik[length(prms$loglik)]) prms_best <- prms
				}
			}
			if (length(prms_best)==1) return(1)
			t <- prms_best$posterior
		}
		else {
			t <- t(rmultinom(N,1,rep(1/K,K)))
			compteur=1
			while(min(colSums(t))<1 && (compteur <- compteur+1)<5) t <- t(rmultinom(N,1,rep(1/K,K)))
			if(min(colSums(t))<1) stop("Random initialization failed because of too many classes and too few observations")
		}
	}
	else t <- matrix(1,N,1)
	
	likely <- c()
	I <- 0
	test <- Inf
	while ((I <- I+1)<=itermax && test>=eps){
		if (algo!='EM' && I!=1) t <- t2
		if (K>1 && (any(is.na(t)) || any(colSums(t>1/K)<min.individuals))) return(1)
		m <- hddc_m_step(DATA,K,t,model,threshold,method,noise.ctrl,com_dim)
		t <- hddc_e_step(DATA,m)
		L <- t$L
		t <- t$t
		if (algo=='CEM') {
			t2 <- matrix(0,N,K)
			t2[cbind(1:N,max.col(t))] <- 1
		}
		else if(algo=='SEM') { 
			t2 <- matrix(0,N,K)
			for (i in 1:N)	t2[i,] <- t(rmultinom(1,1,t[i,]))
		}
		likely[I] <- L
		if (I!=1) test <- abs(likely[I]-likely[I-1])
	}
	
	if ( model%in%c('AKBKQKDK','AKBQKDK','AKBKQKD','AKBQKD') ) {
		a <- matrix(m$a[,1],1,m$K,dimnames=list(c("Ak:"),1:m$K))
	}
	else if(model=='AJBQD') {
		a <- matrix(m$a[1,],1,m$d[1],dimnames=list(c('Aj:'),paste('a',1:m$d[1],sep='')))
	}
	else if ( model%in%c('ABKQKDK','ABQKDK','ABKQKD','ABQKD',"ABQD") ) {
		a <- matrix(m$a[1],dimnames=list(c('A:'),c('')))
	}
	else a <- matrix(m$a,m$K,max(m$d),dimnames=list('Class'=1:m$K,paste('a',1:max(m$d),sep='')))
	
	if ( model%in%c('AKJBQKDK','AKBQKDK','ABQKDK','AKJBQKD','AKBQKD','ABQKD','AJBQD',"ABQD") ) {
		b <- matrix(m$b[1],dimnames=list(c('B:'),c('')))
	}
	else b <- matrix(m$b,1,m$K,dimnames=list(c("Bk:"),1:m$K))
	
	d <- matrix(m$d,1,m$K,dimnames=list(c('dim:'),"Intrinsic dimensions of the classes:"=1:m$K))
	mu <- matrix(m$mu,m$K,p,dimnames=list('Class'=1:m$K,'Posterior group means:'=paste('V',1:p,sep='')))
	prop <- matrix(m$prop,1,m$K,dimnames=list(c(''),'Posterior probabilities of groups'=1:m$K))
	class(b) <- class(a) <- class(d) <- class(prop) <- class(mu) <- class(t) <- 'hd'
	cls <- max.col(t)
	list(model=model,K=K,d=d,a=a,b=b,mu=mu,prop=prop,ev=m$ev,Q=m$Q,loglik=likely,posterior=t,class=cls,com_ev=com_ev,N=N)
}

hddc_e_step  <- function(x,par){
	p <- ncol(x)
	N <- nrow(x)
	K <- par$K
	a <- par$a
	b <- par$b
	mu <- par$mu
	d <- par$d
	prop <- par$prop
	Q <- par$Q
	
	b[b<1e-6] <- 1e-6

	if(par$model=="AJBQD") {
		K_pen <- diag((mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(x))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)+2*(mu%*%Q)%*%(t(Q)%*%t(x))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
	}
	else if(par$model=="ABQD") {
		K_pen <- diag(1/a[1]*(mu%*%Q)%*%(t(Q)%*%t(mu)))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))+2*(1/b[1]-1/a[1])*(mu%*%Q)%*%(t(Q)%*%t(x))
	}
	else{
		K_pen <- matrix(0,K,N)
		for (i in 1:K) {
			s <- sum(log(a[i,1:d[i]]))
			X <- x-matrix(mu[i,],N,p,byrow=TRUE)
			proj <- (X%*%Q[[i]])%*%t(Q[[i]])
			A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
			B <- X-proj
			K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
		}
	}
	A <- -1/2*t(K_pen)
	L <- sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
	
	t <- matrix(0,N,K)
	for (i in 1:K) t[,i] <- 1/rowSums(exp((K_pen[i,]-t(K_pen))/2))
	list(t=t,L=L)
}

hddc_m_step  <- function(x,K,t,model,threshold,method,noise.ctrl,com_dim){
	N <- nrow(x)
	p <- ncol(x)
	prop <- c()
	n <- colSums(t)
	prop <- n/N
	mu <- matrix(,K,p)
	for (i in 1:K) mu[i,] <- colSums(x*t[,i])/n[i]
	
	ind <- apply(t>0,2,which)
	n_bis <- c()
	for(i in 1:K) n_bis[i] <- length(ind[[i]])
	
	#Calculation on Var/Covar matrices
	
	if (N<p) {
		if( model%in%c("AJBQD","ABQD") ){
			Y <- matrix(0,N,p)
			for (i in 1:K) Y <- Y+(x-matrix(mu[i,],N,p,byrow=TRUE))/sqrt(N)*sqrt(t[,i])
			donnees <- eigen(tcrossprod(Y),symmetric=TRUE)
			ev <- donnees$values
		}
		else{
			Y <- vector(mode='list',length=K)
			ev <- matrix(0,K,N)
			Q <- vector(mode='list',length=K)
			for (i in 1:K){ 
				Y[[i]] <- (x-matrix(mu[i,],N,p,byrow=TRUE))/sqrt(n[i])*sqrt(t[,i])
				donnees <- eigen(tcrossprod(Y[[i]]),symmetric=TRUE)
				ev[i,1:N] <- donnees$values
				Q[[i]] <- donnees$vectors
			}
		}
	}
	else if ( model%in%c("AJBQD","ABQD") ){
		W <- matrix(0,p,p)
		for (i in 1:K) W <- W + crossprod((x-matrix(mu[i,],N,p,byrow=TRUE))*sqrt(t[,i]))/N
		donnees <- eigen(W,symmetric=TRUE)
		ev <- donnees$values
	}
	else {
		ev <- matrix(0,K,p)
		Q <- vector(mode='list',length=K)
		for (i in 1:K){ 
			donnees <- eigen(crossprod((x-matrix(mu[i,],N,p,byrow=TRUE))*sqrt(t[,i]))/n[i],symmetric=TRUE)
			ev[i,] <- donnees$values
			Q[[i]] <- donnees$vectors
		}
	}	
	
	#Intrinsic dimensions selection
	
	if (model%in%c("AJBQD","ABQD")) d <- rep(com_dim,length=K)
	else if ( model%in%c("AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD") ){
		dmax <- min(apply((ev>noise.ctrl)*rep(1:ncol(ev),each=K),1,which.max))-1
		if(com_dim>dmax) com_dim <- max(dmax,1)
		d <- rep(com_dim,length=K)
	}
	else d <- hdclassif_dim_choice(ev,n,method,threshold,FALSE,noise.ctrl)
	
	#Setup of the Qi matrices	
	
	if ( model%in%c("AJBQD","ABQD") ){
		if (N>=p) Q <- matrix(donnees$vectors[,1:d[1]],p,d[1])
		else {
			Q <- matrix(t(Y)%*%donnees$vectors[,1:d[1]],p,d[1])
			normalise <- c()
			for(i in 1:d[1]) normalise[i] <- as.double(crossprod(Q[,i]))
			Q <- Q/matrix(sqrt(normalise),p,d,byrow=TRUE)
		}
	}
	else if (N>=p) for(i in 1:K) Q[[i]] <- matrix(Q[[i]][,1:d[i]],p,d[i])
	else{
		for (i in 1:K){ 
			Q[[i]] <- t(Y[[i]])%*%(Q[[i]][,1:d[i]])
			normalise <- c()
			for (j in 1:d[i]) normalise[j] <- as.double(crossprod(as.matrix(Q[[i]][,j])))
			Q[[i]] <- Q[[i]]/matrix(sqrt(normalise),p,d[i],byrow=TRUE)
		}
	}
	
	#Calculation of the remaining parameters of the selected model	
	
	ai <- matrix(NA,K,max(d))
	if ( model%in%c('AKJBKQKDK','AKJBQKDK','AKJBKQKD','AKJBQKD') ){
		for (i in 1:K) ai[i,1:d[i]] <- ev[i,1:d[i]]
	}
	else if ( model%in%c('AKBKQKDK','AKBQKDK' ,'AKBKQKD','AKBQKD') ){
		for (i in 1:K) ai[i,] <- rep(sum(ev[i,1:d[i]])/d[i],length=max(d))
	}
	else if (model=="AJBQD") for (i in 1:K) ai[i,] <- ev[1:d[1]]
	else if (model=="ABQD")	ai[] <- sum(ev[1:d[1]])/d[1]
	else {
		a <- 0
		eps <- sum(prop*d)
		for (i in 1:K) a <- a + sum(ev[i,1:d[i]])*prop[i]
		ai <- matrix(a/eps,K,max(d))
	}

	bi <- c()
	if ( model%in%c('AKJBKQKDK','AKBKQKDK','ABKQKDK','AKJBKQKD','AKBKQKD','ABKQKD') ){
		for(i in 1:K) bi[i] <- sum(ev[i,(d[i]+1):min(N,p)])/(p-d[i])
	}
	else if ( model%in%c("ABQD","AJBQD") ){
		bi[1:K] <- sum(ev[(d[1]+1):min(N,p)])/(min(N,p)-d[1])
	}
	else{		
		b <- 0
		eps <- sum(prop*d)
		for(i in 1:K) b <- b + sum(ev[i,(d[i]+1):min(N,p)])*prop[i]
		bi[1:K] <- b/(min(N,p)-eps)
	}

	list(model=model,K=K,d=d,a=ai,b=bi,mu=mu,prop=prop,ev=ev,Q=Q)
}

plot.hdc  <- function(x,method=NULL,threshold=NULL,noise.ctrl=1e-8,...){
	method <- if(!is.null(method)) method else if(!is.null(x$threshold)) "C" else "B"
	method <- toupper(method)
	method <- if (method%in%c("CATTELL","C")) "C" else if (method%in%c("BIC","B")) "B" 
	threshold <- if(!is.null(threshold)) threshold else if(!is.null(x$threshold)) x$threshold else  0.2
	
	if(!method%in%c("C","B")) stop("Wrong method name.\n",call.=FALSE)
	
	k <- x$K
	N <- x$N
	n <- x$prop*N
	
	if(is.null(x$com_ev)) d <- hdclassif_dim_choice(x$ev,n,method,threshold,TRUE,noise.ctrl)
	else d <- hdclassif_dim_choice(x$com_ev,n,method,threshold,TRUE,noise.ctrl)
}

predict.hdc  <- function(object,data,cls=NULL,...){
	#Extract variables:
	p <- ncol(data)
	N <- nrow(data)
	K <- object$K
	a <- object$a
	b <- object$b
	mu <- object$mu
	d <- object$d
	prop <- object$prop
	Q <- object$Q
	x <- as.matrix(data)
	confusion <- NULL
	ARI <- NULL
	if (length(N)==0) {
		N <- 1
		p <- length(data)
		x <- matrix(data,N,p)
	}
	if (length(object$scaling)!=0){
		x <- scale(x,center=object$scaling$mu,scale=object$scaling$sd)
	}
	
	if(length(b)==1) b <- rep(b,length=K)
	if (length(a)==1) a <- matrix(a,K,max(d))
	else if (length(a)==K) a <- matrix(a,K,max(d))
	else if (object$model=='AJBQD') a <- matrix(a,K,d[1],byrow=TRUE)
	
	if(min(a,na.rm=TRUE)<=0 | min(b)<=0) stop("Some parameters A or B are negative. Prediction can't be done.\nThe reduction of the intrinsic dimensions or a more constrained model can be a solution.\nAlso you can change the value of A's and B's manually by accessing the paramaters.\n",call.=FALSE)
	
	
	if(object$model=="AJBQD") {
		K_pen <- diag((mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(x))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)+2*(mu%*%Q)%*%(t(Q)%*%t(x))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
	}
	else if (object$model=="ABQD") {
		K_pen <- diag(1/a[1]*(mu%*%Q)%*%(t(Q)%*%t(mu)))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))+2*(1/b[1]-1/a[1])*(mu%*%Q)%*%(t(Q)%*%t(x))
	}
	else{
		K_pen <- matrix(0,K,N)
		for (i in 1:K) {
			s <- sum(log(a[i,1:d[i]]))
			X <- x-matrix(mu[i,],N,p,byrow=TRUE)
			proj <- (X%*%Q[[i]])%*%t(Q[[i]])
			A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
			B <- X-proj
			K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])
		}
	}
	
	t <- matrix(0,N,K,dimnames=list(1:N,1:K))
	for (i in 1:K) t[,i] <- 1/rowSums(exp((K_pen[i,]-t(K_pen))/2))
	result <- max.col(t)
	
	if (!is.null(object$kname)){
		result <- factor(result,labels=object$kname,levels=seq_along(object$kname))
		colnames(t) <- object$kname
	}
	if (!is.null(cls)){
		if(is.null(object$kname)){
			ARI <- hddc_ari(cls,result)
			cat("Adjusted Rand Index: ",ARI,".\n")
		}
		else{
			confusion <- table(result,cls)
			dimnames(confusion) <- list('Predicted class'=object$kname,'Initial class'=object$kname)
			cat("Correct classification rate: ",sum(diag(confusion))/N,".\n",sep="")
			print(confusion)
		}
	}
	class(t) <- 'hd'
	list(class=result,posterior=t,confusion=confusion,ARI=ARI)
}

print.hd <- function(x,...){
	class(x) <- NULL
	print.default(x,digits=3,na.print='.')
	class(x) <- 'hd'
}

print.hdc <- function(x,...){
	if(length(x$kname)!=0) cat ("HIGH DIMENSIONAL DISCRIMINANT ANALYSIS\nMODEL: ",x$model,"\n",sep='')
	else cat ("HIGH DIMENSIONAL DATA CLUSTERING\nMODEL: ",x$model,"\n",sep='')
	print(x$prop)
	print(x$d)
	print(x$a)
	print(x$b)
	cat("BIC: ",x$BIC,"\n")
	if(!is.null(x$info)) cat(x$info,"\n")
	if(min(x$a,na.rm=TRUE)<0) cat("Information: a < 0\n")
	if(min(x$b)<10e-6) cat("Information: b < 10e-6\n")
}

simuldata <- function(nlearn,ntest,p,K=3,prop=NULL,d=NULL,a=NULL,b=NULL){
	N=nlearn+ntest
	if (length(prop)==0) prop<-rep(1/K,K)
	else if (length(prop)!=K) stop("Proportions don't fit with the number of classes.")
	else prop<-prop/sum(prop)
	
	# Class sizes
	n<-floor(prop*N)
	N<-sum(n)

	#MEANS
	mu<-matrix(0,K,p)
	j<-sample(p,K)
	mu[cbind(1:K,j)]<-10
	
	# Intrinsic dimensions
	if ( length(d)==0 )	d<-sort(ceiling(runif(K,0,12*(p>20)+5*(p<=20 && p>=6)+(p<6)*(p-1))),decreasing=TRUE)
	else if ( length(d)!=K || !any(is.numeric(d)) ) stop("Wrong value of d.")
	
	# Orientation matrices
	Q<-vector(mode='list',length=K)
	for (i in 1:K) Q[[i]]<-qr.Q(qr(mvrnorm(p,mu=rep(0,p),Sigma=diag(1,p))))
	
	# Variance in the class-specific subspace
	if ( length(a)==0 ) a<-sort(ceiling(runif(K,30,350)))
	else if ( length(a)!=K || !any(is.numeric(a)) ) stop("Wrong value of a.")
	if ( length(b)==0 )b<-sort(ceiling(runif(K,0,25)))
	else if ( length(b)!=K || !any(is.numeric(b)) ) stop("Wrong value of b.")
	
	# Simulation
	S<-vector(mode='list',length=K)
	for (i in 1:K)	S[[i]]<-crossprod(Q[[i]]%*%sqrt(diag(c(rep(a[i],d[i]),rep(b[i],p-d[i])))))
	
	cls<-X<-NULL
	for (i in 1:K)	X<-rbind(X,mvrnorm(n[i],mu=mu[i,],Sigma=S[[i]]))
	
	for (i in 1:K) cls<-c(cls,rep(i,n[i]))
	
	ind<-sample(1:N,N)
	prms<-list(a=a,b=b,prop=prop,d=d,mu=mu)
	data <- list(X=X[ind[1:nlearn],],clx=cls[ind[1:nlearn]],Y=X[ind[(nlearn+1):N],],cly=cls[ind[(nlearn+1):N]],prms=prms)
	
}

hddc_ari <- function(x,y){
	#This function is drawn from the mclust package
	x <- as.vector(x)
	y <- as.vector(y)
	tab <- table(x, y)
	if (all(dim(tab) == c(1, 1))) return(1)
	a <- sum(choose(tab, 2))
	b <- sum(choose(rowSums(tab), 2)) - a
	c <- sum(choose(colSums(tab), 2)) - a
	d <- choose(sum(tab), 2) - a - b - c
	ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
	return(ARI)
}

