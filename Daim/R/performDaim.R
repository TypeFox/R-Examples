

performDaim <- function(x, ...) UseMethod("performDaim")


performDaim.default <- function(x, ...) {
  stop(paste("Do not know how to handle objects of class", class(x)))
}


performDaim.matrix <- function(x, labels, prediction=NULL,
				thres=seq(0,1,by=0.01), cutoff=0.5, labpos="1", ...)
{
	Xdim <- dim(x)
	est.method <- "obs"
	if(!is.null(names(list(...)))){
		if(names(list(...)) == "est.method")
			est.method <- "prob.mean"
	}
	if(Xdim[1] != length(labels))
		stop("\n Number of labels must be equal to the number of observations in x.\n")
	if(Xdim[2] < 1)
		stop("\n Number of runs must be > 0.\n")
	if(is.character(labels))
		labels <- factor(labels)
	rsp.labels <- levels(labels)
	id <- as.logical(match(rsp.labels,labpos,nomatch=0))
	if(sum(id) < 1)
		stop("\n wrong level specified !\n")
	rsp.labels[id] <- "pos"
	rsp.labels[!id] <- "neg"
	levels(labels) <- rsp.labels
	labels <- factor(as.character(labels))
	labpos <- 1
	labneg <- 0
	method <- list(est.method=est.method, cutoff=cutoff)
  	if(est.method == "obs")
		est.method <- 0
	else if(est.method == "prob.mean")
		est.method <- 1
	else
		stop("\n est.method must be one of 'obs' or 'prob.mean'! \n")
	
	if(sum(is.na(x)) > 0){
		testind <- as.list(apply(x,2,function(x) which(!is.na(x))))
		lab.oob <- lapply(testind,function(x,y) y[x],y=as.numeric(labels) - 1)
		prob.oob <- as.list(apply(x,2,function(x) x[complete.cases(x)]))
	}
	else{
		x <- as.data.frame(x)
		testind <- rep(list(1:nrow(x)),ncol(x))
		lab.oob <- lapply(testind,function(x,y) y[x],y=as.numeric(labels) - 1)
		prob.oob <- as.list(x)
	}
	ans <- Daim.performance(testind, lab.oob, prob.oob, labels, prediction,
				est.method, thres, cutoff, labneg, method)
	ans
}



performDaim.data.frame <- function(x, labels, prediction=NULL,
				thres=seq(0,1,by=0.01), cutoff=0.5, labpos="1", ...)
{
	Xdim <- dim(x)
	est.method <- "obs"
	if(!is.null(names(list(...)))){
		if(names(list(...)) == "est.method")
			est.method <- "prob.mean"
	}
	if(Xdim[1] != length(labels))
		stop("\n Number of labels must be equal to the number of observations in x.\n")
	if(Xdim[2] < 1)
		stop("\n Number of runs must be > 0.\n")
	if(is.character(labels))
		labels <- factor(labels)
	rsp.labels <- levels(labels)
	id <- as.logical(match(rsp.labels,labpos,nomatch=0))
	if(sum(id) < 1)
		stop("\n wrong level specified !\n")
	rsp.labels[id] <- "1"
	rsp.labels[!id] <- "0"
	levels(labels) <- rsp.labels
	labels <- factor(as.character(labels))
	labpos <- 1
	labneg <- 0
	method <- list(est.method=est.method, cutoff=cutoff)
  	if(est.method == "obs")
		est.method <- 0
	else if(est.method == "prob.mean")
		est.method <- 1
	else
		stop("\n est.method must be one of 'obs' or 'prob.mean'! \n")	
	if(sum(is.na(x)) > 0){
		testind <- as.list(apply(x,2,function(x) which(!is.na(x))))
		lab.oob <- lapply(testind,function(x,y) y[x],y=as.numeric(labels) - 1)
		prob.oob <- as.list(apply(x,2,function(x) x[complete.cases(x)]))
	}
	else{
		testind <- rep(list(1:nrow(x)),ncol(x))
		lab.oob <- lapply(testind,function(x,y) y[x],y=as.numeric(labels) - 1)
		prob.oob <- as.list(x)
	}
	ans <- Daim.performance(testind, lab.oob, prob.oob, labels, prediction, 
				est.method, thres, cutoff, labneg, method)
	ans
}








Daim.performance <- function(testind, lab.oob, prob.oob, labels, prob.app,
				est.method, thres, cutoff, labneg, method)
{
  	N.thres <- length(thres)
  	all.roc  <-  .Call("roc_value", 
						prob.oob, lab.oob, 
						thres, as.numeric(labneg),
						PACKAGE="Daim")

  	testind <- unlist(testind,use.names = FALSE)
  	prob.oob <- unlist(prob.oob,use.names = FALSE)
	loob <- split(prob.oob,testind)
	ind.rsp <- unique(testind)
	if(length(ind.rsp) < length(labels)){
		warning(paste("\n Only ", length(ind.rsp),
			"observations  have a probability records.\n
			 Increase runs and try again!"))
		labels <- labels[ind.rsp]
	}
	if(is.null(prob.app))
		prob.app <- rep(0,length(labels))
	if(length(prob.app) != length(loob))
		stop("\n Number of predictions and observations must be equal \n")
	if(!est.method){
		if(is.character(cutoff))
			ans <- .Call("sens_spez_obs_cut",loob, 
				as.numeric(prob.app),
				as.numeric(thres),
				as.numeric(labneg),
				as.numeric(0),
				as.numeric(labels),
				PACKAGE="Daim")
		else
			ans <- .Call("sens_spez_obs",loob, 
				as.numeric(prob.app),
				as.numeric(thres),
				as.numeric(labneg),
				as.numeric(cutoff),
				as.numeric(labels),
				PACKAGE="Daim")	
	}
	else{
		ans <- .Call("sens_spez_none",loob, 
				as.numeric(prob.app),
				as.numeric(thres),
				as.numeric(labneg),
				as.numeric(cutoff),
				as.numeric(labels),
				PACKAGE="Daim")
	}
	sensloob <- ans$sensloob
	spezloob <- ans$spezloob
	sensapp <- ans$sensapp
	spezapp <- ans$spezapp
	myq <- ans$myq
	myp <- ans$myp

	if(!is.character(cutoff)){
		id.cut <- which(abs(thres - cutoff) < .Machine$double.eps)
		errapp <- ans$errapp
		errloob <- ans$errloob
	}
	else{
		errapp <- ans$errapp2
		errloob <- ans$errloob2
	}

  	ffrapp <- 1-sensapp
  	fprapp <- 1-spezapp
  	ffrloob <- 1-sensloob
  	fprloob <- 1-spezloob

  	#### .632
  	ffr632 <- 0.368*ffrapp + 0.632*ffrloob
  	fpr632 <- 0.368*fprapp + 0.632*fprloob
  	sens632 <- 1 - ffr632
  	spez632 <- 1 - fpr632
  	err632 <- 0.368*errapp + 0.632*errloob

 	#### .632+
  	fpr632p <- ffr632p <- vector(mode="numeric",N.thres)
  	
  	ffrloob1 <- pmin(ffrloob,(1-myq))
  	fprloob1 <- pmin(fprloob,myq)
  	Rffr <- (ffrloob1-ffrapp)/((1-myq)-ffrapp)
  	Rfpr <- (fprloob1-fprapp)/(myq-fprapp)
  	Rffr[ffrloob1 <= ffrapp] <- 0
  	Rfpr[fprloob1 <= fprapp] <- 0
  	ffr632p <- ffr632 + (ffrloob1-ffrapp) * (0.368*0.632*Rffr)/(1-0.368*Rffr)
  	fpr632p <- fpr632 + (fprloob1-fprapp) * (0.368*0.632*Rfpr)/(1-0.368*Rfpr)
  	sens632p <- 1 - ffr632p
  	spez632p <- 1 - fpr632p
	
	gamm <- myp*(1-myq) + (1-myp)*myq
	if(!is.character(cutoff)){
		errloob1 <- pmin(errloob, gamm[id.cut])
	}else{
		errloob1 <- pmin(errloob, gamm)
	}
	
  	
	Rerr <- rep(0,length(errapp))
	for(i in 1:length(Rerr)){
		if(!(errloob[i] <= errapp[i] | gamm[i] <= errapp[i]))
			Rerr[i] <- (errloob[i]-errapp[i])/(gamm[i]-errapp[i])
	}
	
  	err632p <- err632 + (errloob1 - errapp) * (0.368*0.632*Rerr)/(1 - 0.368*Rerr)

	if(cutoff == "loob"){
		id.cut <- which.max(sensloob+spezloob-1)
		control$best.cutoff <- thres[id.cut]
	}	
	if(cutoff == "0.632"){
		id.cut <- which.max(sens632+spez632-1)
		method$best.cutoff <- thres[id.cut]
	}
	if(cutoff == "0.632+"){
		id.cut <- which.max(sens632p+spez632p-1)
		mehod$best.cutoff <- thres[id.cut]
	}
	if(is.character(cutoff)){
		err632p <- err632p[id.cut]
		err632 <- err632[id.cut]
		errloob <- errloob[id.cut]
		errapp <- errapp[id.cut]
	}

  	output <- list(method=method, err632p=err632p,err632=err632,
					errloob=errloob,errapp=errapp,
					sens632p=sens632p[id.cut], spec632p=spez632p[id.cut],
					sens632=sens632[id.cut], spec632=spez632[id.cut],
					sensloob=sensloob[id.cut], specloob=spezloob[id.cut],
					sensapp=sensapp[id.cut], specapp=spezapp[id.cut],
  					roc=data.frame(sens632p=sens632p,spec632p=spez632p,
								sens632=sens632,spec632=spez632,
								sensloob=sensloob,specloob=spezloob,
								sensapp=sensapp,specapp=spezapp,cut.points=thres),
  					sample.roc=all.roc)
  	class(output) <- c("Daim","predictions")
  	output
}

 
