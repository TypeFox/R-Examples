

############################################################################
###       function to evaluate the accuracy of classification models
###       with binary outcome
############################################################################


Daim <- function(formula, model = NULL, data = NULL, control = Daim.control(),			
			thres = seq(0, 1, by=0.01), cutoff = 0.5, labpos = "1", returnSample = FALSE, 
			cluster = NULL, seed.cluster = NULL, multicore = FALSE, ...)
{
	
	call <- match.call()
	form <- terms(formula,data=data)
	data <- get_all_vars(formula,data)
	est.method <- "obs"
	if(!is.null(names(list(...)))){
		if(names(list(...)) == "est.method")
		est.method <- "prob.mean"
	}
	if(is.character(cutoff)){
		def.cutoff <- c("0.632", "0.632+", "cv", "loob")
		cutoff <- charmatch(cutoff, def.cutoff)
		if(is.na(cutoff))
			stop("\n 'cutoff' must be one of 'loob', '0.632', '0.632+' or 'cv'! \n")
		cutoff <- def.cutoff[cutoff]
	}
	control$cutoff <- cutoff
	control$est.method <- est.method
	method <- control$method
	if(method == "cv" & is.character(cutoff) & cutoff != "cv")
		stop("\n 'cutoff' must be 'cv'! \n")
	if(method != "cv" & is.character(cutoff) & cutoff == "cv")
		stop("\n 'cutoff' must be 'loob', '0.632' or '0.632+'! \n")
	rsp.labels <- levels(data[[1]])
	id <- as.logical(match(rsp.labels, labpos, nomatch=0))
	if(sum(id) < 1)
		stop("\n wrong level specified !\n")
	rsp.labels[id] <- "pos"
	rsp.labels[!id] <- "neg"
	levels(data[[1]]) <- rsp.labels
	data[[1]] <- factor(as.character(data[[1]]))
	labpos <- 1
	labneg <- 0
	N <- nrow(data)
	boot.meth <- charmatch(method,c("boot", "cv", "bcv"))
  	if(est.method == "obs")
		est.meth <- 0
	else if(est.method == "prob.mean")
		est.meth <- 1
	else
		stop("\n 'est.method' must be one of 'obs' or 'prob.mean'! \n")
	if(!is.null(control$dependency)){
		dep.ind <- Daim.depend(data, control$dependency)
		data <- dep.ind[[1]]
		N <- dep.ind[[2]]
		dep.ind <- dep.ind[[3]]
	}
	else{
		dep.ind <- NULL
	}
	if(is.null(cluster) && !multicore){
		if(boot.meth == 1){
			nboot <- control$nboot
			replace <- control$replace
			if(!replace){
				boot.size <- control$boot.size
			}else{
				boot.size <- 1
			}
			prob.oob <- lab.oob <- testind <- inbaggind <- vector(mode="list",length=nboot)
			for(i in 1:nboot){
				mylist <- Daim.boot.index(N, boot.size, replace, dep.ind)
				train <- data[mylist,]
				testid <- unique(mylist)
				test <- data[-testid,]
				inbaggind[[i]] <- mylist
				testind[[i]] <- (1:nrow(data))[-testid]
				prob.oob[[i]] <- model(formula, train,test)
				lab.oob[[i]] <- as.numeric(test[[1]]) - 1
			}			
		}
		if(boot.meth == 2){
			xval <- control$k
			xval.runs <- control$k.runs
			if(xval > N){
				xval <- N
				control$k <- xval
			}
			prob.oob <- lab.oob <- testind <- vector(mode="list", length=xval.runs*xval)
			k <- 1
			for(j in 1:xval.runs){
				xgr <- 1:xval
				id <- sample(rep(xgr, length = N), N)
				for(i in xgr){
					test.id <- id == i
					if(!is.null(dep.ind)){
						train <- data[unlist(dep.ind[!test.id], use.names=FALSE),]
						test.id <- unlist(dep.ind[test.id], use.names=FALSE)
						test <- data[test.id,]
						testind[[k]] <- test.id
					}
					else{
						train <- data[!test.id,]
						test <- data[test.id,]
						testind[[k]] <- which(test.id)
					}
					prob.oob[[k]] <- model(formula, train, test)
					lab.oob[[k]] <- as.numeric(test[[1]]) - 1
					k <- k+1
				}
			}
		}
		if(boot.meth == 3){
			xval <- control$k
			nboot <- control$nboot
			replace <- control$replace
			if(!replace){
				boot.size <- control$boot.size
			}else{
				boot.size <- 1
			}
			prob.oob <- lab.oob <- testind <- vector(mode="list",length=nboot)
			for(i in 1:nboot){
				ID.B <- sample(1:N, N*boot.size, replace=replace)
				ID.CV <- split(ID.B, ID.B)
				xgr <- 1:xval
				id <- sample(rep(xgr, length = length(ID.CV)), length(ID.CV))
				prob.oob.cv <- id.oob.cv <- lab.oob.cv <- NULL
				for(j in xgr){
					test.id <- id == j
					train <- data[unlist(ID.CV[!test.id], use.names = FALSE),]
					test <- data[unlist(ID.CV[test.id], use.names = FALSE),]
					prob.oob.cv <- c(prob.oob.cv,model(formula, train, test))
					id.oob.cv <- c(id.oob.cv, unlist(ID.CV[test.id], use.names = FALSE))
					lab.oob.cv <- c(lab.oob.cv, as.numeric(test[[1]])-1)
				}
				testind[[i]] <- id.oob.cv
				prob.oob[[i]] <- prob.oob.cv
				lab.oob[[i]] <- lab.oob.cv
			}
		}
	}
	else { 
		if(multicore){
			if(boot.meth == 1){
				nboot <- control$nboot
				replace <- control$replace
				if(!replace){
					boot.size <- control$boot.size
				}else{
					boot.size <- 1
				}
				if(!is.null(seed.cluster)){
					RNGkind("L'Ecuyer-CMRG")
					set.seed(seed.cluster)
					BG <- mclapply(1:nboot, Daim.cluster.boot,
								   formula=formula, model=model, data=data, N=N,
								   boot.size=boot.size, replace=replace, 
								   mc.set.seed=TRUE, ...)
				}else{
					BG <- mclapply(1:nboot, Daim.cluster.boot,
								   formula=formula, model=model, data=data, N=N,
								   boot.size=boot.size, replace=replace, ...)
				}
				testind <- lapply(BG,function(x) x$testind)
				inbaggind <- lapply(BG,function(x) x$inbaggind)
				prob.oob <- lapply(BG,function(x) x$prob.oob)	
				lab.oob <- lapply(BG,function(x) x$lab.oob)
			}
			if(boot.meth == 2){
				xval <- control$k
				xval.runs <- control$k.runs
				if(xval > N){
					xval <- N
					control$k <- xval
				}
				if(xval.runs > 1){
					BG <- mclapply(1:xval.runs, Daim.cluster.cv0,
								   formula=formula, model=model,
								   data=data, N=N, xval=xval)					
					testind <- unlist(lapply(BG,function(x) x$testind),FALSE)
					prob.oob <- unlist(lapply(BG,function(x) x$prob.oob),FALSE)	
					lab.oob <- unlist(lapply(BG,function(x) x$lab.oob),FALSE)
				}
				if(xval.runs == 1){
					xgr <- 1:xval
					id <- sample(rep(xgr, length = N), N)
					nboot <- length(xgr)
					BG <- mclapply(xgr, Daim.cluster.cv,
								   formula=formula, model=model,
								   data=data, N=N, xval=xval, id=id)
					testind <- lapply(BG,function(x) x$testind)
					prob.oob <- lapply(BG,function(x) x$prob.oob)	
					lab.oob <- lapply(BG,function(x) x$lab.oob)
				}
			}
			if(boot.meth == 3){
				xval <- control$k
				nboot <- control$nboot
				replace <- control$replace
				if(!replace){
					boot.size <- control$boot.size
				}else{
					boot.size <- 1
				}
				BG <- mclapply(1:nboot, Daim.cluster.bcv,
							   formula=formula, model=model, data=data, N=N,
							   boot.size=boot.size, replace=replace, xval=xval)	
				testind <- lapply(BG,function(x) x$testind)
				prob.oob <- lapply(BG,function(x) x$prob.oob)	
				lab.oob <- lapply(BG,function(x) x$lab.oob)			
			}
		}
		else{
			if(is.null(seed.cluster)){
				clusterSetRNGStream(cluster, iseed=sample(1:9999,1))
			}
			else{
				clusterSetRNGStream(cluster, iseed=seed.cluster)
			}
			clusterEvalQ(cluster, library(Daim))
			if(boot.meth == 1){
				nboot <- control$nboot
				replace <- control$replace
				if(!replace){
					boot.size <- control$boot.size
				}else{
					boot.size <- 1
				}
				BG <- clusterApplyLB(cluster, 1:nboot, Daim.cluster.boot,
									 formula=formula, model=model, data=data, N=N,
									 boot.size=boot.size, replace=replace, ...)
				testind <- lapply(BG,function(x) x$testind)
				inbaggind <- lapply(BG,function(x) x$inbaggind)
				prob.oob <- lapply(BG,function(x) x$prob.oob)	
				lab.oob <- lapply(BG,function(x) x$lab.oob)
			}
			if(boot.meth == 2){
				xval <- control$k
				xval.runs <- control$k.runs
				if(xval > N){
					xval <- N
					control$k <- xval
				}
				if(xval.runs > 1){
					BG <- clusterApplyLB(cluster, 1:xval.runs, Daim.cluster.cv0,
										 formula=formula, model=model, data=data, 
										 N=N, xval=xval)
					testind <- unlist(lapply(BG,function(x) x$testind),FALSE)
					prob.oob <- unlist(lapply(BG,function(x) x$prob.oob),FALSE)	
					lab.oob <- unlist(lapply(BG,function(x) x$lab.oob),FALSE)
				}
				if(xval.runs == 1){
					xgr <- 1:xval
					id <- sample(rep(xgr, length = N), N)
					nboot <- length(xgr)
					BG <- clusterApplyLB(cluster, xgr, Daim.cluster.cv,
										 formula=formula, model=model, 
										 data=data, N=N, xval=xval, id=id)
					testind <- lapply(BG,function(x) x$testind)
					prob.oob <- lapply(BG,function(x) x$prob.oob)	
					lab.oob <- lapply(BG,function(x) x$lab.oob)
				}
			}
			if(boot.meth == 3){
				xval <- control$k
				nboot <- control$nboot
				replace <- control$replace
				if(!replace){
					boot.size <- control$boot.size
				}else{
					boot.size <- 1
				}
				BG <- clusterApplyLB(cluster, 1:nboot, Daim.cluster.bcv,
									 formula=formula, model=model, data=data, N=N,
									 boot.size=boot.size, replace=replace, xval=xval)
				testind <- lapply(BG,function(x) x$testind)
				prob.oob <- lapply(BG,function(x) x$prob.oob)	
				lab.oob <- lapply(BG,function(x) x$lab.oob)
			}
		}
	}
	prob.app <- model(formula,data,data)
	
	rsp <- data[[1]]
  	N.thres <- length(thres)
  	all.roc  <-  .Call("roc_value", 
					   prob.oob, lab.oob, 
					   thres, as.numeric(labneg),
					   PACKAGE="Daim")
	if(returnSample){
		if(boot.meth == 1){
			all.data <- list(testind=testind, inbaggind=inbaggind,
							 prob.oob=prob.oob, lab.oob=lab.oob, prob.app=prob.app)
		}else{
			all.data <- list(testind=testind, prob.oob=prob.oob, 
							 lab.oob=lab.oob, prob.app=prob.app)
		}
	}else{
		all.data <- NULL
	}
  	testind <- unlist(testind, use.names = FALSE)
  	prob.oob <- unlist(prob.oob, use.names = FALSE)
	loob <- split(prob.oob, testind)
	
	if(!est.meth){
		if(is.character(cutoff)){
			ans <- .Call("sens_spez_obs_cut",loob, 
					 as.numeric(prob.app),
					 as.numeric(thres),
					 as.numeric(labneg),
					 as.numeric(0),
					 as.numeric(rsp),
					 PACKAGE="Daim")
		}else{
			ans <- .Call("sens_spez_obs",loob, 
					 as.numeric(prob.app),
					 as.numeric(thres),
					 as.numeric(labneg),
					 as.numeric(cutoff),
					 as.numeric(rsp),
					 PACKAGE="Daim")
		}
	}
	else{
		ans <- .Call("sens_spez_none",loob, 
					 as.numeric(prob.app),
					 as.numeric(thres),
					 as.numeric(labneg),
					 as.numeric(cutoff),
					 as.numeric(rsp),
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
	if(boot.meth != 2){
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
		ffrloob1 <- pmin(ffrloob, (1-myq))
		fprloob1 <- pmin(fprloob, myq)
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
		
		Rerr <- vector("numeric",length(errapp))
		for(i in 1:length(Rerr)){
			if(!(errloob[i] <= errapp[i] | gamm[i] <= errapp[i]))
				Rerr[i] <- (errloob[i]-errapp[i])/(gamm[i]-errapp[i])
		}
		err632p <- err632 + (errloob1-errapp) * (0.368*0.632*Rerr)/(1-0.368*Rerr)
		
		if(cutoff == "loob"){
			id.cut <- which.max(sensloob+spezloob-1)
			control$best.cutoff <- thres[id.cut]
		}		
		if(cutoff == "0.632"){
			id.cut <- which.max(sens632+spez632-1)
			control$best.cutoff <- thres[id.cut]
		}
		if(cutoff == "0.632+"){
			id.cut <- which.max(sens632p+spez632p-1)
			control$best.cutoff <- thres[id.cut]
		}
		if(is.character(cutoff)){
			err632p <- err632p[id.cut]
			err632 <- err632[id.cut]
			errloob <- errloob[id.cut]
			errapp <- errapp[id.cut]
		}
		output <- list(call=call, formula=form, method=control,
					   err632p=err632p,err632=err632,errloob=errloob,errapp=errapp,
					   sens632p=sens632p[id.cut], spec632p=spez632p[id.cut],
					   sens632=sens632[id.cut], spec632=spez632[id.cut],
					   sensloob=sensloob[id.cut], specloob=spezloob[id.cut],
					   sensapp=sensapp[id.cut], specapp=spezapp[id.cut],
					   roc=data.frame(sens632p=sens632p,spec632p=spez632p,
									  sens632=sens632,spec632=spez632,
									  sensloob=sensloob,specloob=spezloob,
									  sensapp=sensapp,specapp=spezapp, cut.points=thres),
					   ans=ans,
					   sample.roc=all.roc, sample.data=all.data)
		class(output) <- c("Daim","boot",est.method)
		output
	}
	else{
		if(cutoff == "cv"){
			id.cut <- which.max(sensloob + spezloob - 1)
			control$best.cutoff <- thres[id.cut]
			errloob <- errloob[id.cut]
			errapp <- errapp[id.cut]
		}
		output <- list(call=call, formula=form, method=control,
					   err632p=NA, err632=NA, errloob=errloob, errapp=errapp,
					   sens632p=NA, sens632=NA,
					   sensloob=sensloob[id.cut], specloob=spezloob[id.cut],
					   sensapp=sensapp[id.cut], specapp=spezapp[id.cut],
					   roc=data.frame(sensloob=sensloob, specloob=spezloob,
									  sensapp=sensapp, specapp=spezapp, cut.points=thres),
					   sample.roc=all.roc, sample.data=all.data)
		class(output) <- c("Daim", "cv", est.method)
		output
	}
}









