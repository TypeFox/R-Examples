skatOMeta <- function(..., SNPInfo=NULL, skat.wts = function(maf){dbeta(maf,1,25)}, burden.wts = function(maf){as.numeric(maf <= 0.01) }, rho=c(0,1), method = c("integration", "saddlepoint", "liu"), snpNames = "Name", aggregateBy = "gene", mafRange = c(0,0.5), verbose=FALSE){
	cl <- match.call(expand.dots = FALSE)
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	}
	
	if(any(rho >1 | rho < 0 ) ) stop("rho must be between 0 and 1")
	method <- match.arg(method)
   #if( !(method %in% c("davies","farebrother","imhof","liu")) ) stop("Method specified is not valid! See documentation")
	
	genelist <- na.omit(unique(SNPInfo[,aggregateBy]))
	cohorts <- list(...)	
	ncohort <- length(cohorts)
	
	if(any(unlist(lapply(cohorts,class)) != "skatCohort")){
	 	stop("argument to ... is not a skatCohort object!")
	}
	
	res.strings <- data.frame("gene"=genelist,stringsAsFactors=F)
	res.numeric <- matrix(NA, nrow= nrow(res.strings),ncol =  length(c("p","pmin","rho","cmaf","nmiss", "nsnps", "errflag")))
	colnames(res.numeric) <- c("p","pmin","rho","cmaf","nmiss", "nsnps","errflag")		
	
    if(verbose){
    	cat("\n Meta Analyzing... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = length(genelist), style = 3)
    	pb.i <- 0
    }
    ri <- 0
    snp.names.list <- split(SNPInfo[,snpNames],SNPInfo[,aggregateBy])
	for(gene in genelist){		
		
		ri <- ri+1
		nsnps.sub <- length(snp.names.list[[gene]])
		
		mscores <- maf <- numeric(nsnps.sub)
		big.cov <- matrix(0, nsnps.sub,nsnps.sub)
		n.total <- numeric(nsnps.sub)
		n.miss <- numeric(nsnps.sub)
		
		vary.ave <- 0
		for(cohort.k in 1:ncohort){
			cohort.gene <- cohorts[[cohort.k]][[gene]]
			
			if(!is.null(cohort.gene)){
				sub <- match(snp.names.list[[gene]],colnames(cohort.gene$cov))
				if(any(is.na(sub)) | any(sub != 1:length(sub), na.rm=TRUE) | length(cohort.gene$maf) > nsnps.sub){
							#if(any(is.na(sub))) warning("Some SNPs were not in SNPInfo file for gene ", gene," and cohort ",names(cohorts)[cohort.k])
							cohort.gene$cov <- cohort.gene$cov[sub,sub,drop=FALSE]
							cohort.gene$cov[is.na(sub),] <- cohort.gene$cov[,is.na(sub)] <- 0
							
							cohort.gene$maf <- cohort.gene$maf[sub]
							cohort.gene$maf[is.na(sub)] <- -1
							
							cohort.gene$scores <- cohort.gene$scores[sub]
							cohort.gene$scores[is.na(sub)] <- 0
					}				
					
					n.total[cohort.gene$maf >= 0] <- n.total[cohort.gene$maf >= 0]+cohort.gene$n
					n.miss[cohort.gene$maf < 0] <- n.miss[cohort.gene$maf < 0] + cohort.gene$n
					cohort.gene$maf[cohort.gene$maf < 0] <- 0
					
					mscores <- mscores + cohort.gene$scores/cohort.gene$sey^2
					maf <- maf + 2*cohort.gene$maf*(cohort.gene$n)
					big.cov <- big.cov + cohort.gene$cov/cohort.gene$sey^2
					vary.ave <- vary.ave + max(cohort.gene$n,na.rm=T)*cohort.gene$sey^2
			}else{
				n.miss <- n.miss + cohorts[[cohort.k]][[1]]$n
			} 
		}
		if(any(maf >0)){ 
			maf <- maf/(2*n.total)
			maf[is.nan(maf)] <- 0
			maf <- sapply(maf, function(x){min(x,1-x)})
		
			if( !all(mafRange == c(0,0.5))){
				keep <- (maf >= min(mafRange)) & (maf <= max(mafRange))
				
				big.cov <- big.cov[keep,keep]
				mscores <- mscores[keep]
				maf <- maf[keep]
			}
		}
		
		if(is.function(skat.wts)){
			w1 <- skat.wts(maf)
		} else if(is.character(skat.wts)){
			w1 <- as.numeric(SNPInfo[SNPInfo[,aggregateBy]==gene,skat.wts])
		} else {
			w1 <- rep(1,length(maf))
		}
		
		if(is.function(burden.wts)){
			w2 <- burden.wts(maf)
		} else if(is.character(burden.wts)){
			w2 <- as.numeric(SNPInfo[SNPInfo[,aggregateBy]==gene,burden.wts])
		} else {
			w2 <- rep(1,length(maf))
		}
		
		w1 <- ifelse(maf >0, w1,0)
		w2 <- ifelse(maf >0, w2,0)
		
		##
		Q.skat <- sum((w1*mscores)^2, na.rm=TRUE)
		V.skat <- (w1)*t(t(big.cov)*as.vector(w1))

		Q.burden <- sum(w2*mscores, na.rm=TRUE)^2
		V.burden <- t(w2)%*%big.cov%*%w2
		
		#If burden test is 0, or only 1 SNP in the gene, do SKAT:
		if(sum(maf > 0) ==1 | V.burden ==0){
			lambda <- eigen(zapsmall(V.skat), symmetric = TRUE)$values
        	if(any(lambda > 0) & length(lambda) >1) {
        		tmpP <- pchisqsum2(Q.skat,lambda=lambda,method=method, acc=1e-7)
        		if(tmpP$errflag !=0 ){ 
        			res.numeric[ri,"errflag"] = 1
        		} else {
        			res.numeric[ri,"errflag"] = 0
        		}
        		p <- tmpP$p
        	} else {
            	p <- ifelse(length(lambda) == 1 & all(lambda) > 0, pchisq(Q.skat/lambda,df=1,lower.tail=FALSE),1)
            	res.numeric[ri,"errflag"] = 0
        	}
 			res.numeric[ri,"pmin"] = res.numeric[ri,"p"] = p
			res.numeric[ri,"rho"] = 0
			
		
		#Else do SKAT-O
		} else {
			skato.res <- skatO_getp(mscores, big.cov, diag(w1), w2, rho, method= method, gene=gene)

	    	res.numeric[ri,"p"] <- skato.res$actualp
			res.numeric[ri,"pmin"] = skato.res$minp
			res.numeric[ri,"rho"] = skato.res$rho    	
    		res.numeric[ri, "errflag"] = skato.res$errflag
    	}
    			
		res.numeric[ri,"cmaf"] = sum(maf,na.rm=TRUE)
		res.numeric[ri,"nsnps"] = sum(maf!= 0, na.rm =T)
		res.numeric[ri,"nmiss"] = sum(n.miss, na.rm =T)
		if(verbose){
			pb.i <- pb.i+1
			setTxtProgressBar(pb, pb.i)
		}
	}
	if(verbose) close(pb)
	return(cbind(res.strings,res.numeric))
}


skatO_getp <- function(U,V, R, w, rho,method = "davies", gene=NULL){
	##Input:
	#U: score vector (length p)
	#R: p x p weight matrix for skat
	#w: burden weights
	#rho: vector of rhos in [0,1] 
	#method: method for calculating Normal quadratic form distribution
	#gene: The name of the region - used for error reporting
	
	##Output: a list with elemeents
	#minp: the minimum p-value
	#actualp: the actual p-value
	#rho: the value of rho which gave the minp
	#ps: the whole vector of p-values
	#errflag: 0 if no problem, 1 if quantile issue, 2 if integration issue
	
	satterthwaite <- function(a, df) {
        if (any(df > 1)) {
        	a <- rep(a, df)
        }
        tr <- mean(a)
        tr2 <- mean(a^2)/(tr^2)
        list(scale = tr * tr2, df = length(a)/tr2)
    }
	
	errflag = 0
	Q.skat <- crossprod(R%*%U) # SKAT 
	Q.burden <- (t(w)%*%U)^2 # burden
	Qs <- (1-rho)*Q.skat + rho*Q.burden
			
	lambdas <- ps <- NULL
	ps <- numeric(length(rho))
	for(i in 1:length(rho)){
		PC <- eigen((1-rho[i])*crossprod(R)+ rho[i]*outer(w,w),symmetric=TRUE)
		v.sqrt <- with(PC,{ values[values < 0] <- 0; (vectors)%*%diag(sqrt(values))%*%t(vectors) })
		lam <- eigen( zapsmall(v.sqrt%*%V%*%v.sqrt),only.values=TRUE )$values
		lam <- lam[lam != 0]
			
		lambdas <- c(lambdas, list( lam ))
		tmpP <- pchisqsum2(Qs[i],lambda=lambdas[[i]],method=method, acc=1e-7)
		if(tmpP$errflag != 0){
			errflag <- 1
			ps[i] <- pchisqsum2(Qs[i],lambda=lambdas[[i]],method="liu")$p
		} else {
			ps[i] <- tmpP$p	
		}
	}
			
	minp <- min(ps)
	Ts <- numeric(length(rho))
			
	for(i in 1:length(rho)){					
		sat <- satterthwaite(lambdas[[i]],rep(1,length(lambdas[[i]])))
		upper <- qchisq(minp/20,df=sat$df,lower.tail=FALSE)*sat$scale		
		tmpT <- try(uniroot(function(x){pchisqsum2(x,lambda=lambdas[[i]],method=method,acc=1e-5)$p- minp }, interval=c(1e-10,upper))$root, silent = TRUE)
		if(class(tmpT) == "try-error"){
			#warning(paste0("Problem finding quantiles in gene ", gene, ", p-value may not be accurate"))
			Ts[i] <- Qs[i]
			errflag <- 2
		} else {
			Ts[i] <- tmpT
		}
	}
	
	v11 <- R%*%V%*%R
	v12 <- R%*%V%*%w	
	v22 <- t(w)%*%V%*%w
	V.cond <- v11 - outer( c(v12), c(v12) )/c(v22)
	lambda.cond <- eigen(V.cond,only.values=TRUE,symmetric=TRUE)$values
			
	EDec <- eigen(V.cond,symmetric=TRUE)
	D <- diag(EDec$values)
	diag(D)[zapsmall(diag(D)) > 0] <- 1/sqrt(diag(D)[zapsmall(diag(D)) > 0])
	diag(D)[diag(D) <= 0 ] <- 0 
	meanvec <- D%*%(EDec$vectors)%*%(v12)/c(v22)
			
	Fcond <- function(x,method){
		pp <- qmax <- numeric(length(x))
		for(i in 1:length(x)){
				qmax[i] <- min( ( (Ts[rho !=1 ] - rho[rho != 1]*x[i])/(1-rho)[rho !=1]) )
			if(any(x[i] > Ts[rho == 1]) ){ 
				pp[i] <- 1
			} else {
				p.tmp <- pchisqsum2(qmax[i], lambda=lambda.cond, delta = c(meanvec)^2*x[i], method = method, acc=min(minp,1e-5) )
				if(p.tmp$errflag != 0) stop("Error in integration! using Liu p-value")
				pp[i] = p.tmp$p
			}
		}			
		return(pp)
	}
			
	if(any(lambda.cond > 0)){
		integrand <- function(x){dchisq(x,1)*Fcond(x*v22,method=method)}
		integral <- try(integrate(Vectorize(integrand),lower=0,upper=Inf, rel.tol=min(minp,1e-4)), silent = TRUE)
		if (class(integral) == "try-error" ) {
       		integrand <- function(x){dchisq(x,1)*Fcond(x*v22,method="liu")}
       		integral <- integrate(Vectorize(integrand),lower=0,upper=Inf)
       		errflag <- 3
       	} else {
       		if(integral$message != "OK") errflag <- 2
       	}
    	actualp <- integral[1]$value
    } else {
    	#cat(".")
    	actualp = minp
    }	
	return(list("actualp"= actualp, "minp" = minp, "rho" = rho[which.min(ps)], "ps" = ps, "errflag" = errflag))
}


