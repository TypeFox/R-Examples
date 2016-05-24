#NOTES
#Need to write up manuals!!!!
#Reliant on matrixStats and mvnfast (for speed improvement)
#Needs user-error catches (modelnames done, Qs done, dfupdate done)

mmtfaEM <- function(x, Gs, Qs, clas, init, scale, models, 
        dfstart, dfupdate, gauss, eps, known){
 
  modgen <- modelgen()
  modold <- modgen$modold
  p <- ncol(x)
  n <- nrow(x)
  
  zlist3 <- ll <- dff <- it <- store <- meanlist <-  siglist <- list()
	if(clas>0){
			testindex <- sample(1:n, ceiling(n*(clas/100)))
			#print(testindex)
			kno <- vector(mode="numeric", length=n)
			kno[testindex] <- 1
			unkno <- (kno-1)*(-1)
			Gs <- length(unique(known))
	}
	
	#BIC and ICL arrays
  gvec <- 1:max(Gs)
  qvec <- 1:max(Qs)
  gstuff <- paste("G=",gvec,sep="")
  qstuff <- paste("Q=",qvec,sep="")
  bic <- rands <- icl <- logls <- array(-Inf, dim=c(length(models), max(Qs), max(Gs)))
  meansave <- meansave2 <- sigsave <- zmatsave <- sigsave2 <- zmatsave2 <- NA
	

  oldmodvec <- modold[match(models, modgen$allmodels)]
	#####ZMAT INITIALIZATION
  if(class(init)!="list" && !(init %in% c("kmeans", "hard", "disc", "cont", "soft", "uniform"))){
    stop("'init' must be one of 'kmeans', 'hard', 'soft', 'uniform' or a list. See ?mmtfa.")
  }
	zmatin <- list()
	for(G in Gs){
		if(G==1){
			zmatin[[G]] <- matrix(1,n,1)
		}
		else{
		  if(is.character(init)){
  			if(init == "hard"){
  				zmatin[[G]] <- discrandz(n,G)
  			}
  			if(init == "soft"){
  				zmatin[[G]] <- contrandz(n,G)
  			}
  # 			if(init == "agglom"){
  # 				zmatin[[G]] <- agglomz(x,n,G,clas,kno,known,testindex)
  # 			}
		    if(init == "uniform"){
		      if(clas>0){
    				zmatin[[G]] <- uniformz(n,G,clas,kno,known)
    			}
  		    else{
  		      stop("Uniform initialization not available for clustering.")
  		      return(NULL)
  		    }
		    }
  # 			if(init == "given"){
  # 				zmatin[[G]] <- givenz(n,G,known)
  # 			}
  			if(init == "kmeans"){
  				zmatin[[G]] <- kmeansz(x,n,G)
  			}
		  }
		  else{
		    zmatin[[G]] <- givenz(n,G,init[[G]])
		  }
		}
	}
		
############## LOOP THROUGH MODELS, Q, and G ##############
for(modnum in 1:length(models)){
	modnew <- models[modnum]
  mod <- modold[which(modgen$allmodels==modnew)]
  
	#print(mod)
	#zlist2 <- it2 <- df2 <- ll2 <- meanlist2 <-  siglist2 <- list()
	for(G in Gs){
	  #OBJECT INITIALIZATION
	  delta <- matrix(0,n,G)
	  mug <- matrix(0,G,p)
	  om <- rep(0,G)
	  yg <- sigma <- sigmainv <- sg <- array(0, dim=c(p, p, G))
	  tri <- array(0,dim=c(p,p,G))
	  w <- matrix(0,n,G)
	  
	  
		#print(G)
		#df1 <- it1 <- zlist1 <- ll1 <- meanlist1 <-  siglist1 <- list()
		#where the z-stuff was before
		for(q in Qs){
      #OBJECT INITIALIZATION
		  lg <- array(0, dim=c(p, q, G))
		  betag <- array(0, dim=c(q, p, G))
		  thetag <- array(0, dim=c(q, q, G))
      
			#print(q)
###########################################################			
#################START OF ALGORITHM########################
singular <- 0
breakit <- 0
if(G==1){
  #will need to fix
  CCCCgroup <- c("UUCU", "UUCC", "UCCU", "UCCC", "CUCU", "CUCC", "CCCU", "CCCC")
  #CCCCgroup <- modgen$allmodels[match(CCCColdgroup, modold)]
  #
  if(any(mod==CCCCgroup)){
    cccdum <- oldmodvec[oldmodvec %in% CCCCgroup]
    if(length(cccdum)>0){
      if(mod!=cccdum[1]){
        breakit <- 1
      }
    }
  }
  #will need to fix
  CCUCgroup <- c("UUUU", "UUUC", "UCUU","UCUC", "CUUU", "CUUC","CCUU", "CCUC",
                 "Mt1U", "Mt1C", "Mt2U", "Mt2C", "Mt3U", "Mt3C", "Mt4U", "Mt4C")
  #CCUCgroup <- modgen$allmodels[match(CCUColdgroup, modold)]
  #
  if(any(mod==CCUCgroup)){
    ccudum <- oldmodvec[oldmodvec %in% CCUCgroup]
    if(length(ccudum)>0){
      if(mod!=ccudum[1]){
        breakit <- 1
      }
    }
  }
}
			
if(breakit==0){
#OTHER PARAMETER INITIALIZATION
	zmat <- zmatin[[G]]
	vg <- vginit(dfstart,G)
	ng <- ngupdate(zmat)
	pig <- pigupdate(ng,n)
	mug <- muginit(G,p,x,zmat,ng)
	sg <- sginit(p,G,x,mug,zmat,n,ng)
	sgc <- sginitc(G,sg,pig,p,n,x)
	if(substring(mod,1,3)=="CCC"){
		for(g in 1:G){
			sg[,,g] <- sgc 
		}
	}

#LAMBDA INITIALIZATION
	if(substring(mod,1,1) == "U"|substring(mod,3,3) == "2"|substring(mod,3,3) == "4"){
		lg <- lginitu(p,q,G,sg)
	}
	if(substring(mod,1,1) == "C"|substring(mod,3,3) == "1"|substring(mod,3,3) == "3"){
		lg <- lginitc(p,q,G,sgc)
		dumg <- lginitu(p,q,G,sg)	
	}

#PSI INITIALIZATION
	if(substring(mod,3,3) == "U"){
		yg <- yginitu(p,G,sg,lg,mod,pig,dumg)
	}
	if(substring(mod,3,3) == "C"){
		yg <- yginitc(p,G,sg,lg,sgc,mod,pig)	
	}
	if(substring(mod,1,1) == "M"){
		yg <- array(0,dim=c(p,p,G))
		om <- rep(0,G)
		tri <- array(0,dim=c(p,p,G))
		ygst <- yginitu(p,G,sg,lg,mod,pig,dumg)
		if(substring(mod,3,3)=="1"|substring(mod,3,3)=="2"){
			for(g in 1:G){
				om[g] <- det(ygst[,,g])^(1/p)
			}
		}
		if(substring(mod,3,3)=="3"|substring(mod,3,3)=="4"){
			dom <- 0
			for(g in 1:G){
				dom <- dom + pig[g]*det(ygst[,,g])^(1/p)
			}
			om[] <- dom
		}
		for(g in 1:G){
			tri[,,g] <- ygst[,,g]/(det(ygst[,,g])^(1/p))
		}
		if(substring(mod,3,3)=="1"|substring(mod,3,3)=="2"){
			av <- diag(p)-diag(p)
			for(g in 1:G){
				av <- av + pig[g]*tri[,,g]
			}
			tri[,,] <- av
		}
		for(g in 1:G){
			yg[,,g] <- om[g] * tri[,,g]
		}
	}
	yginv <- yginvup(p,G,yg)

#OTHER PARAMETER INITIALIZATION
	sigma <- sigmaup(p,G,lg,yg,sigma)
	###CHECK FOR SINGULARITY
				#for(g in 1:G){
				#	test <- rcond(diag(q) + t(lg[,,g]) %*% yginv[,,g] %*% lg[,,g])
				#	if(test <= .Machine$double.eps){
				#			message("Warning: Singularity at G=");
				#			print(G);
				#			singular <- 1
				#	}
				#}
				#if(singular==1){break}
			testing <- try(sigmainv <- sigmainvup(p,G,yginv,lg,q,sigmainv),silent=TRUE)
				#if(all(is.finite(testing))){
				#	sigmainv <- sigmainvup(p,G,yginv,lg,q)
				#}
				#else{break}
			if(!all(is.finite(testing))){
        break
			}

	#sigmainv <- sigmainvup(p,G,yginv,lg,q)
	betag <- betagup(q,p,G,lg,sigmainv,betag)
	thetag <- thetagup(q,G,betag,lg,sg,thetag)
	w <- winit(x,n,G,mug,sigmainv,vg,p,sg,zmat)

### START OF AECM
	cycle <- 0
	dfnewg <- vg
}
	conv <- 0
	num <- matrix(0,n,G)
	ft <- matrix(0,n,G)
	logl <- NaN

	while(conv != 1){
		if(breakit==1){break}
		#CM-STEP 1
		#UPDATE PI, MU, DF
			ng <- ngupdate(zmat)
			pig <- pigupdate(ng,n)
			mug <- mugupdate(G,zmat,w,x,p,mug,n)
			if(dfupdate=="approx"){
				testing <- try(dfnewg <- dfupdatefun2(mod,dfnewg,ng,zmat,w,G,p,n,x,mug,sigmainv),silent=TRUE)
				if(!all(is.finite(testing))){	
					break
				}
#				else{break}
			}
			if(dfupdate=="numeric"){
			  testing <- try(dfnewg <- dfupdatefun(mod,dfnewg,ng,zmat,w,G,p,n,x,mug,sigmainv),silent=TRUE)
			  if(!all(is.finite(testing))){	
			    break
			  }
			  #				else{break}
			}
# 		#E-STEP 1
# 		#IF NOT 1st ITERATION
# 		#UPDATE Z, W
# 			#if(cycle>0){
#         delta <- deltaup(x,mug,sigma,sigmainv,G,n)
# 				zmat <- zupdate(x,G,pig,dfnewg,p,yg,q,betag,lg,mug,sigmainv,n,clas,kno,known,unkno,delta)
# #print(zmat)
# 				if(any(is.nan(zmat))){
# 					break
# 				}
# 				w <- wupdate(x,n,G,mug,sigmainv,dfnewg,p,delta)
# 			#}

		#INTERMEDIATE STEP
		#UPDATE S, BETA, THETA
			ng <- ngupdate(zmat)
			sg <- sgupdate(p,G,n,x,mug,zmat,w,ng,mod,pig,sg)
			betag <- betagup(q,p,G,lg,sigmainv,betag)
			thetag <- thetagup(q,G,betag,lg,sg,thetag)

		#CM-STEP 2
		#UPDATE LAMBDA, PSI, BETA, THETA
			testing <- try(lg <- lgupdate(mod,p,q,G,ng,yginv,sg,betag,thetag,om,tri,lg),silent=TRUE)
				if(!all(is.finite(testing))){	
					break
				}
		if(substring(mod,1,1)=="M"){
			om <- omupdate(mod,q,G,yg,p,sg,lg,betag,thetag,pig,om,tri)
			tri <- triupdate(mod,q,G,yg,p,sg,lg,betag,thetag,pig,om,tri,ng)
			for(g in 1:G){
				yg[,,g] <- om[g] * tri[,,g]
			}
		}
		else{
			yg <- ygupdate(mod,q,G,yg,p,sg,lg,betag,thetag,pig)
		}
			yginv <- yginvup(p,G,yg)
			sigma <- sigmaup(p,G,lg,yg,sigma)
			###CHECK FOR SINGULARITY
				#for(g in 1:G){
				#	test <- rcond(diag(q) + t(lg[,,g]) %*% yginv[,,g] %*% lg[,,g])
				#	if(test <= .Machine$double.eps){
				#			message("Warning: Singularity at G=");
				#			print(G);
				#			singular <- 1
				#	}
				#}
				#if(singular==1){break}
      
			testing <- try(sigmainv <- sigmainvup(p,G,yginv,lg,q,sigmainv),silent=TRUE)
				#if(all(is.finite(testing))){
				#	sigmainv <- sigmainvup(p,G,yginv,lg,q)
				#}
				#else{break}
			if(!all(is.finite(testing))){
        break
			}
#			else{break}
			betag <- betagup(q,p,G,lg,sigmainv,betag)
			thetag <- thetagup(q,G,betag,lg,sg,thetag)

		
		#E-STEP 2
		#UPDATE Z, W
      delta <- deltaup(x,mug,sigma,sigmainv,G,n,delta)
			suppressWarnings(zup <- zupdate(x,G,pig,dfnewg,p,yg,q,betag,lg,mug,sigmainv,n,clas,kno,known,unkno,delta))
#print(zmat)
      zmat <- zup$zmat
			if(any(is.nan(zmat))){
				break
			}
			w <- wupdate(x,n,G,mug,sigmainv,dfnewg,p,delta,w)
			

		#CONVERGENCE CHECK
		#	ng <- ngupdate(zmat)
			cycle <- cycle + 1
			#print(cycle)

# 			for(g in 1:G){
# 				#ft[,g] <- pig[g]*gamma((dfnewg[g]+p)/2)*(det(yg[,,g])/
# 				#	det(diag(q)-betag[,,g]%*%lg[,,g]))^(-1/2)/((pi*dfnewg[g])^(p/2)*gamma(dfnewg[g]/2)*(1
# 				#	+ mahalanobis(x, mug[g,], sigmainv[,,g], inverted=TRUE)/dfnewg[g])^((dfnewg[g]+p)/2))
# 				#log((det(yg[,,g])/det(diag(q)-betag[,,g]%*%lg[,,g])))
# 				ft[,g]<-log(pig[g])+lgamma((dfnewg[g]+p)/2)-(1/2)*(sum(log(diag(yg[,,g])))-log(det(diag(q)-betag[,,g]%*%lg[,,g])))-
# 					   ((p/2)*(log(pi)+log(dfnewg[g]))+lgamma(dfnewg[g]/2)+
#              ((dfnewg[g]+p)/2)*(log(1+ (delta[,g]/dfnewg[g]))))
# 			}
# 		#	kcon <- -apply(ft,1,max)
#     #  kcon <- -rowMaxs(ft)
# 			kcon <- zup$kcon
# 			ft <- ft + kcon
			
			logl[cycle]<- sum(log(rowSums(zup$num))) - sum(zup$kcon)
		#	logl[cycle]<- sum(log(rowSums(zup$num))) - sum(zup$kcon)
			if(is.na(logl[cycle])){break}
			if(cycle>3){
				if(is.finite(logl[cycle-2])){
					ak <- (logl[cycle]-logl[cycle-1])/(logl[cycle-1]-logl[cycle-2])
					linf <- logl[cycle-1] + (logl[cycle]-logl[cycle-1])/(1-ak)
					if(abs(linf-logl[cycle-1]) < eps){
						conv<-1
					}
					if((logl[cycle]-logl[cycle-1])<0){
						#message("log likelihood decrease at G=");
						#print(G);
						#message("and Q=");
						#print(q);
						break
					}
				}
				else{break}
   		}
		#print(paste("cycle ",cycle,": ",logl[cycle],sep=""))
	}
	#Store z-matrix for this q
	#BIC and ICL calculations
	if(conv==1){
		bic[modnum,q,G] <- bicdum <- BICcalc(conv,G,p,mod,q,logl,n,gauss)
		icl[modnum,q,G] <- icldum <- ICLcalc(conv,n,zmat,bic,modnum,q,G)
    #meanlist1[[q]] <- mug
		#zlist1[[q]] <- zmat
		#siglist1[[q]] <-  sigma
    if(bicdum==max(bic)){
      meansave <- mug
      sigsave <- sigma
      zmatsave <- zmat
      dfsave <- dfnewg
      itsave <- cycle
      llsave <- logl[cycle]
    }
    if(icldum==max(icl)){
      meansave2 <- mug
      sigsave2 <- sigma
      zmatsave2 <- zmat
      dfsave2 <- dfnewg
      itsave2 <- cycle
      llsave2 <- logl[cycle]
    }
		#df1[[q]] <- dfnewg
		#it1[[q]] <- cycle
		#ll1[[q]] <- logl[[cycle]]
	}
###########################################################
###########################################################
			}
			#Store a list of all q z-matrices for this G
			#it2[[G]] <- it1
			#df2[[G]] <- df1
			#zlist2[[G]] <- zlist1
			#ll2[[G]] <- ll1
      #meanlist2[[G]] <- meanlist1
      #siglist2[[G]] <- siglist1
		}
		#Store a list of all G and q z-matrices for this model
		#it[[modnum]] <- it2
		#zlist3[[modnum]] <- zlist2
		#dff[[modnum]] <- df2
		#ll[[modnum]] <- ll2 
    #meanlist[[modnum]] <- meanlist2
    #siglist[[modnum]] <- siglist2
	}
	#rownames(bic) <- models
	#colnames(bic) <- qstuff[1:maxQ]
	dimnames(bic) <- list(models,qstuff,gstuff)
	dimnames(icl) <- list(models,qstuff,gstuff)
	#rands <- array(-Inf, dim=c(length(models), maxQ, maxG))
	#dimnames(rands) <- list(models,qstuff,gstuff)
					
	#Return maximum
	maxes <- which(bic==max(bic), arr.ind=TRUE)
	maxicl <- which(icl==max(icl), arr.ind=TRUE)
	if(nrow(maxes)>1){
		message("WARNING: Maximum BIC tie between two or more models")
		bestmodnum <- maxes[1:nrow(maxes),1]
		bestmod <- models[bestmodnum]
		bestq <- maxes[1:nrow(maxes),2]
		bestg <- maxes[1:nrow(maxes),3]
		itf <- "MULTIPLE"
		dff1 <- "MULTIPLE" 
		bestz <- "MULTIPLE"
		bestzmap <- "MULTIPLE"
#		adjrand <- "MULTIPLE"
		tab <- "MULTIPLE"
		blogl <- "Multiple"
	}
	if(nrow(maxes)==1){
		bestmodnum <- maxes[1]
		bestmod <- models[bestmodnum]
		bestq <- maxes[2]
		bestg <- maxes[3]
		bestz <- zmatsave
		dff1 <- dfsave
		itf <- itsave
		blogl <- llsave
		bestzmap <- apply(bestz, 1, which.max)
		if(clas>0){
			newmap <- bestzmap
			newmap[testindex] <- NA
			newknown <- known
			newknown[testindex] <- NA
			tab <- table(known,newmap)
		}
		else{
		  if(!is.null(known)){ tab <- table(known,bestzmap)}
		  else{tab <- NULL}
		}
		#adjrand <- classAgreement(tab)$crand
	}
	if(nrow(maxicl)>1){
		message("WARNING: Maximum ICL tie between two or more models")
		bestmodnumicl <- maxicl[1:nrow(maxicl),1]
		bestmodicl <- models[bestmodnumicl]
		bestqicl <- maxicl[1:nrow(maxicl),2]
		bestgicl <- maxicl[1:nrow(maxicl),3]
		dff1icl <- "MULTIPLE"
		bestzicl <- "MULTIPLE"
		bestzmapicl <- "MULTIPLE"
	#	adjrandicl <- "MULTIPLE"
		itficl <- "MULTIPLE"
		tabicl <- "MULTIPLE" 
		bloglicl <- "MULTIPLE"
	}
	if(nrow(maxicl)==1){
		bestmodnumicl <- maxicl[1]
		bestmodicl <- models[bestmodnumicl]
		bestqicl <- maxicl[2]
		bestgicl <- maxicl[3]
#		bestzicl <- zlist3[[bestmodnumicl]][[bestgicl]][[bestqicl]]
    bestzicl <- zmatsave2
		dff1icl <- dfsave2
		itficl <- itsave2
		bloglicl <- llsave2
		bestzmapicl <- apply(bestzicl, 1, which.max)
		if(clas>0){
			newmapicl <- bestzmapicl
			newmapicl[testindex] <- NA
			newknown <- known
			newknown[testindex] <- NA
			tabicl <- table(known,newmapicl)
		}
		else{
		  if(!is.null(known)){ tabicl <- table(known,bestzmapicl)}
		  else{ tabicl <- NULL}
		}
	#	adjrandicl <- classAgreement(tabicl)$crand  
	}
# 	if(bestmod!=bestmodicl){
# 		message("WARNING: disagreement between BIC and ICL")
# 	}
#   else{
#   	if(bestg!=bestgicl){
#   		message("WARNING: disagreement between BIC and ICL")
#   	}
#     else{
#       if(bestq!=bestqicl){
#         message("WARNING: disagreement between BIC and ICL")
#       }
#     }
#   }
  iclresults <- list()
  par <- list()
  paricl <- list()
  par[["mean"]] <- meansave
# par[["sigma"]] <- siglist[[bestmodnumicl]][[bestgicl]][[bestqicl]]
  par[["sigma"]] <- sigsave
  par[["df"]] <- dff1
  paricl[["sigma"]] <- sigsave2
  paricl[["mean"]] <- meansave2
  paricl[["df"]] <- dff1icl
  store[["parameters"]] <- par
#	store[["z"]] <- zlist3
# store[["z"]] <- zmatsave
	store[["allbic"]] <- bic[, Qs, Gs]
	iclresults[["allicl"]] <- icl[, Qs, Gs]
	store[["bic"]] <- max(bic)
	iclresults[["icl"]] <- max(icl)
	store[["modelname"]] <- bestmod
	store[["bestmodel"]] <- paste("The best model (BIC of ",round(max(bic),2),") is ",bestmod," with G=",bestg,sep="")
	store[["Q"]] <- bestq
	store[["G"]] <- bestg
	store[["classification"]] <- bestzmap
	iclresults[["bestmodel"]] <- paste("The best model (ICL of ",round(max(icl),2),") is ",bestmodicl," with G=",bestgicl,sep="")
	iclresults[["modelname"]] <- bestmodicl
	iclresults[["Q"]] <- bestqicl
	iclresults[["G"]] <- bestgicl
	iclresults[["fuzzy"]] <- bestzicl
	iclresults[["logl"]] <- bloglicl
	iclresults[["classification"]] <- bestzmapicl
	#iclresults[["rand"]] <- adjrandicl
  iclresults[["parameters"]] <- paricl 
#	store[["rand"]] <- adjrand
	store[["tab"]] <- tab
	iclresults[["tab"]] <- tabicl
	store[["iter"]] <- itf
	iclresults[["iter"]] <- itficl
#	store[["alldf"]] <- dff
	store[["x"]] <- x
	store[["fuzzy"]] <- bestz
	store[["logl"]] <- blogl
  store[["iclresults"]] <- iclresults
  #store[["index"]] <- testindex
	store
}
