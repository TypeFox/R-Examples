bivariate.density <- function(data, ID = NULL, pilotH = NULL, globalH = pilotH, adaptive = TRUE, edgeCorrect = TRUE, res = 50, WIN = NULL, counts = NULL, intensity = FALSE, xrange = NULL, yrange = NULL, trim = 5, gamma = NULL, pdef = NULL, atExtraCoords = NULL, use.ppp.methods = TRUE, comment = TRUE){
    
	if(comment) cat(date())
	if(comment) cat("\nperforming basic check of arguments...\n")
	
	data.org.cls <- class(data)
	
	if(data.org.cls=="data.frame"){
		if(ncol(data)!=2&&ncol(data)!=3) stop("data.frame must have exactly 2 or 3 columns")
		if(nrow(data)<10) warning("less than 10 observations!")
	} else if(data.org.cls=="list"){
		if(length(data)==2){
			if(is.null(data$x)||is.null(data$y)) stop("2-component data list must have labels 'x' and 'y'")
			if(length(data$x)!=length(data$y)) stop("data components 'x' and 'y' are of unequal lengths")
		} else if(length(data)==3){
			if(is.null(data$ID)) stop("3-component data list must have labels 'x', 'y' and 'ID'")
			if(length(data$x)!=length(data$y)||length(data$x)!=length(data$ID)) stop("data components 'x', 'y' and 'ID' are of unequal lengths")
		} else {
			stop("data list must have exactly 2 or 3 components")
		}
		
		if(length(data$x)<10) warning("less than 10 observations!")
		data <- as.data.frame(data)
		
	} else if(data.org.cls=="matrix"){
		if(ncol(data)!=2&&ncol(data)!=3) stop("data matrix must have exactly 2 or 3 columns")
		
		if(nrow(data)<10) warning("less than 10 observations!")
		data <- as.data.frame(data)
	} else if(data.org.cls=="ppp"){
		if(is.null(data$x)||is.null(data$y)) stop("data ppp.object contains no observations!")
		if(length(data$x)<10) warning("less than 10 observations!")
		
		WIN <- data$window
		
		if(is.null(data$marks)){
			data <- data.frame(cbind(data$x,data$y))
		} else {
			org.levels <- levels(as.factor(data$marks))
			num.levels <- as.numeric(levels(as.factor(as.numeric(as.factor(data$marks)))))
			data <- data.frame(cbind(data$x,data$y,data$marks))
		}
	} else {
		stop("'data' must be an object of class 'data.frame', 'list', 'matrix', or 'ppp'")
	}
	
	if((!is.null(ID))&&(ncol(data)<3)) stop("data must have 3 columns if 'ID' argument specified")
	
	if(is.null(ID)){
		if(ncol(data)==3){
			data <- data[,1:2]
		}
	} else {
		if(ncol(data)==2){
			warning("no ID information provided in 'data' argument; cannot match 'ID'. Fitting on all present observations")
		} else {
			if(length(unique(data[,3]))>2) warning("more than two distinct ID values detected in 'data'")
			if(data.org.cls=="ppp") ID <- num.levels[which(org.levels==ID)]
			if(nrow(data[data[,3]==ID,])==0) stop("no observations match given 'ID'!")
			data <- data[data[,3]==ID,1:2]
			if(nrow(data)<10) warning("less than 10 observations!")
		}
	}
	
	
	kType <- "gaus"
	
	if(is.null(WIN)&&is.null(xrange)&&data.org.cls!="ppp") stop("if 'data' is not of class 'ppp' then user must supply either 'WIN' or 'xrange' AND 'yrange'")
	
	
	if(is.null(pdef)){
		if(pilotH<=0||is.na(pilotH)||is.null(pilotH)) stop("'pilotH' must be > 0")
	} else {
		if(!adaptive){
			warning("'pdef' object ignored for fixed-bandwidth estimates")
			if(pilotH<=0||is.na(pilotH)||is.null(pilotH)) stop("'pilotH' must be > 0")	
		} else {
			if(class(pdef)!="bivden") stop("'pdef' must be an object of class 'bivden'")
		}
	}
	if(globalH<=0||is.na(globalH)||is.null(globalH)) stop("'globalH' must be > 0")
	if(res<=0||is.na(res)||is.null(res)) stop("invalid resolution")
	if(length(trim)>1){
		trim <- trim[1]
		warning("trimming argument has more than one component... attempting to interpret first only")
	}
	if(!is.na(trim)){
		if(trim<=0) trim <- NA
	}
	
	if(!is.null(counts)){
		counts <- round(counts) ##warning here?
		if(any(counts<=0)) stop("counts must be positive integers")
	}
	
	if(!is.null(WIN)){
		if(class(WIN)!="owin"){
			stop("WIN must be an object of class 'owin'")
		} else {
			xrange <- WIN$xrange
			yrange <- WIN$yrange
		}
	} else {
		if(length(xrange)!=2||length(yrange)!=2){
			stop("no WIN argument supplied - xrange and yrange must both be vectors of length 2")
		} else {
			WIN <- owin(xrange,yrange)
			warning("WIN not specified - using rectangular region defined by xrange and yrange arguments")
		}
	}
	
	
	duplicates <- dupli.data.frame(data,WIN,comment)
	if(!is.null(duplicates$counts)&&!is.null(counts)){
		warning("duplicated coords detected - ignoring user supplied arg 'counts' and creating own")
		counts <- duplicates$counts
		data <- duplicates$data
	} else if(!is.null(duplicates$counts)&&is.null(counts)){
		counts <- duplicates$counts
		data <- duplicates$data
	} else if(is.null(duplicates$counts)&&is.null(counts)){
		counts <- rep(1,nrow(data))
	}
	
	n <- sum(counts)
	if(intensity) nfac <- n
	else nfac <- 1
	
	if(length(counts)!=nrow(data)) stop("'counts' arg must be of equal length to no. of obervations")

    xrg <- seq(xrange[1],xrange[2],length=res)
    yrg <- seq(yrange[1],yrange[2],length=res)
	xdatarange <- sort(rep(xrg,res))
	ydatarange <- rep(yrg,res)
	
	if(sum(!inside.owin(data[,1],data[,2],WIN))>0){
		if(comment) warning(paste("data contain",sum(!inside.owin(data[,1],data[,2],WIN)),"observations outside study region WIN - these were removed"))
		data <- data[-which(!inside.owin(data[,1],data[,2],WIN)),]
	}
	 
	hypoH <- spec_pilot_f_values <- total_pilot_f_values <- extra_pilot_f_values <- NULL
	datarange <- data.frame(cbind(xdatarange,ydatarange))
	datarange.list <- list(x=xrg,y=yrg)
    datarangeNA <- data.frame(cbind(xdatarange,ydatarange))
    datarangeNA[!inside.owin(datarangeNA[,1],datarangeNA[,2],WIN),] <- NA
    
	if(comment) cat("setting up bandwidth(s)...\n")
	
	if(adaptive){
		if(is.null(pdef)){
			if(!use.ppp.methods){
				spec_pilot_f_values <- apply(as.matrix(data),1,KSPEC,data=data,h=pilotH,type=kType,counts=counts)
				total_pilot_f_values <- apply(as.matrix(datarangeNA),1,KSPEC,data=data,h=pilotH,type=kType,counts=counts)
				if(edgeCorrect){
					spec_pilot_f_values <- spec_pilot_f_values/getQhz_Fixed(Xseq=data[,1],Yseq=data[,2],kType=kType,WIN=WIN,h=pilotH,both=F)$qhz
					total_pilot_f_values <- total_pilot_f_values/getQhz_Fixed(Xseq=datarange[,1],Yseq=datarange[,2],kType=kType,WIN=WIN,h=pilotH,both=F)$qhz
				}
			} else {
				pilot.ppp <- ppp(x=data[,1],y=data[,2],window=WIN,check=F)
				pilot.density <- density(pilot.ppp,sigma=pilotH,xy=datarange.list,weights=counts*rep(1/n,nrow(data)),edge=edgeCorrect)
			
				corrGridSpec <- apply(data,1,getNearest,gridx=xdatarange,gridy=ydatarange,WIN=WIN)
				total_pilot_f_values <- as.vector(pilot.density$v)
				total_pilot_f_values[which(is.nan(total_pilot_f_values))] <- NA
				total_pilot_f_values[which(is.infinite(total_pilot_f_values))] <- NA
				total_min <- min(na.omit(total_pilot_f_values)[na.omit(total_pilot_f_values)>0])
				total_pilot_f_values[!is.na(total_pilot_f_values)][total_pilot_f_values[!is.na(total_pilot_f_values)]<=0] <- total_min
				spec_pilot_f_values <- total_pilot_f_values[corrGridSpec]
			}
		
			if(!is.null(atExtraCoords)){
				if(!use.ppp.methods){
					if(edgeCorrect){
						extra_pilot_f_values <- apply(as.matrix(atExtraCoords),1,KSPEC,data=data,h=pilotH,type=kType,counts=counts)/getQhz_Fixed(Xseq=atExtraCoords[,1],Yseq=atExtraCoords[,2],kType=kType,WIN=WIN,h=pilotH,both=F)$qhz
					} else {
						extra_pilot_f_values <- apply(as.matrix(atExtraCoords),1,KSPEC,data=data,h=pilotH,type=kType,counts=counts)
					}
				} else {
					corrGridExtra <- apply(atExtraCoords,1,getNearest,gridx=xdatarange,gridy=ydatarange,WIN=WIN)
					extra_pilot_f_values <- total_pilot_f_values[corrGridExtra]
				}
			}
		} else {
			if(length(xrg)!=length(pdef$X)||length(yrg)!=length(pdef$Y)){
				stop("'pdef' grid resolution must be identical to current value of 'res'")
			} else {
				if(!all(c(xrg,yrg)==c(pdef$X,pdef$Y))) stop("'pdef' must have identical grid coordinates to the current estimate")
			}
			
			if(!identical_windows(WIN,pdef$WIN)) stop("'pdef' window object must be identical to current 'WIN'")
			
			pilotH <- pdef$pilotH
			total_pilot_f_values <- as.vector(t(pdef$Zm))
			corrGridSpec <- apply(data,1,getNearest,gridx=xdatarange,gridy=ydatarange,WIN=WIN)
			spec_pilot_f_values <- total_pilot_f_values[corrGridSpec]
			
			if(!is.null(atExtraCoords)){
				corrGridExtra <- apply(atExtraCoords,1,getNearest,gridx=xdatarange,gridy=ydatarange,WIN=WIN)
				extra_pilot_f_values <- total_pilot_f_values[corrGridExtra]
			}
		} 
	
		if(is.null(gamma)) gamma <- exp(mean(log(1/sqrt(spec_pilot_f_values))))
		
		h <- globalH*sqrt(1/spec_pilot_f_values)*(1/gamma)
		hypoH <- globalH*sqrt(1/total_pilot_f_values)*(1/gamma)
		
		if(!is.na(trim)){
			if(is.null(pdef)){
				beta.hM <- trim*median(h[!is.na(h)])
			} else {
				beta.hM <- trim*median((globalH*sqrt(1/pdef$zSpec)*(1/gamma))[!is.na(globalH*sqrt(1/pdef$zSpec)*(1/gamma))])
			}
			h[!is.na(h)][h[!is.na(h)] > beta.hM] <- beta.hM
			hypoH[!is.na(hypoH)][hypoH[!is.na(hypoH)] > beta.hM] <- beta.hM 
		}
		
		if(!is.null(atExtraCoords)){
			extraH <- globalH*sqrt(1/extra_pilot_f_values)*(1/gamma)
			if(!is.na(trim)){
				extraH[!is.na(extraH)][extraH[!is.na(extraH)] > beta.hM] <- beta.hM 
			}
		}
	} else {
        h <- rep(pilotH,nrow(data))
		if(use.ppp.methods){
			if(comment) cat("calculating density and edge-correcting if elected...\n")
			ppp.den.QEX <- run_ppp(data,datarange.list,pilotH,WIN,counts)
			corrGridSpec <- apply(data,1,getNearest,gridx=xdatarange,gridy=ydatarange,WIN=WIN)
			if(!is.null(atExtraCoords))	corrGridExtra <- apply(atExtraCoords,1,getNearest,gridx=xdatarange,gridy=ydatarange,WIN=WIN)
			
			if(comment) cat("returning...\n")
			
			edg.v <- qhzSpec <- qhzExtra <- 1
			if(edgeCorrect){
				edg.v <- as.vector(ppp.den.QEX$edg$v)
				edg.v[!inside.owin(xdatarange,ydatarange,WIN)] <- NA
				qhzSpec <- edg.v[corrGridSpec]
				if(!is.null(atExtraCoords)){
					qhzExtra <- edg.v[corrGridExtra]
				} else {
					qhzExtra <- NA
				}
			}
			raw.v <- as.vector(ppp.den.QEX$raw$v)
			raw.v[!inside.owin(xdatarange,ydatarange,WIN)] <- NA
						
			zv <- raw.v/edg.v
			zv[zv<=0] <- NA
			Zmat <- matrix(zv,res,res,byrow=T)
			extra <- NA
			if(!is.null(atExtraCoords)){
				extra <- raw.v[corrGridExtra]/qhzExtra
				extra[extra<=0] <- NA
			}
			
			result <- list(	Zm=nfac*Zmat,
						X=xrg,
						Y=yrg,
						kType=kType,
						h=pilotH,
						pilotH=pilotH,
						globalH=NA,
						hypoH=NULL,
						zSpec=nfac*raw.v[corrGridSpec]/qhzSpec,
						zExtra=nfac*extra,
						WIN=WIN,
						qhz=matrix(edg.v,res,res,byrow=T),
						qhzSpec=qhzSpec,
						qhzExtra=qhzExtra,
						pdef=NULL,
						pilotvals=NULL,
						gamma=NA,
						counts=counts,
						data=data)	
			class(result) <- "bivden"
			if(comment) cat(date(),"\n")
			return(result)
			
		}
    }

	if(comment) cat("calculating density and edge-correcting if elected...\n")
    surfA <- apply(as.matrix(datarange),1,compute.coord,data=data,h=h,n=n,WIN=WIN,type=kType,counts=counts)
    surfB <- apply(as.matrix(data),1,compute.coord,data=data,h=h,n=n,WIN=WIN,type=kType,counts=counts)
    surfC <- NA
	QC <- list(qhz=NA,qhz_sq=NA)
    if(!is.null(atExtraCoords)){
		surfC <- apply(as.matrix(atExtraCoords),1,compute.coord,data=data,h=h,n=n,WIN=WIN,type=kType,counts=counts)
		QC <- list(qhz=1,qhz_sq=1)
	}

    QA <- QB <- list(qhz=1,qhz_sq=1)
	
	
	if(edgeCorrect&&adaptive&&use.ppp.methods){
		hypoQuan <- unique(quantile(hypoH,(1:100)/100,na.rm=T))
		corrQuan <- apply(as.matrix(hypoH),1,idQuan,q=hypoQuan)
		qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=data,xy=datarange.list,WIN=WIN,counts=counts)
		qhz <- rep(-1,res*res)
		qhz[is.na(hypoH)] <- NA
		for(i in 1:length(hypoQuan)) qhz[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
		qhzSpec <- qhz[corrGridSpec]
		
		qhzExtra <- NA
		if(!is.null(atExtraCoords)){
			qhzExtra <- qhz[corrGridExtra]
			surfC <- surfC/qhzExtra
			surfC[surfC<=0] <- NA
		}
		
		surfA <- surfA/qhz
		surfA[surfA<=0] <- NA
		surfB <- surfB/qhzSpec
		surfB[surfB<=0] <- NA
		
		if(comment) cat("returning...\n")
		
		result <- list(	Zm=nfac*matrix(surfA,res,res,byrow=T),
					X=xrg,
					Y=yrg,
					kType=kType,
					h=h,
					pilotH=pilotH,
					globalH=globalH,
					hypoH=matrix(hypoH,res,res,byrow=T),
					zSpec=nfac*surfB,
					zExtra=nfac*surfC,
					WIN=WIN,
					qhz=matrix(qhz,res,res,byrow=T),
					qhzSpec=qhzSpec,
					qhzExtra=qhzExtra,
					pdef=pdef,
					pilotvals=spec_pilot_f_values,
					gamma=gamma,
					counts=counts,
					data=data)	
		class(result) <- "bivden"
		if(comment) cat(date(),"\n")
		return(result)
	} else {
		if(edgeCorrect){
			datarangeA <- as.matrix(data.frame(cbind(1:nrow(datarange),hypoH)))
			datarangeB <- as.matrix(data.frame(cbind(1:nrow(data),h)))

			if(adaptive){
				datarangeA <- as.matrix(data.frame(cbind(xdatarange,ydatarange,hypoH)))
				datarangeB <- as.matrix(data.frame(cbind(data,h)))

				QA <- apply(datarangeA,1,getQhz_Adaptive,kType=kType,WIN=WIN,both=F)
				QA <- list(qhz=QA[1,])

				QB <- apply(datarangeB,1,getQhz_Adaptive,kType=kType,WIN=WIN,both=F)
				QB <- list(qhz=QB[1,])
					
				if(!is.null(atExtraCoords)){
					datarangeC <- as.matrix(data.frame(cbind(atExtraCoords,extraH)))
					QC <- apply(datarangeC,1,getQhz_Adaptive,kType=kType,WIN=WIN,both=F)
					QC <- list(qhz=QC[1,])
				}
			} else {
				QA <- getQhz_Fixed(Xseq=xdatarange,Yseq=ydatarange,kType=kType,WIN=WIN,h=pilotH,both=F)
				QB <- getQhz_Fixed(Xseq=data[,1],Yseq=data[,2],kType=kType,WIN=WIN,h=pilotH,both=F)
				
				if(!is.null(atExtraCoords)) QC <- getQhz_Fixed(Xseq=atExtraCoords[,1],Yseq=atExtraCoords[,2],kType=kType,WIN=WIN,h=pilotH,both=F)
			}

			surfA <- surfA/QA$qhz
			surfB <- surfB/QB$qhz

			if(!is.null(atExtraCoords)) surfC <- surfC/QC$qhz
		}
	}
	
	if(comment) cat("returning...\n")
	if(comment) cat(date(),"\n")
  
    Zm <- matrix(surfA,res,res,byrow=T)
	if(!adaptive){
		hypoH <- NA
		h <- pilotH
	} else {
		hypoH <- matrix(hypoH,res,res,byrow=T)
	}
	
    result <- list(Zm=nfac*Zm,
                X=xrg,
                Y=yrg,
                kType=kType,
                h=h,
                pilotH=pilotH,
                globalH=globalH,
                hypoH=hypoH,
                zSpec=nfac*surfB,
                zExtra=nfac*surfC,
                WIN=WIN,
                qhz=matrix(QA$qhz,res,res,byrow=T),
                qhzSpec=QB$qhz,
                qhzExtra=QC$qhz,
                pdef=pdef,
                pilotvals=spec_pilot_f_values,
                gamma=gamma,
				counts=counts,
                data=data)
	class(result) <- "bivden"
	return(result)
}