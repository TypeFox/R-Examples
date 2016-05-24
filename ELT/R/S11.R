

	## ------------------------------------------------------------------------ ##
	##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
	##                 Adjustment to a reference"                               ##
	## ------------------------------------------------------------------------ ##
	##  Script        S11.R                                                     ##
	## ------------------------------------------------------------------------ ##
	##  Description   Computation of confidence intervals and                   ##
	##                relative dispersion of the periodic life expectancies     ##
	## ------------------------------------------------------------------------ ##
	##  Authors       Tomas Julien, Frederic Planchet and Wassim Youssef        ##
	##                julien.tomas@univ-lyon1.fr                                ##
	##                frederic.planchet@univ-lyon1.fr                           ##
	##                wassim.g.youssef@gmail.com                                ##
	## ------------------------------------------------------------------------ ##
	##  Version       01 - 2013/11/06                                           ##
	## ------------------------------------------------------------------------ ##

	## ------------------------------------------------------------------------ ##
	##  Definition of the functions                                             ##
	## ------------------------------------------------------------------------ ##

	## ------------------------------------------------------------------------ ##
	##  Simualtion des décès                                                    ##
	## ------------------------------------------------------------------------ ##

	.SimDxt = function(q, e, x1, x2, t1, nsim){
		DxtSim <- vector("list", nsim)
		for(i in 1:nsim){
			DxtSim[[i]] <- matrix(rpois(length(x1) * length(t1), e[x1+1-min(x2),as.character(t1)] * q[x1+1-min(as.numeric(rownames(q))),as.character(t1)]), length(x1), length(t1))
			colnames(DxtSim[[i]]) <- t1; rownames(DxtSim[[i]]) <- x1
		}
		return(DxtSim)
	}

	## ------------------------------------------------------------------------ ##
	##  Coefficient of variation                                                ##
	## ------------------------------------------------------------------------ ##

	.GetCV = function(q, Age, t1, t2, nsim){
		CvMat <- matrix(,length(Age), length(min(t1):max(t2)))
		q2 <- vector("list",length(min(t1):max(t2)))
		for (y in 1:length(min(t1):max(t2))){
			q2[[y]] <- matrix(,length(Age),nsim)
			for (k in 1:nsim){
				for (x in 1:length(Age)){
					q2[[y]][x,k] <- q[[k]][x,y] 
				}
			}
		}
		for(y in 1:length(min(t1):max(t2))){
			sq <- matrix(,length(Age),nsim)
			sum_sim <- apply(q2[[y]],1,sum)
			for(x in 1:length(Age)){
				for(k in 1:nsim){
					sq[x,k] <- (q2[[y]][x,k] - ((sum_sim/nsim)[x]))^2 
				} 
			}
			CvMat[,y] <- sqrt((1/(nsim-1)) * apply(sq,1,sum)) / (sum_sim/nsim)
		}
		return(CvMat)
	}

	## ------------------------------------------------------------------------ ##
	##  Simulated periodic life expectancies                                    ##
	## ------------------------------------------------------------------------ ##

	.GetSimExp = function(q, Age, t1, t2, nsim){
		LifeExp <- vector("list", nsim)
		for(k in 1:nsim){
			LifeExp[[k]] <- matrix(,length(Age),length(min(t1):max(t2)))
			colnames(LifeExp[[k]]) <- as.character(min(t1):max(t2))
			for (j in 1:(length(min(t1):max(t2)))){
				for (x in 1:length(Age)){
					age.vec <- ((Age[x]):max(Age))+1-min(Age)
					LifeExp[[k]][x,j] <- sum(cumprod(1-q[[k]][age.vec, j])) } } }
					LifeExp2 <-vector("list",length(min(t1):max(t2)))
					for (y in 1:(length(min(t1):max(t2)))){
						LifeExp2[[y]] <- matrix(,length(Age),nsim)
						for (k in 1:nsim){
							for (x in 1:length(Age)){
								LifeExp2[[y]][x,k] <- LifeExp[[k]][x,y] } } }
								return(LifeExp2)
							}

	## ------------------------------------------------------------------------ ##
	##  Computation of the quantiles                                            ##
	## ------------------------------------------------------------------------ ##

	.GetQtiles = function(LifeExp, Age, t1, t2, qval){
		Qtile <- vector("list",length(min(t1):max(t2)))
		for (y in 1:(length(min(t1):max(t2)))){
			Qtile[[y]] <- matrix(,length(Age),length(qval))
			colnames(Qtile[[y]]) <- as.character(qval)
			for (i in 1:(length(Age))) {
				Qtile[[y]][i,] <- quantile(LifeExp[[y]][i,],  probs = qval/100) } }
				return(Qtile)
			}

	## ------------------------------------------------------------------------ ##
	##  Computation of the relative dispersion                                  ##
	## ------------------------------------------------------------------------ ##

	.GetRelDisp = function(LifeExp, Age, t1, t2, qval){
		Qtile <- .GetQtiles(LifeExp, Age, t1, t2, qval)
		RelDisp <- matrix(,length(Age),length(min(t1):max(t2)))
		colnames(RelDisp) <- as.character(min(t1):max(t2))
		for (yy in 1:length(min(t1):max(t2))){
			for (xx in 1:length(Age)){
				qvec <- Qtile[[yy]][xx,]		
				RelDisp[xx,yy] <- (qvec[3] - qvec[1]) / qvec[2] 
			}
		}
		RelDisp[RelDisp == Inf] <- NaN; RelDisp[RelDisp == -Inf] <- NaN
		return(RelDisp)
	}

	## ------------------------------------------------------------------------ ##
	##  .GetFitSim function                                                      ##
	## ------------------------------------------------------------------------ ##


	.GetFitSim = function(DxtSim, MyData, AgeMethod, NamesMyData, NameMethod, NbSim, CompletionTable, AgeRangeOpt=NULL,BegAgeComp=NULL,P.Opt=NULL,h.Opt=NULL){
		QxtFittedSim <- QxtFinalSim <- vector("list", NbSim)
		for(i in 1:NbSim){
	## ---------- If Method 1
	if(NameMethod == "Method1"){
		QxtFittedSim[[i]] <- FctMethod1(DxtSim[[i]], MyData$Ext, MyData$QxtRef, AgeMethod, MyData$AgeRef, MyData$YearCom, MyData$YearRef)
	}
	## ---------- If Method 2
	if(NameMethod == "Method2"){
		QxtFittedSim[[i]] <- FctMethod2(DxtSim[[i]], MyData$Ext, MyData$QxtRef, AgeMethod, MyData$AgeRef, MyData$YearCom, MyData$YearRef)
	}
	## ---------- If Method 3
	if(NameMethod == "Method3"){
		QxtFittedSim[[i]] <- FctMethod3(DxtSim[[i]], MyData$Ext, MyData$QxtRef, AgeMethod, MyData$AgeRef, MyData$YearCom, MyData$YearRef)
	}
	## ---------- If Method 4
	if(NameMethod == "Method4"){					
		QxtFittedSim[[i]] <- FctMethod4_2ndPart(DxtSim[[i]], MyData$Ext, MyData$QxtRef, AgeMethod, MyData$AgeRef, MyData$YearCom, MyData$YearRef, P.Opt, h.Opt)
	}	
	## ---------- Completion
	if(CompletionTable == T){
		QxtFinalSim[[i]] <- .CompletionDG2005(QxtFittedSim[[i]]$QxtFitted, as.numeric(rownames(QxtFittedSim[[1]]$QxtFitted)), min(MyData$YearCom):max(MyData$YearRef), AgeRangeOpt, c(BegAgeComp, 130), NameMethod)$QxtFinal
	}
	if(CompletionTable == F){
		QxtFittedSim[[i]]$QxtFitted[QxtFittedSim[[i]]$QxtFitted > 1] <- 1
		QxtFittedSim[[i]]$QxtFitted[is.na(QxtFittedSim[[i]]$QxtFitted)] <- 1
		QxtFinalSim[[i]] <- QxtFittedSim[[i]]$QxtFitted
	}
	}
	return(QxtFinalSim)
	}

	## ------------------------------------------------------------------------ ##
	##  Dispersion function                                                     ##
	## ------------------------------------------------------------------------ ##
									


	Dispersion = function(FinalMethod, MyData, NbSim, CompletionTable=T, Plot = F, Color = MyData$Param$Color){
		AgeMethod=FinalMethod[[1]]$AgeRange
		print("Simulate the number of deaths ...")
		DxtSim <- vector("list",length(MyData)-1)
		names(DxtSim) <- names(MyData)[1:(length(MyData)-1)]
		for (i in 1:(length(MyData)-1)){
			DxtSim[[i]] <- .SimDxt(FinalMethod[[i]]$QxtFinal, MyData[[i]]$Ext, AgeMethod, MyData[[i]]$AgeRef, MyData[[i]]$YearCom, NbSim)
			}
		print("Fit of the simulated data ...")
		QxtFinalSim <- vector("list",length(MyData)-1)
		names(QxtFinalSim) <- names(MyData)[1:(length(MyData)-1)]
		for (i in 1:(length(MyData)-1)){
			QxtFinalSim[[i]] <- .GetFitSim(DxtSim[[i]], MyData[[i]], AgeMethod, names(MyData)[i], FinalMethod[[1]]$NameMethod, NbSim, CompletionTable, FinalMethod[[i]]$AgeRangeOpt, FinalMethod[[i]]$BegAgeComp, FinalMethod[[i]]$P.Opt, FinalMethod[[i]]$h.Opt)
			}
		AgeFinal <- as.numeric(rownames(FinalMethod[[1]]$QxtFinal))
		if(Plot == T){
			Path <- "Results/Graphics/Dispersion"
			.CreateDirectory(paste("/",Path,sep=""))
			print(paste("Create graphics of the fit of the simulated data in .../",Path," ...", sep=""))
			for(j in MyData[[1]]$YearCom){
				if(length(MyData) == 3){
					png(filename=paste(Path,"/",FinalMethod[[1]]$NameMethod,"-GraphSim-",j,".png", sep=""), width  = 3800, height = 2100, res=300, pointsize= 12)
					print(.SimPlot(FinalMethod, QxtFinalSim, MyData, min(AgeFinal):95, j, c(paste("Fit of the simulated data -",names(MyData)[1:(length(MyData)-1)],"- year",j)), NbSim, Color))
					dev.off()
					}
				if(length(MyData) == 2){
					png(filename=paste(Path,"/",FinalMethod[[1]]$NameMethod,"-GraphSim-",j,".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
					print(.SimPlot(FinalMethod, QxtFinalSim, MyData, min(AgeFinal):95, j, paste("Fit of the simulated data -",names(MyData)[i],"- year",j), NbSim, Color))
					dev.off()
					}
				}
			}
		print("Compute the coefficients of variation ...")
		Cv <- vector("list",length(MyData)-1)
		names(Cv) <- names(MyData)[1:(length(MyData)-1)]
		for(i in 1:(length(MyData)-1)){
			Cv[[i]] <- .GetCV(QxtFinalSim[[i]], AgeFinal, MyData[[i]]$YearCom, MyData[[i]]$YearRef, NbSim)
			}
		if(Plot == T){
			print(paste("Create graphics of the coefficients of variation in .../",Path," ...", sep=""))
			for(i in 1:(length(MyData)-1)){
				png(filename=paste(Path,"/",FinalMethod[[1]]$NameMethod,"-CoefVar-",names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
				print(SurfacePlot(as.matrix(Cv[[i]]),expression(cv[xt]), paste("Coefficient of variation, ",FinalMethod[[1]]$NameMethod,", ",names(MyData)[i]," pop.", sep=""), c(min(AgeFinal),130,min(MyData[[i]]$YearCom),max(MyData[[i]]$YearRef)),Color))
				dev.off()
				}
			}
		print("Compute the simulated periodic life expectancies ...")
		LifeExpSim <- vector("list",length(MyData)-1)
		names(LifeExpSim) <- names(MyData)[1:(length(MyData)-1)]
		for(i in 1:(length(MyData)-1)){
			LifeExpSim[[i]] <- .GetSimExp(QxtFinalSim[[i]], AgeFinal, MyData[[i]]$YearCom, MyData[[i]]$YearRef, NbSim)
			}
		print("Compute the quantiles of the simulated periodic life expectancies ...")
		QtilesVal <- c(2.5, 50, 97.5)
		Qtile <- vector("list",length(MyData)-1)
		names(Qtile) <- names(MyData)[1:(length(MyData)-1)]
		for(i in 1:(length(MyData)-1)){	
			Qtile[[i]] <- .GetQtiles(LifeExpSim[[i]], AgeFinal, MyData[[i]]$YearCom, MyData[[i]]$YearRef, QtilesVal)
			}
		if(Plot == T){
			print(paste("Create graphics of the quantiles of the simulated periodic life expectancies in .../",Path," ...", sep=""))
			for(i in 1:(length(MyData)-1)){
				png(filename=paste(Path,"/",FinalMethod[[1]]$NameMethod,"-QtileLifeExp-",names(MyData)[i],".png", sep=""), width  = 8400, height = 1200, res=300, pointsize= 12)
				print(.PlotExpQtle(Qtile[[i]], min(MyData[[i]]$YearCom):max(MyData[[i]]$YearRef) ,(ceiling(min(AgeFinal)/10)*10):100, paste("Quantiles of the simulated periodic life expectancies, ",FinalMethod[[1]]$NameMethod,", ",names(MyData)[i]," pop.",sep=""), Color)) 
				dev.off()
				}
			}
		print("Compute the relative dispersion ...")
		QtilesVal <- c(5, 50, 95)
		RelDisp <- vector("list",length(MyData)-1)
		names(RelDisp) <- names(MyData)[1:(length(MyData)-1)]
		for(i in 1:(length(MyData)-1)){			
			RelDisp[[i]] <- .GetRelDisp(LifeExpSim[[i]], AgeFinal, MyData[[i]]$YearCom, MyData[[i]]$YearRef, QtilesVal)
			}
		if(Plot == T){
			print(paste("Create graphics of the relative dispersion in .../",Path," ...", sep=""))
			for(i in 1:(length(MyData)-1)){
				png(filename=paste(Path,"/",FinalMethod[[1]]$NameMethod,"-RelDisp-",names(MyData)[i],".png", sep=""), width  = 8400, height = 1200, res=300, pointsize= 12)
				print(.PlotRelDisp(RelDisp[[i]], min(MyData[[i]]$YearCom):max(MyData[[i]]$YearRef) , (ceiling(min(AgeFinal)/10)*10):100, paste("Relative dispersion, ",FinalMethod[[1]]$NameMethod,", ",names(MyData)[i]," pop.",sep=""), Color, c(670,50))) 
				dev.off()
				}
			}			
		return(list(Cv = Cv, LifeExpSim = LifeExpSim, Qtile = Qtile, RelDisp = RelDisp))
		}
