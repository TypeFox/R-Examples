#############################################
# arf3DS4 S4 VAR/COVAR FUNCTIONS	  		#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#makeDerivs			[user]	
#makeResiduals
#makeWeightsBin
#varcov				[user]
#BIC				[user]
#RMSEA				[user]
#detSigmaDeriv
#wald				[user]
#mcpCorrect			[user]
#readDerivs			[user]
#approxHessianR		
#getlocsquare
#makeBWlocations		


makeDerivs <- 
function(arfmodel,method=c('viaR','direct'))
#make derivs calculates analytical first order derivatives based on the estimated model parameters
{
	
	method=match.arg(method[1],c('viaR','direct'))
	
	if(.model.valid(arfmodel)) {
		#load headinfo for dims
		headinf <- readHeader(getFileInfo(.model.avgdatfile(arfmodel)))
	
		#set the filename
		fn <- paste(.model.modeldatapath(arfmodel),.Platform$file.sep,.model.derivativeFile(arfmodel),sep='')
		
				
		#calculate and write the derivatives
		if(method=='viaR') {
			derivs = .C('dfgauss',as.integer(.model.regions(arfmodel)*.model.params(arfmodel)),as.integer(.model.mask(arfmodel)),as.integer(.nifti.header.dims(headinf)[2]),as.integer(.nifti.header.dims(headinf)[3]),as.integer(.nifti.header.dims(headinf)[4]),as.double(.model.estimates(arfmodel)),as.double(numeric(.model.regions(arfmodel)*.model.params(arfmodel)*length(which(.model.mask(arfmodel)!=0)))))[[7]]
 			con = file(fn,'wb')
			writeBin(derivs,con,double(),endian=.Platform$endian)
			close(con)
		}
		
		if(method=='direct') {
			invisible(.C('dfgaussFile',as.integer(.model.regions(arfmodel)*.model.params(arfmodel)),as.integer(.model.mask(arfmodel)),as.integer(.nifti.header.dims(headinf)[2]),as.integer(.nifti.header.dims(headinf)[3]),as.integer(.nifti.header.dims(headinf)[4]),as.double(.model.estimates(arfmodel)),as.character(fn)))
		}
		
		return(invisible(TRUE))
	
	} else return(invisible(FALSE))
		
}


makeResiduals <- 
function(arfmodel)
#make residuals makes residuals for each run
{
	if(.model.valid(arfmodel)) {
		#set separator
		sp <- .Platform$file.sep
		
		#load in the model data
		model <- .fmri.data.datavec(readData(paste(.model.modeldatapath(arfmodel),sp,.model.fullmodelDataFile(arfmodel),sep='')))
		
		#open a binary connection to the residualfile
		con <- file(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep=''),'wb') 
	
		nonbrain = which(.model.mask(arfmodel)==0)
		if(length(nonbrain)>0) model = model[-nonbrain]
		
		#write the residuals per run
		for(bfile in .model.betafiles(arfmodel)) {
			
			data = .fmri.data.datavec(readData(bfile))
			if(length(nonbrain)>0) data = data[-nonbrain]
			writeBin(data-model,con,double(),endian=.Platform$endian)
	
		}
		#close the connection
		close(con)
	
		return(invisible(TRUE))
		
	} else return(invisible(FALSE))

}

makeWeightsBin <-
function(arfmodel) 
#make binary weightsfile
{

	weights = readData(.model.avgWfile(arfmodel))
	weightdata = .fmri.data.datavec(weights)
	
	rem = which(.model.mask(arfmodel)==0)
	if(length(rem)>0) weightdata = weightdata[-rem]
	
	con <- file(paste(.model.modeldatapath(arfmodel),.Platform$file.sep,.model.weightFile(arfmodel),sep=''),'wb')
	writeBin(weightdata,con,double(),endian=.Platform$endian)
	close(con)
	
}


varcov <- 
function(arfmodel)
#calculate sandwich estimate
{
	#set separator
	sp <- .Platform$file.sep
	
	if(.model.valid(arfmodel)) {
	
		#check for deriv and resid
		if(!file.exists(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep=''))) makeResiduals(arfmodel) 
		if(!file.exists(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep=''))) makeDerivs(arfmodel)
		if(!file.exists(paste(.model.modeldatapath(arfmodel),sp,.model.weightFile(arfmodel),sep=''))) makeWeightsBin(arfmodel)
		
		st_time <- Sys.time()
		
		#try to invert hessian
		hessian <- try(solve(.model.hessian(arfmodel)),silent=T)
				
		#if Hessian is not ok, try pseudoinverse
		if(!is.null(attr(hessian,'class')))  {
			#library(corpcor)
			hessian <- try(pseudoinverse(.model.hessian(arfmodel)),silent=T)
			.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[varcov] Using Pseudoinverse (hessian probably singular)',sep=''))
			pseudo = T
		} else { pseudo = F}
				
		#check if hessian is good
		if(is.null(attr(hessian,'class')))  {
			weights <- readData(.model.avgWfile(arfmodel))
			n = .model.n(arfmodel)
			
			#perform the inner_sandwich procedure
			if(.model.sandwichmethod(arfmodel)[1]=='diag') B <- try(.C('innerSWdiag',as.integer(n),as.integer(.model.regions(arfmodel)*.model.params(arfmodel)),as.integer(.model.runs(arfmodel)),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.weightFile(arfmodel),sep='')),as.double(numeric((.model.regions(arfmodel)*.model.params(arfmodel))^2)))[[7]],silent=T)
			if(.model.sandwichmethod(arfmodel)[1]=='fast') B <- try(.C('innerSWfast',as.integer(n),as.integer(.model.regions(arfmodel)*.model.params(arfmodel)),as.integer(.model.runs(arfmodel)),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.weightFile(arfmodel),sep='')),as.double(numeric((.model.regions(arfmodel)*.model.params(arfmodel))^2)))[[7]],silent=T)
			if(.model.sandwichmethod(arfmodel)[1]=='full') B <- try(.C('innerSWfull',as.integer(n),as.integer(.model.regions(arfmodel)*.model.params(arfmodel)),as.integer(.model.runs(arfmodel)),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.weightFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,'mean',.model.residualFile(arfmodel),sep='')),as.double(numeric((.model.regions(arfmodel)*.model.params(arfmodel))^2)))[[8]],silent=T)
			if(.model.sandwichmethod(arfmodel)[1]=='bw') {
				#save original model values
				old_n = .model.n(arfmodel)
				old_mask = .model.mask(arfmodel)
				
				#set mask and n to full
				.model.n(arfmodel) <- length(.fmri.data.datavec(weights))
				.model.mask(arfmodel) <- rep(1,.model.n(arfmodel))
				
				#remake Resids,Derivs,Weights
				makeResiduals(arfmodel)
				makeDerivs(arfmodel)
				makeWeightsBin(arfmodel)
				
				#make bandwidth and locations matrices
				bw = as.numeric(.model.sandwichmethod(arfmodel)[2])
				escapevar = .model.n(arfmodel)+1
				Lv = makeBWLocations(arfmodel,escapevar)	
				
				#run BandWidth SW
				B <- try(.C('innerSWbw',as.integer(.model.n(arfmodel)),as.integer(.model.regions(arfmodel)*.model.params(arfmodel)),as.integer(.model.runs(arfmodel)),as.integer(((bw*2)+1)^3),as.integer(escapevar),as.integer(Lv),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.weightFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,'mean',.model.residualFile(arfmodel),sep='')),as.double(numeric((.model.regions(arfmodel)*.model.params(arfmodel))^2)))[[11]],silent=T)
				
				#reset original values
				.model.n(arfmodel) = old_n
				.model.mask(arfmodel) = old_mask
			}
			
			if(.model.sandwichmethod(arfmodel)[1]=='bwf') {
				#save original model values
				old_n = .model.n(arfmodel)
				old_mask = .model.mask(arfmodel)
				
				#set mask and n to full
				.model.n(arfmodel) <- length(.fmri.data.datavec(weights))
				.model.mask(arfmodel) <- rep(1,.model.n(arfmodel))
				
				#remake Resids,Derivs,Weights
				makeResiduals(arfmodel)
				makeDerivs(arfmodel)
				makeWeightsBin(arfmodel)
			
				#make bandwidth and locations matrices
				bw = as.numeric(.model.sandwichmethod(arfmodel)[2])
				escapevar = .model.n(arfmodel)+1
				Lv = makeBWLocations(arfmodel,escapevar)	
				
				#run BandWidth SW
				B <- try(.C('innerSWbwfast',as.integer(.model.n(arfmodel)),as.integer(.model.regions(arfmodel)*.model.params(arfmodel)),as.integer(.model.runs(arfmodel)),as.integer(((bw*2)+1)^3),as.integer(escapevar),as.integer(Lv),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,.model.weightFile(arfmodel),sep='')),as.character(paste(.model.modeldatapath(arfmodel),sp,'mean',.model.residualFile(arfmodel),sep='')),as.double(numeric((.model.regions(arfmodel)*.model.params(arfmodel))^2)))[[11]],silent=T)
				
				#reset original values
				.model.n(arfmodel) = old_n
				.model.mask(arfmodel) = old_mask
			}
			
			#check if innersandwich works
			if(is.null(attr(B,'class'))) {
				#set B to be a matrix
				dim(B) <- c(.model.regions(arfmodel)*.model.params(arfmodel),.model.regions(arfmodel)*.model.params(arfmodel))
				
				#add the outer sandwich parts
				if(!pseudo) SW <- try(solve(.5*.model.hessian(arfmodel))%*%B%*%solve(.5*.model.hessian(arfmodel)),silent=T)
				if(pseudo) SW <- try(pseudoinverse(.5*.model.hessian(arfmodel))%*%B%*%pseudoinverse(.5*.model.hessian(arfmodel)),silent=T)
				
				
				#check if outersandwich works
				if(is.null(attr(SW,'class'))) {
					.model.varcov(arfmodel) <- SW
					.model.valid(arfmodel) <- TRUE
				} else { #outersandiwch not good
					.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[varcov] outerSandwich did not compute: ',gsub('\n','',SW),sep='')) 
					.model.valid(arfmodel) <- FALSE
				}
			} else { #innersandiwch not good
				.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[varcov] innerSandwich did not compute: ',gsub('\n','',B),sep='')) 
				.model.valid(arfmodel) <- FALSE
			}
		} else { #hessian not good
			.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[varcov] Hessian is singular: ',gsub('\n','',hessian),sep='')) 
			.model.valid(arfmodel) <- FALSE
		}
	
		en_time <- Sys.time()
		.model.proctime(arfmodel)[1,2] <- as.numeric(difftime(en_time,st_time,units='sec'))
		
	} else { #model not valid
		.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'[varcov] No valid model. var/cov not calculated')
		.model.valid(arfmodel) <- FALSE
	}
	
	
	#remove derivatives and residuals
	#fn <- paste(.model.modeldatapath(arfmodel),.Platform$file.sep,.model.derivativeFile(arfmodel),sep='')
	#if(file.exists(fn)) file.remove(fn)	
	#fn <- paste(.model.modeldatapath(arfmodel),.Platform$file.sep,.model.derivativeFile(arfmodel),sep='')
	#if(file.exists(fn)) file.remove(fn)
	#fn <- paste(.model.modeldatapath(arfmodel),.Platform$file.sep,.model.residualFile(arfmodel),sep='')
	#if(file.exists(fn)) file.remove(fn)
	#fn <- paste(.model.modeldatapath(arfmodel),.Platform$file.sep,'mean',.model.residualFile(arfmodel),sep='')
	#if(file.exists(fn)) file.remove(fn)

	#save the modelInfo
	saveModel(arfmodel)
	
	#return the varcov
	return(invisible(arfmodel))
}


BIC <- 
function(arfmodel,options=loadOptions(arfmodel)) 
#calculate the BIC
{
	if(.model.valid(arfmodel)) {
		
		#read in weights, used in constant
		Wdata <- readData(.model.avgWfile(arfmodel))
		n <- .model.n(arfmodel) 
		W <- .fmri.data.datavec(Wdata)[1:(.fmri.data.dims(Wdata)[2]*.fmri.data.dims(Wdata)[3]*.fmri.data.dims(Wdata)[4])]
	
		#calculate the determinant of the weights
		dtm <- prod(W)
		
		#check if determinant is valid
		if(!is.na(dtm) & !is.nan(dtm)) {
			if(is.numeric(try(log(n))) & is.numeric(try(log(dtm))) & is.numeric(log(.model.minimum(arfmodel)))) {
				cons <- try((2*(((n/2)*log(2*pi))+((1/2)*log(dtm))+((1/2)*(.model.minimum(arfmodel))))),silent=T)
				
				if(cons!=0) { 
					if(log(dtm)==-Inf) {
						dtm=1e-323
						#.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'Determinant (in BIC) was -Inf, set to minimum value 1e-323.')
					}
					if(log(dtm)==Inf) {
						dtm=1e308
						#.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'Determinant (in BIC) was Inf, set to maximum value 1e308.')
					}
					cons <- (2*(((n/2)*log(2*pi))+((1/2)*log(dtm))+((1/2)*(.model.minimum(arfmodel)))))
				}
				
			} else { #logs not valid
				.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'Error calculating logs. BIC not calculated')
				.model.valid(arfmodel) <- FALSE
			}
		} else { #determinant of W is not valid
			.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'Invalid determinant. BIC not calculated')
			.model.valid(arfmodel) <- FALSE
		} 
		
		#check if constant is a number and calculate BIC
		if(is.numeric(cons)) {
			.model.fit(arfmodel)[1,1]  <- cons + (((.model.regions(arfmodel)*.model.params(arfmodel)))*log(n))
			
		} else { #constant is not a number
			.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'Constant invalid. BIC not calculated')
			.model.valid(arfmodel) <- FALSE
		} 
		
	} else	{
		.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'No valid model. BIC not calculated')
		.model.valid(arfmodel) <- FALSE
	}
	
	#save the modelInfo
	saveModel(arfmodel)
	
	return(invisible(arfmodel))
	
}

RMSEA <- 
function(arfmodel,options=loadOptions(arfmodel)) 
#RMSEA calculates root mean square errors of models
{
	
	#check model validity
	if(.model.valid(arfmodel)) {
		
		#set number of voxels
		n <- .model.n(arfmodel) 
			
		#Hotellings T
		HTs = .model.runs(arfmodel)*.model.minimum(arfmodel)
		
		#noncentrality parameter	
		ncp = max(c((HTs-((n-(.model.regions(arfmodel)*.model.params(arfmodel)))/.model.runs(arfmodel))),0))
		
		#check if ncp < 0
		if(ncp<0) .model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'Noncentrality parameter smaller than zero, ncp is set to zero')
		
		#calculate RMSEA
		eps = sqrt(ncp/(n-(.model.regions(arfmodel)*.model.params(arfmodel))))

		if(checkVersion(.model.version(arfmodel),1,2,7)) .model.fit(arfmodel)[1,2] = eps else .model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'RMSEA not calculated for models prior to v1.2-7')
		
	} else {
		#if not good warn
		.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'No valid model. RMSEA not calculated')
		.model.valid(arfmodel) <- FALSE
	}
	
	#save the modelInfo
	saveModel(arfmodel)
	
	return(invisible(arfmodel))
	
}

 
detSigmaDeriv <- 
function(theta) 
#get determinant of sigma with derivatives
{
	
	#set relevant parameters
	theta3=theta[1]
	theta4=theta[2]
	theta5=theta[3]
	theta6=theta[4]
	theta7=theta[5]
	theta8=theta[6]
	
	#expression for det(sigma)
	det_sigma <- expression((theta3^2)*(theta4^2)*(theta5^2)-(theta3^2)*(theta8*theta4*theta5)*(theta8*theta4*theta5)-(theta6*theta4*theta3)*(theta6*theta4*theta3)*(theta5^2)+(theta6*theta4*theta3)*(theta7*theta5*theta3)*(theta8*theta4*theta5)+(theta7*theta5*theta3)*(theta6*theta4*theta3)*(theta8*theta4*theta5)-(theta7*theta5*theta3)*(theta4^2)*(theta7*theta5*theta3))			

	#return derivative and value
	return(list(value=eval(det_sigma),gradient=attr(eval(deriv(det_sigma,c('theta3','theta4','theta5','theta6','theta7','theta8'))),'gradient')))
	
}

wald <- 
function(arfmodel,waldobject=new('wald'),options=loadOptions(arfmodel)) 
#calculate Wald statistics 
{

	if(.model.modeltype(arfmodel)!='gauss') stop('wald statistics can only be calculated for the full model')
	
	if(length(.model.varcov(arfmodel))==0) stop('(co)variance matrix not yet calculated, cannot compute Wald statistics!')
	
	#define function to calculate Wald statistic 
	W <- function(a,A,C) t(a)%*%solve(A%*%C%*%t(A))%*%a
		
	if(.model.valid(arfmodel)) {
		
		#if no design matrix is specified in the waldobject, make the default matrix (zero-filled) 
		if(dim(.wald.consts(waldobject))[1]==0) .wald.consts(waldobject) <- matrix(0,.model.regions(arfmodel),5)
		
		# get dimensions (set number of voxels)
		n <- .model.n(arfmodel) 
		
		#set relevant matrix sizes and dfs
		.wald.stats(waldobject) <- matrix(0,.model.regions(arfmodel),5)
		.wald.pvalues(waldobject) <- matrix(0,.model.regions(arfmodel),5)
		.wald.df1(waldobject) <- rep(n,5)
		.wald.df2(waldobject) <- .wald.df1(waldobject)-rep(.model.regions(arfmodel)*.model.params(arfmodel),5)
		
		#perform hypothesis tests for each region and for locations, extent and amplitude
		for(region in 1:.model.regions(arfmodel)) {
			
			#select the 10*10 vcov matrix and estimates for each region, with determinant and deriv
			theta <- .model.estimates(arfmodel)[((1+(region-1)*.model.params(arfmodel)):(region*.model.params(arfmodel)))]
			C <- .model.varcov(arfmodel)[((1+(region-1)*.model.params(arfmodel)):(region*.model.params(arfmodel))),((1+(region-1)*.model.params(arfmodel)):(region*.model.params(arfmodel)))]
			sigma <- detSigmaDeriv(theta[4:9])
			
			#define the a matrix (containing hypotheses), uses info from the designmatrix
			a <- c(theta[1]-.wald.consts(waldobject)[region,1],theta[2]-.wald.consts(waldobject)[region,2],theta[3]-.wald.consts(waldobject)[region,3],sigma$value-.wald.consts(waldobject)[region,4],theta[10]-.wald.consts(waldobject)[region,5])

			#define the A matrix (containing the derivatives of a)
			A <- matrix(0,5,10)
			A[1,1] <- 1
			A[2,2] <- 1
			A[3,3] <- 1
			A[4,4] <- sigma$gradient[1]
			A[4,5] <- sigma$gradient[2]
			A[4,6] <- sigma$gradient[3]
			A[4,7] <- sigma$gradient[4]
			A[4,8] <- sigma$gradient[5]
			A[4,9] <- sigma$gradient[6]
			A[5,10] <- 1
			
			#perform tests for locations (calc stats and p-values)
			.wald.stats(waldobject)[region,1] <- W(a[1],A[1,1],C[1,1])
			.wald.pvalues(waldobject)[region,1] <- 1-pf(.wald.stats(waldobject)[region,1],.wald.df1(waldobject)[1],.wald.df2(waldobject)[1])
			.wald.stats(waldobject)[region,2] <- W(a[2],A[2,2],C[2,2])
			.wald.pvalues(waldobject)[region,2] <- 1-pf(.wald.stats(waldobject)[region,2],.wald.df1(waldobject)[2],.wald.df2(waldobject)[2])
			.wald.stats(waldobject)[region,3] <- W(a[3],A[3,3],C[3,3])
			.wald.pvalues(waldobject)[region,3] <- 1-pf(.wald.stats(waldobject)[region,3],.wald.df1(waldobject)[3],.wald.df2(waldobject)[3])
			
			#perform tests for spatial extent (calc stats and p-values)
			.wald.stats(waldobject)[region,4] <- W(a[4],matrix(A[4,(4:9)],1,6),C[(4:9),(4:9)])
			.wald.pvalues(waldobject)[region,4] <- 1-pf(.wald.stats(waldobject)[region,4],.wald.df1(waldobject)[4],.wald.df2(waldobject)[4])
			
			#perform tests for amplitude (calc stats and p-values)
			.wald.stats(waldobject)[region,5] <- W(a[5],A[5,10],C[10,10])
			.wald.pvalues(waldobject)[region,5] <- 1-pf(.wald.stats(waldobject)[region,5],.wald.df1(waldobject)[5],.wald.df2(waldobject)[5])
		}	
		
		.model.wald(arfmodel) <- waldobject	
		
		#save the model file
		#saveModel(arfmodel)
		
	} else	warning('No valid model. wald statistics not calculated.')
	
	#save the modelInfo
	saveModel(arfmodel)
	
	return(invisible(arfmodel))
	
}


mcpCorrect <- 
function(fmridata,type=c('uncorrected','bonferroni','FDR'),alpha=.05,q=.05,cv=1,df=100,sig.steps=1,adj.n=T) 
#calulate the multiple comparison correction on fmri data (uncorrected,bonferroni or FDR)
{
	veclen <- length(.fmri.data.datavec(fmridata))
	if(adj.n) n <- length(.fmri.data.datavec(fmridata)[.fmri.data.datavec(fmridata)!=0]) else n <- length(.fmri.data.datavec(fmridata))
	
	pseq = seq(sig.steps,1,-1)
	
	which <- match.arg(type)
		
	if(which=='bonferroni') {
		adj.p = alpha/n
		pvec = adj.p / seq(1,sig.steps)
		sigvec = numeric(veclen)
		for(i in 1:sig.steps) 	sigvec[pvec[i]>(1-pt(abs(.fmri.data.datavec(fmridata)),df))]=1/pseq[i]
		thres.t = abs(qt(adj.p,df))
		thres.p = adj.p
		
	}
	
	if(which=='FDR') {
		fdr.value = sort(1-pt(.fmri.data.datavec(fmridata),df))
		i = 1
		while(fdr.value[i]<=((i*q)/(n*cv))) {
			i = i + 1
			if(i==length(fdr.value)) break()		
		}
		
		i = max(i-1,i)
		fdr=fdr.value[i]
		pvec = fdr / seq(1,sig.steps)
		sigvec = numeric(veclen)
		for(i in 1:sig.steps) 	sigvec[pvec[i]>(1-pt(abs(.fmri.data.datavec(fmridata)),df))]=1/pseq[i]
		thres.t = abs(qt(fdr,df))
		thres.p = fdr
		
	}
	
	if(which=='uncorrected') {
		adj.p = alpha
		pvec = adj.p / seq(1,sig.steps)
		sigvec = numeric(veclen)
		for(i in 1:sig.steps) 	sigvec[pvec[i]>(1-pt(abs(.fmri.data.datavec(fmridata)),df))]=1/pseq[i]
		thres.t = abs(qt(adj.p,df))
		thres.p = adj.p
		
	}
	
	#make masked and overlay
	masked_fmridata <- overlay_fmridata <- fmridata
	.fmri.data.datavec(overlay_fmridata) <- sigvec
 	.fmri.data.datavec(masked_fmridata)[sigvec==0] = 0 
		
	return(list(masked=masked_fmridata,overlay=overlay_fmridata,threshold=list(t=thres.t,p=thres.p)))
	
}

readDerivs <- 
function(arfmodel) 
#readDerivs reads in a binary derivative file
{
	sp = .Platform$file.sep
	
	fn = paste(.model.modeldatapath(arfmodel),sp,.model.fullmodelDataFile(arfmodel),sep='')
	dat = readData(fn)
	
	n = length(which(.model.mask(arfmodel)!=0))
	rm(dat)
	p = .model.regions(arfmodel)*.model.params(arfmodel)
	fn = paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep='')
	
	if(file.exists(fn)) {
		con <- file(fn,open='rb')		
		dfvec <- readBin(con,double(),n=n*p,size=.Machine$sizeof.longlong,endian=.Platform$endian)
		dim(dfvec) = c(n,p)
		close(con)
		return(dfvec)
	}
	return(NULL)
}


approxHessianR <- 
function(arfmodel) 
#approxHessian calculates the approximate hessian (used only to check C code)
{
	
	df = readDerivs(arfmodel)
	W = .fmri.data.datavec(readData(.model.avgWfile(arfmodel)))
	
	hessian = 2 * (t(df)%*%(df*(1/W)))
	
	return(hessian)
	
}

getlocsquare <-
function(x,y,z,dx,dy,dz,bw,escapevar,mask) 
#return voxel-vector locations from a box with width bw*2 + 1
{
	fw = (bw*2)+1
	lv = numeric(fw^3)
	
	xv = (x-bw):(x+bw)
	yv = (y-bw):(y+bw)
	zv = (z-bw):(z+bw)
	
	i=1;
	for(zz in zv) {
		for(yy in yv) {
			for(xx in xv) {
				if(xx<1 | xx>dx | yy<1 | yy>dy | zz<1 | zz>dz) lv[i]=escapevar else lv[i] = (xx + (yy-1)*(dx) + (zz-1)*(dx)*(dy))-1
				#if(lv[i]!=escapevar) {
					#if(mask[(xx + (yy-1)*(dx) + (zz-1)*(dx)*(dy))]==0) lv[i]=escapevar;
				#}
				i=i+1
			}
		}
	}
	
	#returns location in c-style (starting at zero)
	return(lv)
	
}

makeBWLocations <-
function(arfmodel,escapevar)
#make Bandwidth location matrix
{
	
	bw = as.numeric(.model.sandwichmethod(arfmodel)[2])
	
	dat = readData(.model.avgdatfile(arfmodel))
	dx = .fmri.data.dims(dat)[2]
	dy = .fmri.data.dims(dat)[3]
	dz = .fmri.data.dims(dat)[4]
	
	#make and fill locations matrix
	lv = matrix(NA,length(which(.model.mask(arfmodel)!=0)),((bw*2)+1)^3)
	i=j=1
	for(z in 1:dz) {
		for(y in 1:dy) { 
			for(x in 1:dx) {
				if(.model.mask(arfmodel)[j]!=0) {
					lv[i,] = getlocsquare(x,y,z,dx,dy,dz,bw,escapevar,.model.mask(arfmodel))
					i=i+1
				}
				j=j+1
			}
		}
	}
	
	return(as.vector(lv))
}


 