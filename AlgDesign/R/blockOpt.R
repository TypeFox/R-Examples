# file blockOpt.R
# copyright (C) 2002-2004 by Robert E. Wheeler
#



"optBlock" <-
function (frml,withinData=sys.frame(sys.parent()),blocksizes,rows=NULL,
	wholeBlockData=NULL,center=FALSE,nRepeats=5,criterion="D",args=FALSE)

{

	if (!exists(".Random.seed"))
		set.seed(555111666)

	seed<-.Random.seed

	contr.block<-function (n, contrasts = TRUE) 
	{
		# makes a contrast matrix with all 1's
		# makes the first row all zeros if the first level is zero
		if (length(n) <= 1) {
			if (is.numeric(n) && length(n) == 1 && n > 1) 
				levels <- 1:n
			else stop("contrasts are not defined for 0 degrees of freedom")
		}
		else levels <- n
		lenglev <- length(levels)
		if (contrasts) {
			cont <- array(1, c(lenglev, lenglev - 1), list(levels, 
				NULL))
			if (levels[1]==0)
				cont[1,]<-0
		}
		else {
			cont <- array(0, c(lenglev, lenglev), list(levels, levels))
			cont[col(cont) == row(cont)] <- 1
		}
		cont
	}

	set.factors<-function(A,B,indicator=FALSE) {
		# Turns variables in A into factors corresponding to those in B
		# A and B must have the same number of columns and be data.frames
		# if indicator is true, makes the levels start with 0
		x<-sapply(B,is.factor)
		if (any(x)) {
			n<-ncol(B)
			for (i in (1:n)[x]) {
				m<-length(levels(B[,i]))
				if (indicator)
					levels<-0:(m-1)
				else
					levels<-1:m
				v<-factor(A[,i],levels=levels)
				A[,i]<-C(v,contr.block(levels))
			}

		}
		A
	}

	remove.nesting<-function(frml) {
		# replaces / with * 
		frml<-deparse(frml)
		frml<-gsub("/","*",frml)
		frml<-gsub("%in%","*",frml)
		frml<-as.formula(frml)
		frml
	}


	######### Check the input and make withinData if necessary

	if (missing(frml) || !inherits(frml,c("formula","character"))) {
		if (missing(withinData))
			stop("frml and withinData cannot both be missing.")
		frml<-~.
	}
		
	if (missing(withinData)) { # Create a data matrix from the global variables in frml
		frmla<-formula(paste("~-1+",paste(all.vars(frml),sep="",collapse="+"),sep=""))
		withinData<-data.frame(model.matrix(frmla,withinData))
	}
	else {
		if (!inherits(withinData,"data.frame")) {
			withinData<-data.frame(withinData)   # to insure the columns are named
			if (ncol(withinData)==1)
				colnames(withinData)<-"X1" 
		} 
	} 

	if (any(is.na(withinData)))
		stop("Missing values are not allowed.")

	if (colnames(withinData)[1]=="Rep..") { # Expand data according to Rep.. from optFederov() output.
		ND<-nrow(withinData)		
		reps<-withinData[,1]
		withinData<-withinData[,-1,drop=FALSE]
		withinData<-withinData[rep(1:ND,reps),]
		rownames(withinData)<-1:nrow(withinData)
	}

	########## Center if needed and prepare to expand the formula

	withinBlock<-withinData # to preserve data for output in case of centering

	varNames<-colnames(withinBlock)
	if (!missing(wholeBlockData)) {
		varNames<-c(colnames(wholeBlockData),varNames)
	}


	numericColumns<-sapply(withinBlock,is.numeric) # used to keep expand.formula honest

	if (center) { 
		meansWithin<-apply(withinBlock[,numericColumns,drop=FALSE],2,mean)
		withinBlock[,numericColumns]<-sweep(withinBlock[,numericColumns,drop=FALSE],2,meansWithin)
	}

	if (!missing(wholeBlockData)) {
		if (!inherits(wholeBlockData,"data.frame"))
			wholeBlockData<-data.frame(wholeBlockData)
		if (colnames(wholeBlockData)[1]=="Rep$") { # Expand data according to Rep$ from optFederov() output.
			ND<-nrow(wholeBlockData)		
			reps<-wholeBlockData[,1]
			wholeBlockData<-wholeBlockData[,-1]
			wholeBlockData<-wholeBlockData[rep(1:ND,reps),]
			rownames(wholeBlockData)<-1:nrow(wholeBlockData)
		}
		wholeBlock<-wholeBlockData # to preserve wholeBlockData in case of centering
		numericWhole<-sapply(wholeBlock,is.numeric)
		if (center) {
			meansWhole<-apply(wholeBlock[numericWhole],2,mean)
			wholeBlock[,numericWhole]<-sweep(wholeBlock[,numericWhole,drop=FALSE],2,meansWhole)
		}
		numericColumns<-c(numericWhole,numericColumns)
	}


    if (criterion=="D") crit<-0
    else
      if (criterion=="Dpc") crit<-1
    else
      if (criterion=="Dp") crit<-2
    else
	  if (criterion=="OB") crit<-3
	else
	  if (criterion=="OBS") crit<-4
	else
      stop("Did not recognize the criterion.")

	######## The formula

	frml<-expand.formula(frml,varNames,numerics=numericColumns)
	frml<-remove.nesting(frml)

	#########  Deal with the whole blocks if present 
	# Create matrices (Uwh,Within) and (Whole,Uw), where Uwh, and Uw are appropriately 
	# dimensioned matrices of 1's, in such a fashion that model.matrix() will produce x=xr*delta.
	# (see the documentation). Note, the following does not work with nesting, but
	# remove.nesting(), above, has changed this to multiplication. The nested and the
	# multiplicative models are related by a non-singular linear transformation, with 
	# first column (1,0,...,0)', which ensures that the D criteria for the two will
	# differ by a constant, and thus have no effect on the optimization.

	if (!missing(wholeBlockData)) {
		if (any(is.na(wholeBlock)))
			stop("Missing values are not allowed.")


		NW<-nrow(wholeBlock)
		kW<-ncol(wholeBlock)

		if (NW!=length(blocksizes))
			stop("blocksizes and wholeBlock must indicate the same number of blocks")
		nV<-length(all.vars(frml))
		if (kW>=nV)
			stop("All variables seem to be block variables")
		ND<-nrow(withinBlock)
		kD<-ncol(withinBlock)

		wholeBlockNames<-colnames(wholeBlock)

			# make NW by kD matrix of 1's to append to wholeBlock
		dataW<-data.frame(matrix(1,NW,kD))
		colnames(dataW)<-colnames(withinBlock)
		dataW<-set.factors(dataW,withinBlock)

			# make ND by kW matrix of 1's to preappend to data
		wholeD<-data.frame(matrix(1,ND,kW))
		colnames(wholeD)<-colnames(wholeBlock)
		wholeD<-set.factors(wholeD,wholeBlock)

		Z<-cbind(wholeBlock,dataW)
		Z1<-data.frame(matrix(c(rep(1,kW),rep(0,kD)),1,(kW+kD)))
		colnames(Z1)<-colnames(Z)
		Z1<-set.factors(Z1,Z,TRUE)
		Z<-model.matrix(frml,Z)

		Z1<-model.matrix(frml,Z1) # 1's mark terms involving whole blocks only
								  # Note, Z1 contains 1 for the intercept


			# After expansion, eliminate intercept and any term involving whole blocks only
		Z<-Z[,!Z1,drop=FALSE]


		withinBlock<-cbind(wholeD,withinBlock)

	}
	else {
		Z<-NULL
	}
	
	######### Expand the model in X for within block data and remove any pure whole block terms

	X<-model.matrix(frml,withinBlock)

	if (!missing(wholeBlockData)) { # eliminate intercept and terms involving whole blocks only
		X<-X[,!Z1,drop=FALSE]
	}
	else {  # use this method of removing constant since -1 in the formula will cause
			# problems when model.matrix() codes nesting factors.
		constPosition<-match("(Intercept)",dimnames(X)[[2]],nomatch=0)
		if (constPosition>0)
			X<-X[,-constPosition,drop=FALSE]
	}

    N <- nrow(X)
    k <- ncol(X)
	
	nB<-length(blocksizes)

	nBlock<-sum(blocksizes)

	if (min(blocksizes)>N)
		stop("The number of trials must be at least as large as the minimum blocksize.")

	if (!missing(rows)) {
		if (length(rows)!=nBlock)
			stop("The length of rows must equal the sum of the blocksizes.")
		if (max(rows)>N)
			stop("Some element of rows is larger than the number of withinData rows.")
	}



	if (missing(wholeBlockData)) {
		if (nBlock<(k+nB))
			stop("The number of withinData rows is not large enough to support the blocked model.")
		doWholeBlock<-FALSE; # used only by Call
	}
	else {
		if (nBlock<(1+k+length(colnames(wholeBlockData))))
			stop("The number of withinData rows is not large enough to support the blocked model.")	
		doWholeBlock<-TRUE; # used only by Call
	}

        initRows<-TRUE
	if (!missing(rows)) {
		rows<-rows-1
		initRows<-TRUE  # used only by Call
	}
	else
		initRows<-FALSE # used only by Call


    value <- .Call("BlockOpt", X,as.integer(initRows),as.integer(rows),as.integer(nB),as.integer(blocksizes),
		as.integer(doWholeBlock),Z,as.integer(nRepeats),as.integer(crit),PACKAGE="AlgDesign")

	if (value$error==13) {
		stop("All repeats produced singular designs.")
	}
	if (value$error==27) {
		stop("The Dp algorithimm has produced a singular design.")
	}
	if (value$error==21)
		stop("Overrun maxArray. Please save setup and report this.")
	if (value$error>0 && value$error<13)
		stop("Memory allocation error.")

	BlockArray<-value$BlockArray
	Blocks<-list(1)

	design<-NULL
	fin<-0
	rnames<-NULL
	for (i in 1:nB) {
		bs<-blocksizes[i]
		strt<-fin+1
		fin<-fin+bs
		arr<-sort(BlockArray[strt:fin])
		if (missing(wholeBlockData)) {
			Blocks[i]<-list(withinData[arr,,drop=FALSE])
			design<-rbind(design,withinData[arr,,drop=FALSE])
			rnames<-c(rnames,rownames(withinData[arr,,drop=FALSE]))
		}
		else {
			dta<-cbind(wholeBlockData[rep(i,bs),,drop=FALSE],withinData[arr,,drop=FALSE])
			colnames(dta)<-varNames
			rnames<-c(rnames,rownames(withinData[arr,,drop=FALSE]))
			Blocks[i]<-list(dta)
			design<-rbind(design,dta)
		}
	}
	names(Blocks)<-paste("B",1:nB,sep="")

	output<-list(D=value$D)
	if (crit==1)
		output<-c(output,list(Dpc=value$Dp))
	else
	if (crit==2)
			output<-c(output,list(Dp=value$Dp))
	else
	if (crit==3 || crit==4)
		output<-c(output,list(SS=value$diagonality))
	else
		output<-c(output,list(diagonality=round(value$diagonality,3)))
	output<-c(output,list(Blocks=Blocks,design=design,rows=as.numeric(rnames)))
	if (value$error==22)
		print("No improvement over initial random design.")
	if (args)
		output<-c(output,list(seed=seed),args=actuals(formals("optBlock")))
	output
}

