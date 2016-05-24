
# file FederovOpt.R
# copyright (C) 2002-2004 by Robert E. Wheeler
#


"optFederov" <-
function (frml,data=sys.frame(sys.parent()),nTrials,center=FALSE,approximate=FALSE,criterion="D",
	evaluateI=FALSE,space=NULL,augment=FALSE,rows=NULL,nullify=0,maxIteration = 100,nRepeats=5,
	DFrac=1,CFrac=1,args=FALSE)
{
	if (!exists(".Random.seed"))
		set.seed(555111666)
	seed<-.Random.seed

	if (missing(frml) || !inherits(frml,c("formula","character"))) {
		if (missing(data))
			stop("frml and data cannot both be missing.")
		frml<-~.
	}

	if (missing(data)) { # Create a data matrix from the global variables in frml
		frmla<-formula(paste("~-1+",paste(all.vars(frml),sep="",collapse="+"),sep=""))
		data<-data.frame(model.matrix(frmla,data))
	}
	else {
		if (!inherits(data,"data.frame")) {
			data<-data.frame(data)   # to insure the columns are named
			if (ncol(data)==1)
				colnames(data)<-"X1"
		}
	}

	if (any(is.na(data)))
		stop("Missing values are not allowed.")

	numericColumn<-sapply(data,is.numeric)

	if (center) {
		means<-apply(data[,numericColumn],2,mean)
		data[,numericColumn]<-sweep(data[,numericColumn,drop=FALSE],2,means)
	}

	frml<-expand.formula(frml,colnames(data),numerics=numericColumn)

	X<-model.matrix(frml,data)

    N <- nrow(X)
    k <- ncol(X)

        doSpace<-TRUE
	if ((criterion=="I" || evaluateI) && !missing(space)) {
		S<-model.matrix(frml,space)
		B<-t(S)%*%S/nrow(S)
		B<-B[lower.tri(B,diag=TRUE)]
		doSpace<-TRUE
	}
	else {
		B<-NULL
		doSpace<-FALSE
	}

	nRound<-0

	if (missing(nTrials) || nTrials==0) {
		if (!augment && !missing(rows)) {
			nTrials<-length(unique(rows))
		}
		else
			nTrials<-k+5
	}
	else {
		if (approximate) {
			nRound<-nTrials # final approximate design will be rounded to this value
			nTrials<-k+5  # the extra 5 is not really needed
		}
	}

	if (nTrials<k)
		stop("nTrials must be greater than or equal to the number of columns in expanded X")

	if (approximate==FALSE && nTrials>N)
		stop("nTrials must not be greater than the number of rows in data")

	nTrials<-as.integer(nTrials) # to be safe


	if (approximate==FALSE) {
		proportions<-NULL
	}
	else {
		proportions<-rep(0,N)
		if (maxIteration<=100)
			maxIteration<-1000
	}

    if (criterion=="D") crit<-0
    else
      if (criterion=="A") crit<-1
    else
      if (criterion=="I") crit<-2
    else
      stop("Did not recognize the criterion.")


        RandomStart<-TRUE
	if (missing(rows) || is.null(rows) || augment) {
		if (augment) {
			rows<-unique(rows)
			augment<-length(rows)
			if (augment>nTrials || augment>N)
				stop("Number of non-duplicated elements in rows is too large.")
			rows<-c(rows-1,rep(0,nTrials-augment))
			RandomStart<-FALSE
		}
		else {
			rows<-rep(0,nTrials)
			if (nullify!=0)
				RandomStart<-FALSE # this has no effect when approximate!=FALSE since nullify is
								   # always used in the C program
			else
				RandomStart<-TRUE
		}
	}
	else {
		RandomStart<-FALSE;
		rows<-unique(rows)
		rows<-rows-1
		if (length(rows)!=nTrials)
			stop("Number of non-duplicated elements in rows must equal nTrials")
	}


	value<-.Call("FederovOpt", X,as.integer(RandomStart),as.integer(rows),as.integer(nullify),
		as.integer(crit),as.integer(evaluateI),as.integer(doSpace),B,as.integer(augment),as.integer(approximate),
		as.double(proportions),as.integer(nTrials),as.integer(maxIteration),as.integer(nRepeats),
		as.double(DFrac),as.double(CFrac),PACKAGE="AlgDesign")

	if (value$error==12) {
		stop("Singular design.")		}
	if (value$error==13) {
		stop("Nullifcation failed to find a starting design.")
		}
	if (value$error!=0)
		stop(value$error)

        proportions<-0
	if (approximate) {
		proportions<-value$proportions;

		if (nRound>0) {
			RowNos<-(1:N)[proportions>(1/(2*maxIteration))]
			if (length(RowNos)>nRound) { # in case nRound is less than the number of support points.
				RowNos<-sort(proportions,index=TRUE,decreasing=TRUE)$ix
				RowNos<-RowNos[1:nRound]
				RowNos<-sort(RowNos)
			}
			proportions<-proportions[RowNos]
			proportions<-proportions/sum(proportions)
			proportions<-efficient.rounding(proportions,nRound)
			RowNos<-RowNos[proportions>0]
			proportions<-proportions[proportions>0]
		}
		else {
			RowNos<-(1:N)[proportions>0]
			proportions<-proportions[RowNos]
			proportions<-round(proportions,4)
		}

	}

	if (approximate==FALSE) {
		RowNos<-sort(1+((value$rows[1:nTrials])%%N))
	}

	if (center) { # reverse centering for ouput
		data[,numericColumn]<-sweep(data[,numericColumn,drop=FALSE],2,-means)
	}

	Design<-data[RowNos,,drop=FALSE]
	if (approximate) {
		Design<-cbind(proportions,Design)
		if (nRound>0)
			colnames(Design)[1]<-"Rep.."
		else
			colnames(Design)[1]<-"Proportion"
	}
	rownames(Design)<-RowNos

	De<-exp(1-1/value$G)


	output<-list(D = value$D, A = value$A)
	if (criterion=="I" || evaluateI) {
		output<-c(output,list( I = value$I))
	}
	output<-c(output,list(Ge = round(value$G,3),Dea = round(De,3),design=Design,rows=RowNos))

	if (args)
		output<-c(output,list(seed=seed),args=actuals(formals("optFederov")))

	output
}

