# file MonteOpt.R
# copyright (C) 2002-2004 by Robert E. Wheeler
#

"optMonteCarlo" <-
function(frml,data,nTrials,approximate=FALSE,criterion="D",evaluateI=FALSE,space=NULL,mixtureSum=1,
	constraints=NULL,RandomStart=TRUE,nRepeats=5,nCand,nCandNull,DFrac=1,CFrac=1,args=FALSE)
{
	if (!exists(".Random.seed"))
		set.seed(555111666)
	seed<-.Random.seed

	RandomLevels<-function(N,s) {
	# returns a N element vector of uniform probabilities on s levels
		s<-s-1
		round(s*runif(N))/s
	}


	RandomMixture<-function(N) {
        # Returns a vector of random Barycentric values.
          v<-rexp(N)
          v<-v/sum(v)
          if (mixtureSum!=0) {
            mixtureSum<-round(mixtureSum,mixRound)
            v=v*mixtureSum
          }
          # Must make sure v sums to the rounded  mixtureSum
          v<-round(v,mixRound)
          d<-mixtureSum-sum(v)
          i<-sample(1:N)[1]
          v[i]<-v[i]+d
          v
	}

	RandomCand<-function (N,nTrials=100)
	# Generates a random matrix of centered variable values using the limits in the first two columns of limits
	# If the original limits on a variable are (a,b,c) with c the centering value. The limits input here
	# should be (a,b-a,c). The Constraint() function expects non-centered values.
	{
		inside<-0
		mat<-matrix(0,nTrials,N)
		q<-rep(0,N)
		while(inside<nTrials) {
			if (N>nMixture)
				q[noMixNos]<-Limits[noMixNos,1]+RandomLevels(N-nMixture,nLevels[noMixNos])*Limits[noMixNos,2]
			if (doMixture)
				q[mixNos]<-RandomMixture(nMixture)
			if (is.null(constraints) ||	constraints(q)) {
				inside<-inside+1
				if (doCenter)
					q[noMixNos]<-q[noMixNos]-Limits[noMixNos,3]
				mat[inside,]<-q
			}
		}
		for (i in noMixNos) {
			mat[,i]<-round(mat[,i],roundDigits[i])
		}
		mat
	}



	randNullify<-function(N,nCand=100) {
	# Nullify using random points
		rc<-RandomCand(N,nCand)
		j<-(sort(apply(rc*rc,1,sum),decreasing=TRUE,ind=TRUE)$ix)[1]
		design<-matrix(rc[j,],1,N)
		failure<-0
		expdes<-t(expand(design))
		while (nrow(design)<k && failure<5) {
			cand<-RandomCand(N,nCand)
			expcand<-t(expand(cand))

			aa<-qr(expdes)
			if (aa$rank==nrow(design)) {
				failure<-0
				j<-(sort(diag(t(expcand)%*%(expcand-expdes%*%qr.coef(aa,expcand))),decreasing=TRUE,ind=TRUE)$ix)[1]
				design<-rbind(design,cand[j,])
				expdes<-cbind(expdes,expcand[,j])
			}
			else
				failure<-failure+1
		}
		if (failure==5)
			stop("Singular Matrix in Nullify")
		design
	}

	expand<-function(M){
		colnames(M)<-varNames
		v<-data.frame(M)
		k<-ncol(v)
		for (i in 1:k) {
			if (isFactor[i]) {
				cl<-v[,i]
				v[,i]<-factor(cl,levels=1:nLevels[i])
			}
		}
		frm<-expand.formula(frml,varNames)
		model.matrix(frm,v)
	}



	more=TRUE
	if (criterion=="D")
		bestCrit<-(-1e8)
	else
		bestCrit<-1e8

	getCenteredDesign<-function(N,nTrials,missTrials) {
		design<-0
		if (!RandomStart) {
			# Uses nullification to find a random design
			design<-randNullify(N,nCandNull)

			q<-N
			expdes<-expand(design)
			M<-solve(t(expdes)%*%expdes)*q
			# Continues adding points til nTrials is reached
			while (q<nTrials) {
				cand<-RandomCand(N,nCandNull)
				expcand<-expand(cand)
				w<-M%*%t(expcand)
				d<-diag(expcand%*%w)
				j<-sort(d,decreasing=TRUE,ind=TRUE)$ix[1]
				design<-rbind(design,cand[j,])
				M<-((q+1)/q)*(M-(w[,j]%*%t(w[,j]))/(q+d[j]))
				q<-q+1
			}

			if (nTrials<q) {
				return(list(D=0,A=0,I=0,G=0,design=0,iter=0,error=15))
			}
		}

		iter<-0
		D<-A<-I<-Ge<-De<-0
		error<-NULL
		maxIter<-nRepeats
		aaBest<-NULL
		nT<-nTrials
		while(more){
 			cand<-RandomCand(N,nCand)
			if (!RandomStart) {
				cand<-rbind(design,cand)
				rows<-1:nTrials
			}
			else
			if (approximate) {
				rows<-1:nTrials # optFedrov() will take nTrials from the length of rows whem missTrials is TRUE
				if (missTrials)
					nT<-0
				else
					nT<-nTrials # optFederov() will round proportions to nTrials
			}
			else {
				rows=NULL
			}

			colnames(cand)<-varNames


			aa<-optFederov(frml,cand,nTrials=nT,approximate=approximate,rows=rows,
				criterion=criterion,evaluateI=evaluateI,space=space,
				DFrac=DFrac,CFrac=CFrac)
			save=FALSE
			if ((criterion=="D") && (aa$D>bestCrit)) {
				bestCrit<-aa$D
				save<-TRUE
			}else
			if ((criterion=="A") && (aa$A<bestCrit)) {
				bestCrit<-aa$A
				save<-TRUE
			}else
			if ((criterion=="I") && (aa$I<bestCrit)) {
				bestCrit<-aa$I
				save<-TRUE
			}

			if (save) {
				aaBest<-aa
			}
			else {
				more=FALSE
			}


			iter<-iter+1
			if (iter>maxIter) {
				more=FALSE
			}


		}
		aaBest
	}


	if (missing(frml))
		stop("A formula is required.")
	if (missing(data))
		stop("data is required.")

    if (criterion!="D" && criterion!="A" && criterion!="I")
      stop("Incorrect criterion")

	if (approximate && !RandomStart)
		stop("Cannot use RandomStart==FALSE with approximate==TRUE")


	N<-nrow(data)


	isFactor<-as.vector(data[,7])
	nLevels<-as.vector(data[,5])
	for (i in 1:N) {
		if (isFactor[i]) { # reset to produce integer levels, just in case
			data[i,c(2,3,4,6)]<-c(1,nLevels[i],0,1)
		}
	}
	varNames<-as.character(data[,1])
	Limits<-data.matrix(data[,c(2,3,4)])
	Limits[,2]<-Limits[,2]-Limits[,1]
	roundDigits<-as.vector(data[,6])

	if (ncol(data)==7) {
		doMixture=FALSE
		mixNos<-NULL
		noMixNos<-1:N
		nMixture<-0
	}
	else {
		mixtureVariables<-as.logical(data[,8])
		doMixture=any(mixtureVariables)
		mixNos<-(1:N)[mixtureVariables]
		nMixture<-length(mixNos)
		noMixNos<-(1:N)[!mixtureVariables]
		if (doMixture) {
			Limits[mixNos,3]<-0 # just to be safe
			mixRound<-max(roundDigits[mixNos],2)
		}
	}

	if (doMixture) { # make sure there is no constant
		frml<-deparse(frml,width.cutoff=500)
		frml<-paste(frml,"+0",sep="")
		frml<-as.formula(frml)
	}

	doCenter<-FALSE # to stop centering in this use of RandomCand
	X<-expand(RandomCand(N,1))
	k<-ncol(X)
    doCenter <- any(as.logical(match(dimnames(X)[[2]],"(Intercept)", nomatch = 0)))

	missTrials<-FALSE # used when approximate
	if (missing(nTrials) || nTrials<=0) {
		nTrials<-k+5
		missTrials<-TRUE
	}
	if (missing(nCand))
		nCand<-10*k

	if (missing(nCandNull))
		nCandNull<-nCand

	output<-getCenteredDesign(N,nTrials,missTrials)


	# remove centering
	design<-output$design
	nRows<-nrow(design)
	if (doCenter)
		output$design<-design+matrix(Limits[,3],nRows,N,byrow=TRUE)
	rownames(output$design)<-1:nRows


		# output the criterion values and the uncentered design
	output<-output[names(output)!="rows"]

	if (args)
		output<-c(output,list(seed=seed),args=actuals(formals("optMonteCarlo")))
	output
}

