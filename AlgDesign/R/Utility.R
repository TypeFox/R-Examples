# file Utility.R
# copyright (C) 2002-2004 by Robert E. Wheeler
#

"expand.formula"<-
function (frml, varNames, const = TRUE, numerics = NULL)
{
    env <- environment(frml)
    noNumerics <- missing(numerics)
    nameargs <- function(...) {
        dots <- as.list(substitute(list(...)))[-1]
        nm <- names(dots)
        fixup <- if (is.null(nm))
            seq(along = dots)
        else nm == ""
        dep <- sapply(dots[fixup], function(x) deparse(x, width.cutoff = 500)[1])
        if (is.null(nm))
            nm <- dep
        else {
            nm[fixup] <- dep
        }
        if (!noNumerics) {
            if ((nm[1] == "." && !all(numerics)) || !all(numerics[is.element(varNames,
                nm)]))
                stop("All arguments to special functions such as quad() must be numeric.")
        }
        nm
    }

	find.dot<-function (strng) {
	# Returns location of dot in strng or -1 if not found
		loc<-regexpr("([^0-9a-zA-Z_]+ *\\.)|(^ *\\.)",strng) #not variable followed by one or more spaces
		locp<-loc[1]
		if (locp!=-1) {	# dot found
			leng<-attr(loc,"match.length")
			subs<-substring(strng,locp,locp+leng-1) # substring containing dot
			locsub<-regexpr("\\.",subs) # loc of dot in substring
			return(locp+locsub[1]-1)
		}
		locp
	}

    quad <- function(...) {
		nms<-nameargs(...)
		if (-1!=find.dot(nms)) {
            nms <- varNames
            nVars <- length(nms)
        }
        else nVars <- nargs()
        strg <- paste(paste("(", paste("X", 1:nVars, sep = "",
            collapse = "+"), ")^2", sep = ""), "+", paste("I(X",
            1:nVars, "^2)", sep = "", collapse = "+"))
        for (i in 1:nVars) {
            ag <- paste("X", i, sep = "")
            strg <- gsub(ag, nms[i], strg)
        }
        strg
    }
    cubic <- function(...) {
		nms<-nameargs(...)
		if (-1!=find.dot(nms)) {
            nms <- varNames
            nVars <- length(nms)
        }
        else nVars <- nargs()
        strg <- paste(paste("(", paste("X", 1:nVars, sep = "",
            collapse = "+"), ")^3", sep = ""), "+", paste("I(X",
            1:nVars, "^2)", sep = "", collapse = "+"), "+", paste("I(X",
            1:nVars, "^3)", sep = "", collapse = "+"))
        for (i in 1:nVars) {
            ag <- paste("X", i, sep = "")
            strg <- gsub(ag, nms[i], strg)
        }
        strg
    }
    cubicS <- function(...) {
		nms<-nameargs(...)
		if (-1!=find.dot(nms)) {
            nms <- varNames
            nVars <- length(nms)
        }
        else nVars <- nargs()
        strg <- paste("(", paste("X", 1:nVars, sep = "", collapse = "+"),
            ")^3", sep = "")
        for (i in 1:(nVars - 1)) {
            var <- paste("X", i, sep = "")
            strg <- paste(strg, "+", paste(paste("I(", var, paste("*X",
                (i + 1):nVars, sep = ""), sep = ""), paste("*(",
                var, paste("-X", (i + 1):nVars, "))", sep = "")),
                collapse = "+"), sep = "", collapse = "+")
        }
        for (i in 1:nVars) {
            ag <- paste("X", i, sep = "")
            strg <- gsub(ag, nms[i], strg)
        }
        strg
    }
    findFunction <- function(name, string) {
        if (-1 == (strt <- regexpr(name, string)))
            return(c(0, 0))
        head <- substr(string, 1, strt - 1)
        tail <- substr(string, strt, nchar(string))
        if (-1 == (fin <- regexpr(")", tail)))
            return(c(0, 0))
        c(strt, strt + fin - 1)
    }
    frml <- deparse(frml, width.cutoff = 500)
    while ((0 != (pos <- findFunction("quad", frml))[1]) || (0 !=
        (pos <- findFunction("cubicS", frml))[1]) || (0 != (pos <- findFunction("cubic",
        frml))[1])) {
        prog <- substr(frml, pos[1], pos[2])
        strHead <- substr(frml, 1, pos[1] - 1)
        strTail <- substr(frml, pos[2] + 1, nchar(frml))
        prog <- eval(parse(text = prog))
        frml <- paste(strHead, prog, strTail, sep = "")
    }
    if (!const)
        frml <- paste(frml, "+0", sep = "")
    frml <- as.formula(frml)
    environment(frml) <- env
    frml
}




"actuals" <-
function(frmals){
	# returns the actual argument list of a function
		nargs<-length(frmals)
		fnames<-names(frmals)
		for (i in 1:nargs) {
			if (i==1)
				actuls=list(eval(parse(text=fnames[i]),sys.frame(-1)))
			else
				actuls<-c(actuls,list(eval(parse(text=fnames[i]),sys.frame(-1))))
		}
		names(actuls)<-fnames
		actuls
	}

"efficient.rounding" <-
function(proportions,n,random=TRUE){
# Efficient rounding, following Pulkesheim and Rieder, (1992). Efficient rounding of
#	approximate designs. Biometrika. 79, 763-770
# proportions: a vector of proportions
# n: the integer sum of the rounded proportions
# random: if TRUE will randomly select the weight to change, otherwise the first available will be changed.

	l<-length(proportions)
	n<-as.integer(n)

	if (n<l) # in case input is too small
		n<-l
	if (any(proportions<0))
		stop("Negative weights are not allowed")
	if (abs(sum(proportions)-1)>1e-12)
		stop("proportions must sum to unity")

	ni<-ceiling(proportions*(n-l/2))
	N<-sum(ni)
	while (N!=n) {
		if (N<n) {
			m<-ni/proportions
			ee<-min(m)
		}
		else {
			m<-(ni-1)/proportions
			ee<-max(m)
		}
		v<-(abs(ee-m)<1e-8)
		i<-(1:l)[v]

		if (random) {
			i<-sample(i)
		}

		if (N<n) {
			ni[i[1]]<-ni[i[1]]+1
		}
		else {
			found=FALSE
			for (j in i) {
				if (ni[j]>1) {
					ni[j]<-ni[j]-1
					found=TRUE
					break
				}
			}
			if (!found) # this will result in a zero frequency
				ni[i[1]]<-ni[i[1]]-1
		}
		N<-sum(ni)
	}
	ni
}

"eval.blockdesign" <-
function(frml,design,blocksizes,rho=1,confounding=FALSE,center=FALSE){
# evaluates a blocked design.
# a formula
# a design
# a vector of blocksizes
# rho is a vector of variance component ratios.
# confounding=TRUE, returns the confounding matrix
# center=TRUE, the numeric columns will be centered before model expansion


	rownames(design)<-1:nrow(design) # prevents warning messages in case of duplicate rownames
	frml<-expand.formula(frml,colnames(design))
	design<-data.frame(design)
	numericColumns<-sapply(design,is.numeric)
	if (center) {
		means<-apply(design[,numericColumns,drop=FALSE],2,mean)
		design[,numericColumns]<-sweep(design[,numericColumns,drop=FALSE],2,means)
	}
	des<-model.matrix(frml,design)
	desBlkCtr<-des

	if (missing(blocksizes))
		stop("blocksizes must be specified")

	if (any(rho<1e-10 || rho>1e10))
		stop("Some rho's are too large or too small")

	nB<-length(blocksizes)
	fn<-0
	S<-matrix(0,nB,ncol(des))
	for (i in 1:nB) {
		st<-fn+1
		fn<-fn+blocksizes[i]
		ad<-desBlkCtr[st:fn,]
		S[i,]<-apply(ad,2,sum)
		ad<-scale(ad,scale=FALSE)
		desBlkCtr[st:fn,]<-ad
	}
	index<-apply((abs(desBlkCtr)<1e-10),2,all) # TRUE corresponds to constant or whole block vars
	desUnblocked<-scale(des,scale=FALSE)
		# the following combines uncenterd whole terms with centered within terms
	desUnblocked<-cbind(des[,index,drop=FALSE],desUnblocked[,!index,drop=FALSE]) # rearrange columns
	des<-cbind(des[,index,drop=FALSE],desBlkCtr[,!index,drop=FALSE]) # to conform to these
	desBlkCtr<-cbind(desBlkCtr[,index,drop=FALSE],desBlkCtr[,!index,drop=FALSE])
	S<-cbind(S[,index,drop=FALSE],S[,!index,drop=FALSE])
	constantPosition<-match("(Intercept)",dimnames(des)[[2]],nomatch=0)

	k.whole<-sum(index)
	k.within<-sum(!index)
	nB<-length(blocksizes)

	wholeIndex<-c(rep(TRUE,k.whole),rep(FALSE,k.within))

	N<-nrow(desUnblocked)
	k<-ncol(desUnblocked)
	output<-NULL

	Mu<-(t(desUnblocked)%*%desUnblocked)/N
	Mw<-(t(desBlkCtr[,!wholeIndex,drop=FALSE])%*%desBlkCtr[,!wholeIndex,drop=FALSE])/N


	determinant<-det(Mu)
	if (determinant<=0)
		stop("Singular design")
	determinant<-(determinant)^(1/k)

	MI<-solve(Mu)
	Mbc<-solve((t(des)%*%des)/N)

	determinant.within<-det(Mw)
	if (determinant.within<=0)
		stop("Singular within design")
		# If the model contains nesting, the determinant.within will differ
		# from that reported by optBlock(), which uses multiplication instead
		# of nesting. This is a constant shift of no importance to the optimization
	determinant.within<-(determinant.within)^(1/k.within)

	if (confounding == TRUE) {
		Mbc<-solve((t(des)%*%des)/N)
		dg<-diag(-1/diag(Mbc))
		output <- list(confounding = round(Mbc%*%dg,4))
	}

		# calcualate the blocking efficiencies
	T<-S[,wholeIndex,drop=FALSE]
	S<-S[,!wholeIndex,drop=FALSE];
	S<-scale(S,scale=FALSE)
	Xa<-desUnblocked[,!wholeIndex,drop=TRUE]
	Xa<-scale(Xa,scale=FALSE)
	Xb<-desUnblocked[,wholeIndex,drop=TRUE]
	M<-t(Xa)%*%Xa
	Nb<-t(Xb)%*%Xb
	H<-t(Xa)%*%Xb
	D<-blocksizes

	dM<-det(M)
	trM<-sum(diag(solve(M)))

	lambda<-matrix(0,3,length(rho))
	lambda[1,]<-rho
	colnames(lambda)<-rep(" ",length(rho))
	rownames(lambda)<-c("rho","lambda.det","lambda.trace")

	lambdaD<-lambdaT<-NULL
	for (r in rho) {
		G<-diag(1/((1/r)+D))
		MM<-M-t(S)%*%G%*%S
		A1<-H-t(S)%*%G%*%T
		A2<-solve(Nb-t(T)%*%G%*%T)
		MM<-MM-A1%*%A2%*%t(A1)
		dMM<-det(MM)
		lambdaD<-c(lambdaD,round((dMM/dM)^(1/k.within),3))
		if (dMM>dM*1e-16) {
			lambdaT<-c(lambdaT,round(trM/sum(diag(solve(MM))),3))
		}
		else {
			lambdaT=c(lambdaT,0)
		}

	}
	lambda[2,]<-lambdaD
	lambda[3,]<-lambdaT

	outs<-data.frame(constant=I(rep(" ",4)),whole=I(rep(" ",4)),within=I(rep(" ",4)))
	rownames(outs)<-c("df","determinant","gmean.variance","gmean.efficiencies")

	whIndex<-wholeIndex
	if (constantPosition && k.whole>0) {
		df<-k.whole-1
		whIndex[constantPosition]<-FALSE
	}
	else
		df<-k.whole

	outs[1,]<-c(sprintf("1"),sprintf("%.0f",df),sprintf("%.0f",k.within))
	outs[2,3]<-sprintf("%f",determinant.within)


	output<-c(output,list(determinant.all.terms.within.terms.centered=determinant))
	output<-c(output,list(within.block.efficiencies=lambda))
	if (k>nB) {
		output<-c(output,list(comment="Too few blocks to recover interblock information."))
	}

	 # variances of blockcentered design
	dg<-diag(Mbc)
	if (k.whole>0) {
		if (constantPosition) {
			s<-dg[constantPosition]
			outs[3,1]<-sprintf("%f",s)
		}
		if (sum(whIndex)>0) {
			s<-dg[whIndex]
			s<-exp(mean(log(s)))
			outs[3,2]<-sprintf("%f",s)
		}
	}
	t<-dg[!wholeIndex]
	t<-exp(mean(log(t)))
	outs[3,3]<-sprintf("%f",t)


	 # efficiencies of un-blockcentered to blockcentered designs
	dg<-(diag(MI)/diag(Mbc))
	if (k.whole>0) {
		if (constantPosition) {
			s<-dg[constantPosition]
			outs[4,1]<-sprintf("%f",s)
		}
		if (sum(whIndex)>0) {
			s<-dg[whIndex]
			s<-exp(mean(log(s)))
			outs[4,2]<-sprintf("%.3f",round(s,3))
		}
	}
	t<-dg[!wholeIndex]
	t<-exp(mean(log(t)))
	outs[4,3]<-sprintf("%.3f",round(t,3))


	output<-c(output,list(block.centered.properties=outs))

	output
}

"eval.design" <-
function (frml, design, confounding = FALSE, variances = TRUE,
     center = FALSE, X=NULL)
{
	rownames(design)<-1:nrow(design) # prevents warning messages in case of duplicated rownames
	proportions<-NULL
	if (colnames(design)[1]=="Rep..") { # Expand data according to Rep.. from optFederov() output.
		ND<-nrow(design)
		reps<-design[,1]
		design<-design[,-1,drop=FALSE]
		design<-design[rep(1:ND,reps),]
		rownames(design)<-1:nrow(design)
	}
	else
	if (colnames(design)[1]=="Proportion") {
		proportions<-design[,1]
		design<-design[,-1,drop=FALSE]
	}


    frml <- expand.formula(frml, colnames(design))
    design <- data.frame(design)
    numericColumns <- sapply(design, is.numeric)
    if (center) {
        means <- apply(design[, numericColumns, drop = FALSE],
            2, mean)
        design[, numericColumns] <- sweep(design[, numericColumns,
            drop = FALSE], 2, means)
    }
    des <- model.matrix(frml, design)

	if (!missing(X)){
		X<-data.frame(X)
		if (center) {
			means <- apply(X[, numericColumns, drop = FALSE],
				2, mean)
			X[, numericColumns] <- sweep(X[, numericColumns,
				drop = FALSE], 2, means)
		}
		Xd <- model.matrix(frml, X)
	}

    output <- NULL
	if (is.vector(proportions)) {
		M<-t(des)%*%diag(proportions)%*%des
	}
	else {
		M <- (t(des) %*% des)/nrow(des)
	}
	k<-ncol(M)
	determinant<-det(M)
	if (determinant<=0)
		stop("Singular design")
    determinant<-(determinant)^(1/k)
	trac<-0
	Ival<-0
	Geff<-0
	Deff<-0
	orthog<-0
	if (determinant>0) {
		MI<-solve(M)
		if (confounding == TRUE) {
			dg<-diag(-1/diag(MI))
			output <- list(confounding = round(MI%*%dg,4))
		}
		trac<-sum(diag(MI))/k
		if (!missing(X)) {
			CC<-t(Xd)%*%Xd/nrow(Xd)
			Ival<-sum(diag(qr.solve(M,CC)))
			Z<-t(qr.solve(M,t(Xd)))
			Geff<-ncol(Xd)/max(apply(Z*Xd,1,sum))
			Deff<-exp(1-1/Geff)
		}
		cp<-match("(Intercept)",dimnames(M)[[2]],nomatch=0)
		if (cp>0) {
			orthog<-exp((log(det(M[-cp,-cp]))-sum(log(diag(M)[-cp])))/(k-1))
		}
		else {
			orthog<-exp(log(determinant)-sum(log(diag(M)))/k)
		}
	}
    output <- c(output, list(determinant = determinant,A = trac))
	if (!missing(X)) {
		output<-c(output,list(I = Ival,Geff = round(Geff,3), Deffbound = round(Deff,3)))
	}
	output<-c(output,list(diagonality = round(orthog,3)))

    if (variances) {
        if (determinant > 0) {
			t <- diag(solve(M))
			if (cp>0)
				t<-t[-cp]
            t <- exp(mean(log(t)))
            output <- c(output, list(gmean.variances = t))
        }
    }
    output
}

"gen.factorial" <-
function(levels,nVars=0,center=TRUE,factors="none",varNames=NULL){
# Generates a full factorial design
# nVars: number of factors. May be missing if length(levels)>1
# levels: a vector of levels
# center: centers non-factorial variables
# factors: variable numbers of the factors or "all" if all are factors
# varNames: the names of the variables

	if (length(levels)==1){
		levels<-rep(levels,nVars)
		if (missing(nVars))
			stop("nVars is needed when levels is a scalar")
	}
	else {
		nVars<-length(levels)
	}

	if (any(levels<1))
		stop("All levels must be greater than 1")

	N<-prod(levels)
	factorVec<-rep(0,nVars)
	if (missing(factors)) {
		factors<-NULL
	}
	else
	if (is.character(factors)) {
		if (factors=="all") {
			factors<-1:nVars
			factorVec<-rep(1,nVars)
		}
		else
			stop("factors error")
	}
	else {
		nFactors<-length(factors);
		for (i in factors)
			factorVec[i]<-1
	}

	design<-matrix(0,N,nVars)
	.Call("GetFactorial",design,as.integer(levels),as.integer(center),as.integer(factorVec),PACKAGE="AlgDesign")

	design<-data.frame(design)
	if (!missing(varNames) && length(varNames==nVars))
		colnames(design)<-varNames
	else
		colnames(design)<-paste("X",1:nVars,sep="")
	if (!missing(factors)) {
		for (i in factors) {
			design[,i]<-factor(design[,i])
		}
	}
	design
}

"model.matrix.formula" <-
function(frml,data=sys.frame(sys.parent()),...){
	if (!missing(data)) {
		if (!inherits(data,"data.frame"))
			stop("data must be a data.frame")
		if (!inherits(frml,"formula"))
			stop("frml must be a formula")
		frml<-expand.formula(frml,colnames(data))
	}
	model.matrix.default(frml,data,...)
}

"model.matrix.terms"<-
function(...){
  model.matrix.default(...)
}


"gen.mixture" <-
function(levels,vars) {

	if (is.numeric(vars)) {
		nvars<-vars
		vars<-paste("X",1:nvars,sep="")
	}
	else {
		nvars<-length(vars)
	}

	levels<-levels[1]-1

	if (levels<1)
		stop("levels error")
	if (nvars<2)
		stop("vars error")

	N<-choose(levels+nvars-1,levels);

	X<-matrix(0,N,nvars)

	.Call("GetMixture", X,as.integer(levels),PACKAGE="AlgDesign")

	X<-data.frame(X)
	colnames(X)<-vars
	X
}



