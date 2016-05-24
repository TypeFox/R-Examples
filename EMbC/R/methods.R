# The EMbC Package for R
#
# Copyright 2013, 2014, 2015 Joan Garriga <jgarriga@ceab.csic.es>, Aitana Oltra
# <aoltra@ceab.csic.es>, John R.B. Palmer <johnrbpalmer@gmail.com>, Frederic
# Bartumeus <fbartu@ceab.csic.es>
#
# EMbC is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or  (at your option) any later
# version.
#
# EMbC is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see http://www.gnu.org/licenses.

# Methods with no generic
# -----------------------

#' @title Select a single path from a \code{binClstStck} instance.
#'
#' @description \code{slct} selects a single path from a \code{binClstStck}
#'   instance.
#'
#' @param stck A \code{binClstStck} instance.
#'
#' @param pathNmbr The number of the single path to be selected.
#'
#' @return Returns the single \link{binClstPath_instance} selected.
#'
#' @export
#' @rdname slct
#'
#' @examples
#' \dontrun{
#' # -- select path number 3 in mybcpstack --
#' bcp3 <- slct(mybcpstack,3)
#' }

slct <- function(stck,pathNmbr){
	if (pathNmbr<=length(stck@bCS))
		return(stck@bCS[[pathNmbr]])
	else
		cat('Error: selected path does not exist !!\n')
	}

#' @title Sets binClst color palette .
#'
#' @description \code{setc} sets the color palette to a color family from the
#'   RColorbrewer package
#'
#' @param bC A \link{binClst_instance}.
#'
#' @param fam The name of a color family from the Rcolorbrewer R-package,
#'   (default color palette is 'RdYlBu' which is colorblind safe and print
#'   friendly up to 6 colors).
#'
#' @export
#' @rdname setc
#'
#' @examples
#' # -- change the color palette of mybc to "PuOr" --
#' \dontrun{
#' setc(mybc,'PuOr')
#' }

setc <- function(bC,fam='RdYlBu'){
	objName <- deparse(substitute(bC))
	bC@C <- getColors(bC@k,fam)
	assign(objName,bC,envir=parent.frame())}

# binClst Methods
# ---------------

#' @title Likelihood profile plots
#'
#' @description \code{lkhp} likelihood optimization plot.
#'
#' @param obj A \code{BinClst_instance} or a list of them.
#'
#' @param offSet A numeric value indicating an offset to avoid the initial
#'   iterations. This is usefull to see the likelihood evolution in the last
#'   iterations where the changes in likelihood are of different order of
#'   magnitud than those at the starting iterations.
#'
#' @export
#' @rdname lkhp
#'
#' @examples
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth)
#' # -- inspect the likelihood evolution --
#' lkhp(mybcp)
#' # -- avoid the initial values --
#' lkhp(mybcp,10)

setGeneric("lkhp",function(obj,offSet=1){standardGeneric("lkhp")})

#' @rdname lkhp
setMethod("lkhp",signature(obj="binClst"),function(obj,offSet){
	bC <- obj
	parDef <- parSet()
	plot(seq(offSet,length(bC@L)),bC@L[offSet:length(bC@L)],col=4,type='l',xlab='itr',ylab="logLikelihood")
	par(parDef)
	})

#' @rdname lkhp
setMethod("lkhp",signature(obj="list"),function(obj,offSet){
	parDef <- parSet()
	lkhLst <- lapply(obj,function(bC) bC@L[offSet:length(bC@L)])
	maxLkh <- max(sapply(lkhLst,max))
	minLkh <- min(sapply(lkhLst,min))
	maxLen <- max(sapply(lkhLst, length))
 	lkhLst <- lapply(lkhLst,function(l) l[1:maxLen])
	plot(seq(maxLen),lkhLst[[1]],type='l',xlab='iters',ylab='logLikelihood',ylim=c(minLkh,maxLkh))
	lapply(seq(2:length(lkhLst)), function(i) lines(seq(maxLen),lkhLst[[i]],col=i))
	grid()
 	cat(as.character(lapply(obj,function(bc) length(which(diff(bc@L)<0)))),'\n')
 	cat(as.character(lapply(obj,function(bc) round(sum(bc@L[which(diff(bc@L)<0)]),4))),'\n')
 	par(parDef)
	})

#' @title Clustering statistics.
#'
#' @description \code{stts} clustering statistics information.
#'
#' @param obj Either a \link{binClst_instance} or a \code{binClstStck} instance.
#'   In the latter case statistics are given at stack level.
#'
#' @param dec The number of decimals for mean/stdv formatting.
#' @param width The number of digits for mean/stdv formatting.
#'
#' @details
#'
#' This method prints a line for each cluster with the following information:
#' the cluster number, the cluster binary label, the cluster mean and variance
#' of each input feature (two columns for each variable), and the size of the
#' cluster in number and proportion of points (the posterior marginal
#' distribution).
#'
#' @export
#' @rdname stts
#'
#' @examples
#' # -- apply EMbC to the example path with solar covariate 'height'--
#' mybcp <- stbc(expth,scv='height',info=-1)
#' # -- show clustering statistics --
#' stts(mybcp,width=5,dec=1)
#' \dontrun{
#' # -- show clustering statistics of mybcpstack at stack level --
#' stts(mybcpstack)
#' # -- show individual statistics for path number 3 in mybcpstack --
#' stts(slct(mybcpstack,3))
#' }

setGeneric("stts",function(obj,dec=2,width=8){standardGeneric("stts")})

#' @rdname stts
setMethod("stts",signature(obj="binClst"),function(obj,dec,width){
	bC <- obj
	width <- max(width,5)
	cat(paste(rep(' ',(bC@m+3)),collapse=''),sapply(seq(bC@m),function(m) paste(paste(rep(' ',(width-5)),collapse=''),'X',m,'.mn',paste(rep(' ',(width-4)),collapse=''),'X',m,'.sd',sep='')),'    kn  kn(%)\n')
	for (k in seq(bC@k)){
		kn <- length(which(bC@A==k))
		if (kn>0){
			bk  <- paste(as.integer(rev(intToBits(k-1)[1:bC@m])),collapse="")
			kLbl <- paste(formatC(k,width=2),gsub('1','H',gsub('0','L',bk)),sep=" ")
			kMrg <- formatC(round(kn/bC@n*100,2),format='f',digits=2,width=6)
			cat(paste(kLbl,frmTht(bC@P[[k]],dec,width),formatC(kn,width=6),kMrg,'\n'))
			}
		}
	})

#' @rdname stts
setMethod("stts",signature(obj="binClstStck"),function(obj,dec,width){
	stts(obj@bC,dec,width)
	})

#' @title Clustering 2D-scatterplot
#'
#' @description \code{sctr} generates a scatterplot from a
#'   \link{binClst_instance}, showing clusters in different colors.
#'
#' @param obj A \link{binClst_instance}.
#'
#' @param ref
#'
#'   A numeric vector with expert/reference labelling for visual validation of
#'   the clustering.
#'
#'   A second \link{binClst_instance} to be compared with the former.
#'
#' @param showVars When the number of variables is greater than two, a length 3
#'   numeric vector indicating one splitting variable and two variables to be
#'   scattered (given in that order).
#'
#' @param showClst When the number of variables is greater than two, a numeric
#'   vector (of variable length) indicating a subset of the clusters that will
#'   be shown in the scatter plot. This is usefull in case of overlapping
#'   clusters.
#'
#' @param bg A valid colour especification to be used as background colour for multivariate scatterplots. By default a light-grey colour is used to enhance data points visibility.
#'
#' @param ... Parameters \code{ref}, \code{showVars} and \code{showClst} are
#'   optional.
#'
#' @export
#' @rdname sctr
#'
#' @examples
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth,info=-1)
#' # -- show the scatterplot compared with expert labelling--
#' sctr(mybcp,expth$lbl)

setGeneric("sctr",function(obj,...){standardGeneric("sctr")})

#' @rdname sctr
setMethod("sctr",signature(obj="binClst"),function(obj,ref=NULL,showVars=NULL,showClst=NULL, bg=NULL, ...){
	if (!is.null(ref)){
		if (class(ref) == 'numeric'){
			bCX <- obj
			bCX@A <- ref
			}
		else{
			bCX <- ref
			}
		if (obj@m==2){
			parDef <- parSet(mtx=matrix(c(1,2),2,1,byrow=TRUE), widths=c(1,1,0.5))
			sctr2D(obj)
			sctr2D(bCX)
			}
		else{
			if (is.null(bg)) bg <- brewer.pal(8,'Greys')[1]
			parDef <- parSet(mtx=matrix(c(1,2,3),1,3,byrow=TRUE), widths=c(1,1,0.5), bg=bg, oma=c(4,1,4,1))
			sctr3D(bCX,showVars,showClst)
			}
		}
	else{
		if (obj@m==2){
			parDef <- parSet()
			sctr2D(obj)
			}
		else{
			if (is.null(bg)) bg <- brewer.pal(8,'Greys')[1]
			parDef <- parSet(mtx=matrix(c(1,2,3),1,3,byrow=TRUE),widths=c(1,1,0.5), bg=bg, oma=c(4,1,4,1))
			sctr3D(obj,showVars,showClst)
			}
		}
	par(parDef)
	})

#' @rdname sctr
setMethod("sctr",signature(obj="binClstStck"),function(obj,ref=NULL,showVars=NULL,showClst=NULL,...){
	sctr(obj@bC,showVars=showVars,showClst=showClst)
	})

#' @title Dynamic 3D-scatterplot
#'
#' @description \code{sct3} generates a dynamic 3D-scatterplot of a multivariate
#'   \link{binClst_instance}, showing clusters in different colors. The scatter
#'   plot can be zoomed/rotated with the mouse.
#'
#' @inheritParams sctr
#'
#' @details This function needs the package "rgl" to be installed.
#'
#' @export
#' @rdname sct3
#'
#' @examples
#' \dontrun{
#' # -- apply EMbC to the example path with scv='height' --
#' mybcp <- stbc(expth,scv='height')
#' # -- show a dynamic 3D-scatterplot --
#' sct3(mybcp)
#' # -- show only a subset of clusters --
#' sct3(mybcp,showClst=c(2,4,6))
#' }

setGeneric("sct3",function(obj,...){standardGeneric("sct3")})

#' @rdname sct3
setMethod("sct3",signature(obj="binClst"),function(obj,showVars=NULL,showClst=NULL,...){
  # check rgl is available
  if (!requireNamespace('rgl',quietly=TRUE)){
    return(cat(' ... this function needs package "rgl" to be installed!\n'))
  }
  # set variable and cluster subsets
  lims <- apply(obj@X,2,range)
	if (!is.null(showVars) && length(showVars)>=3) m <- showVars[1:3]
	else m <- c(1,2,3)
	X <- obj@X
	A <- obj@A
	if (!is.null(showClst)) {
		X <- X[which(A %in% showClst),]
		A <- A[which(A %in% showClst)]}
	# 3d scatterplot
	parDef <- parSet(bg=brewer.pal(4,'Greys')[4])
	rgl::bg3d(brewer.pal(8,'Greys')[4])
#~ 	rgl::par3d('font'=3,'cex'=0.7) # works only on windows !!!
	rgl::par3d('family'='bitmap','font'=1,'cex'=1)
	rgl::plot3d(X[,m[1]],X[,m[2]],X[,m[3]],xlab=paste('X',m[1],sep=''),ylab=paste('X',m[2],sep=''),zlab=paste('X',m[3],sep=''),xlim=lims[,m[1]],ylim=lims[,m[2]],zlim=lims[,m[3]],type="s",size=1,col=obj@C[A],top=TRUE,legend=TRUE)
	rgl::grid3d(c("x", "y+", "z"))
 	par(parDef)
	})

#' @title Labeling profile plot
#'
#' @description \code{lblp} plots the labeling profile of a
#'   \link{binClst_instance}.
#'
#' @param obj A \link{binClst_instance}.
#'
#' @param ref
#'
#'   A numeric vector with an expert's labelling profile.
#'
#'   A second \link{binClst_instance} to be compared with the first.
#'
#' @param lims A numeric vector with lower and upper bounds to limit the plot.
#'
#' @param ... Parameters \code{ref} and \code{lims} are optional.
#'
#' @export
#' @rdname lblp
#'
#' @examples
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth)
#' # -- plot the labeling profile comparing with expert labelling --
#' lblp(mybcp,expth$lbl)
#' # -- compare original and smoothed labeling profiles --
#' lblp(mybcp,smth(mybcp))

setGeneric("lblp",function(obj,ref,...){standardGeneric("lblp")})

#' @rdname lblp
setMethod("lblp",signature(obj="binClst",ref="missing"),function(obj,ref,lims=NULL,...){
	parDef <- parSet(mtx=matrix(c(1,2),2,1,byrow=TRUE),mar=c(2.5,2.5,0.5,0.5))
	if (is.numeric(lims) && (length(lims)==2)) {a <- max(1,lims[1]); b <- min(lims[2],(obj@n-1))}
	else {a <- 1; b <- (obj@n-1)}
	plot(range(a:b),range(0:1),type='n',xlab='loc.',ylab='X',yaxt='n')
 	lapply(1:obj@m, function(m) lines(seq(a,b),maxminScale(obj@X[a:b,m]),col=obj@C[2*m]))
	legend("topleft",legend=colnames(obj@X),col=obj@C[seq(obj@m)*2],cex=0.8,lwd=3,text.font=1,bty='n')
	grid()
	plot(seq(a,b),obj@A[a:b],col=4,type='l',xlab="loc.",ylab="",yaxt="n")
	grid()
	axis(2,at=seq(obj@k),labels=c(getkLbls(obj)),las=2,cex.axis=0.7)
 	par(parDef)
	})

#' @rdname lblp
setMethod("lblp",signature(obj="binClstStck",ref="missing"),function(obj,ref,lims=NULL,...){
	lblp(obj@bC,lims)
	})

#' @rdname lblp
setMethod("lblp",signature(obj="binClst",ref="numeric"),function(obj,ref,lims=NULL,...){
	parDef <- parSet()
	difA <- rep(0,max(obj@n,length(ref)))
	difA[which(obj@A[1:min(obj@n,length(ref))]!=ref[1:min(obj@n,length(ref))])] <- 1
	if (is.numeric(lims) && (length(lims)==2)) {a <- max(1,lims[1]); b <- min(lims[2],(length(difA)-1))}
	else {a <- 1; b <- length(difA) -1}
	plot(range(a:b),range(0:(2*obj@k+1)),type='n',xlab='loc.',ylab='',yaxt='n')
	lines(seq(a,b),(obj@A[a:b]+(obj@k+1)),col=4)
	lines(seq(a,b),(ref[a:b]+1),col=3)
	lines(seq(a,b),difA[a:b],col=2)
	grid()
	axis(2,at=1.0+seq(2*obj@k),labels=rep(getkLbls(obj),2),las=2,cex.axis=0.7)
	par(parDef)
	})

#' @rdname lblp
setMethod("lblp",signature(obj="binClst",ref="binClst"),function(obj,ref,lims=NULL,...){
	lblp(obj,ref@A,lims)
	})

#' @title Check labeling profile
#'
#' @description Plots the labeling profile of a
#'   \link{binClst_instance} against a control variable (e.g. environmental
#'   information) depicted as background coloured bars.
#'
#' @param obj A \link{binClst_instance}.
#'
#' @param ctrlLbls A numeric vector with the control labels or a string
#'   specifying one of 'height', 'azimuth' or 'both' solar covarites. By
#'   default, for a \link{binClstPath_instance} it is set to the solar height
#'   covariate, regardless it has been used or not for the clustering.
#'
#' @param ctrlClrs A vector of colors to depict the control labelling. At least
#'   one colour should be specified for each different control label. By default
#'   white/grey colours are used for the default control labels.
#'
#' @param ctrlLgnd A vector of strings identifying the labels for the legend of
#'   the plot. They are automatically generated for the solar covariates.
#'
#' @param lims A numeric vector with lower and upper bounds to limit the plot.
#'
#' @param ... Parameters \code{ctrLbls}, \code{ctrlClrs}, \code{ctrlLgnd} and
#'   \code{lims} are optional.
#'
#' @export
#' @docType methods
#' @rdname chkp
#'
#' @examples
#' # -- apply EMbC to \code{expth} --
#' mybcp <- stbc(expth)
#' # -- plot the labeling profile against 'both' solar covariates --
#' chkp(mybcp,ctrlLbls='both',ctrlClrs=RColorBrewer::brewer.pal(8,'Oranges')[1:4])

setGeneric("chkp",function(obj,...){standardGeneric("chkp")})

#' @rdname chkp
setMethod("chkp",signature(obj="binClst"),function(obj,ctrlLbls=NULL,ctrlClrs=NULL,ctrlLgnd=NULL,lims=NULL){

	if (is.null(ctrlLbls)){
		if (class(obj)!='binClst') ctrlLbls <- 'height'
		else {
			cat('ctrlLbls should be provided \n')
			return
			}
		}
	if (class(ctrlLbls)=='character' && ctrlLbls %in% c('height','azimuth','both')){
		solP <- getSolarPos(obj@pth,'both')
		if (ctrlLbls=='height'){
			cLbls <- as.integer(solP[,2]<0) +1
			if (is.null(ctrlLgnd)) cLgnd <- c('daytime','nighttime')
			if (is.null(ctrlClrs)) cClrs <- brewer.pal(8,'Greys')[c(2,4)]
			}
		if (ctrlLbls=='azimuth'){
			cLbls <- as.integer(solP[,1]>180) +1
			if (is.null(ctrlLgnd)) cLgnd <- c('forenoon','afternoon')
			if (is.null(ctrlClrs)) cClrs <- brewer.pal(8,'Greys')[c(2,4)]
			}
		if (ctrlLbls=='both'){
			cMtrx <- matrix(c(0,0,1,0,1,1,0,1),c(4,2),byrow=TRUE)
			cLbls <- sapply(seq(obj@n),function(i) which(apply(cMtrx,1,function(l) (l[1]==as.integer(solP[i,1]>180) && l[2]==as.integer(solP[i,2]<0)))))
			if (is.null(ctrlLgnd)) cLgnd <- c('morning','afternoon','evening','night')
			if (is.null(ctrlClrs)) cClrs <- brewer.pal(8,'OrRd')[2:5]
			}
		}
	else {
		cLbls <- ctrlLbls
		if (is.null(ctrlLgnd)) cLgnd <- names(table(cLbls))
		if (is.null(ctrlClrs)) cClrs <- brewer.pal(8,'Greys')[2:min(length(names(table(cLbls))),8)]
		}
	if (!is.null(ctrlClrs)) cClrs <- ctrlClrs
	if (!is.null(ctrlLgnd)) cLgnd <- ctrlLgnd
	if (is.numeric(lims) && (length(lims)==2)) {a <- max(1,lims[1]); b <- min(lims[2],(obj@n-1))}
	else {a <- 1; b <- (obj@n-1)}

	parDef <- parSet()
	plot(range(a:b),range(0:obj@k+0.5),type='n',xlab='loc.',ylab='',yaxt='n')
	rect(seq(a,b-1),1,seq(a+1,b),obj@k,col=cClrs[cLbls[a:b]],border=NA)
	lines(seq(a,b),obj@A[a:b],col=4,type='l')
	grid()
	axis(2,at=seq(obj@k),labels=c(getkLbls(obj)),las=2,cex.axis=0.7)
	legend("bottom",legend=cLgnd,col=cClrs,cex=0.8,lwd=3,text.font=1,bty='n',horiz=TRUE)
	par(parDef)

	})

#' @title Confusion matrix
#'
#' @description \code{cnfm} computes the confusion matrix of the clustering with
#'   respect to an expert/reference labelling of the data. Also, it can be used
#'   to compare the labellings of two different clusterings of the same
#'   trajectory, (see details).
#'
#' @param obj A \link{binClst_instance} or \code{bnClstStck} instance.
#'
#' @param ref
#'
#'   A numeric vector with an expert/reference labelling of the data.
#'
#'   A second \link{binClst_instance} (see details).
#'
#' @param ret A boolean value (defaults to FALSE). If ret=TRUE the confusion
#'   matrix is returned as a matrix object.
#'
#' @param ... Parameters \code{ref} and \code{ret} are optional.
#'
#' @details
#'
#' The confusion matrix yields marginal counts and Recall for each row, and
#' marginal counts, Precission and class F-measure for each column. The 3x2
#' subset of cells at the bottom right show (in this order): the overall
#' Accuracy, the average Recall, the average Precision, NaN, NaN, and the
#' overall Macro-F-Measure. The number of classes (expert/reference labelling)
#' should match or, at least not be greater than the number of clusters. The
#' overall value of the Macro-F-Measure is an average of the class F-measure
#' values, hence it is underestimated if the number of classes is lower than the
#' number of clusters.
#'
#' If \code{obj} is a \link{binClstPath_instance} and there is a column "lbl" in
#' the obj@@pth slot with an expert labelling, this labelling will be used by
#' default.
#'
#' If \code{obj} is a \code{binClstStck} instance and, for all paths in the
#' stack, there is a column "lbl" in the obj@@pth slot of each, this labelling
#' will be used to compute the confusion matrix for the whole stack.
#'
#' If \code{obj} and \code{ref} are both a \link{binClst_instance} (e.g.
#' smoothed versus non-smoothed), the confusion matrix compares both labellings.
#'
#' @return If ret=TRUE returns a matrix with the confusion matrix values.
#'
#' @export
#' @rdname cnfm
#'
#' @examples
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth,info=-1)
#' # -- compute the confusion matrix --
#' cnfm(mybcp,expth$lbl)
#' # -- as we have expth$lbl the following also works --
#' cnfm(mybcp,mybcp@@pth$lbl)
#' # -- or simply --
#' cnfm(mybcp)
#' # -- numerical differences with respect to the smoothed clustering --
#' cnfm(mybcp,smth(mybcp))

setGeneric("cnfm",function(obj,ref,...){standardGeneric("cnfm")})

#' @rdname cnfm
setMethod("cnfm",signature(obj="binClst",ref="numeric"),function(obj,ref,ret=FALSE,...){
	bC <- obj
	cMtx <- matrix(rep(0,(bC@k+3)*(bC@k+2)),c((bC@k+3),(bC@k+2)))
	for (i in seq(bC@k)){
		for (j in seq(bC@k)){
			cMtx[i,j] <- length(which(ref[which(bC@A==i)]==j))
			}
		}
	# marginal sums
	cMtx[,(bC@k+1)] <- apply(cMtx,1,sum)
	cMtx[(bC@k+1),] <- apply(cMtx,2,sum)
	# precision
	for (i in seq(bC@k)) cMtx[i,(bC@k+2)] <- ifelse(cMtx[i,(bC@k+1)]>0,round(cMtx[i,i]/cMtx[i,(bC@k+1)],2),0)
	# recall
	for (j in seq(bC@k)) cMtx[(bC@k+2),j] <- ifelse(cMtx[(bC@k+1),j]>0,round(cMtx[j,j]/cMtx[(bC@k+1),j],2),0)
	# marginal distribution
	cMtx[1:bC@k,(bC@k+1)] <- round(cMtx[1:bC@k,(bC@k+1)]/cMtx[(bC@k+1),(bC@k+1)],2)
	cMtx[(bC@k+1),1:bC@k] <- round(cMtx[(bC@k+1),1:bC@k]/cMtx[(bC@k+1),(bC@k+1)],2)
	# accuracy
	cMtx[(bC@k+1),(bC@k+1)] <- round(sum(unlist(lapply(1:bC@k,function(i) cMtx[i,i])))/(bC@n-1),2)
	# average precision
#	cMtx[(bC@k+1),(bC@k+2)] <- round(mean(cMtx[which(cMtx[1:bC@k,(bC@k+2)]>0),(bC@k+2)]),2)
	cMtx[(bC@k+1),(bC@k+2)] <- round(mean(cMtx[1:bC@k,(bC@k+2)]),2)
	# average recall
#	cMtx[(bC@k+2),(bC@k+1)] <- round(mean(cMtx[(bC@k+2),which(cMtx[(bC@k+2),1:bC@k]>0)]),2)
	cMtx[(bC@k+2),(bC@k+1)] <- round(mean(cMtx[(bC@k+2),1:bC@k]),2)
	# class Fm
	cMtx[(bC@k+3),1:bC@k] <- unlist(lapply(1:bC@k,function(i) ifelse((cMtx[i,bC@k+2]+cMtx[bC@k+2,i])>0,round(2*cMtx[i,bC@k+2]*cMtx[bC@k+2,i]/(cMtx[i,bC@k+2]+cMtx[bC@k+2,i]),2),0)))
	# macro_Fm (average Fms)
#	cMtx[(bC@k+3),(bC@k+2)] <- round(mean(cMtx[(bC@k+3),which(cMtx[(bC@k+3),1:bC@k]>0)]),2)
	cMtx[(bC@k+3),(bC@k+2)] <- round(mean(cMtx[(bC@k+3),1:bC@k]),2)
	# non-significant matrix elements
	cMtx[(bC@k+2),(bC@k+2)] <- NaN
	cMtx[(bC@k+3),(bC@k+1)] <- NaN
	rowlbls <- getkLbls(bC,kNmbrs=TRUE)
 	rownames(cMtx) <- c(rowlbls,formatC(c('mrg ','Rcl ','Fms '),width=nchar(rowlbls[[1]]),flag='.'))
	colnames(cMtx) <- c(paste('cls.',formatC(1:bC@k,width=2,format='d',flag='0'),sep=''),'  mrg.','  Prc.')
	print(cMtx[1:bC@k,])
 	colnames(cMtx) <- rep('------',(bC@k+2))
	print(cMtx[(bC@k+1):(bC@k+3),])
	colnames(cMtx) <- c(paste('cls.',formatC(1:bC@k,width=2,format='d',flag='0'),sep=''),'  mrg.','  Prc.')
	if (ret) return(cMtx)
	})

#' @rdname cnfm
setMethod("cnfm",signature(obj="binClstPath",ref="missing"),function(obj,ref,ret=FALSE,...){
	if ('lbl' %in% names(obj@pth))
		return(cnfm(obj,obj@pth$lbl,ret))
	else
		cat(paste('Error: no reference labels for ',deparse(substitute(obj)),'\n',sep=''))
	})

#' @rdname cnfm
setMethod("cnfm",signature(obj="binClstStck",ref="missing"),function(obj,ref,ret=FALSE,...){
	if (all(unlist(lapply(obj@bCS,function(bci) (class(bci)=='binClstPath' && 'lbl' %in% names(bci@pth))))))
		return(cnfm(obj@bC,as.numeric(unlist(lapply(obj@bCS,function(bci) bci@pth['lbl']))),ret))
	else
		cat(paste('Error: no reference labels for ',deparse(substitute(obj)),'\n',sep=''))
	})

#' @rdname cnfm
setMethod("cnfm",signature(obj="binClst",ref="binClst"),function(obj,ref,ret=FALSE,...){
	if (length(obj@A)==length(ref@A))
		return(cnfm(obj,ref@A,ret))
	else
		cat('Error: instances length mismatch\n')
	})

#' @title Posterior smoothing of single local labels.
#'
#' @description \code{smth} Performs a posterior smoothing of single local
#'   labels (locations that differ from their neighbouring locations while the
#'   later have equal labels).
#'
#' @param obj A \code{binClst_instance}.
#'
#' @param dlta A numeric value in the range (0,1) (default is 1) indicating the
#'   user's will to accept a change of label. The change of label is done
#'   whenever the decrease in likelihood is not greater then \code{dlta}.
#'
#' @return A smoothed copy of the input instance.
#'
#' @export
#' @rdname smth
#'
#' @examples
#' # -- cluster the example path with a prior smooth of 1 hour --
#' mysmoothbcp <- stbc(expth,smth=1,info=-1)
#' # -- apply a posterior smoothing --
#' mysmoothbcpsmoothed <- smth(mysmoothbcp,dlta=0.5)

setGeneric("smth",function(obj,dlta=1){standardGeneric("smth")})

#' @rdname smth
setMethod("smth",signature(obj="binClst"),function(obj,dlta=1){
		if (dlta>1) return(cat('"dlta" must be in the range (0,1)\n'))
		return(postSmth(obj,dlta))})

#' @title Manual relabelling of clusters.
#'
#' @description \code{rlbl} Manual relabelling of clusters (to merge clusters or
#'   relabel merged clusters).
#'
#' @param obj A \link{binClst_instance}.
#'
#' @param old The number of the cluster to be relabelled.
#'
#' @param new The new number of the cluster.
#'
#' @param reset A boolean value (defaults to FALSE). If reset=TRUE the labelling
#'   is reset to the original state.
#'
#' @details
#'
#' Whenever two adjacent clusters are merged, the label identifying the
#' splitting variable between them both is meaningless, and the algorithm ends
#' up assigning either a L or H only depending on how it evolved until reaching
#' the merging point. Thus it can happen that the final labelling of the
#' resulting cluster is not the most intuitive one. With this method the labels
#' can be changed as desired. It can also be used to manually force the merging
#' of two clusters.
#'
#' This method does not return a rellabeled copy of the input \code{obj},
#' instead the \link{binClst_instance} itself is relabelled. However, this is
#' intended only for output and visualization pourposes (sctr(), lblp(),
#' cnfm(), view()) as the \link{binClst_instance} parameters (GMM parameters and
#' binary delimiters) are not recomputed. Thus the input instance can always be
#' reset to its original state.
#'
#' @return This method does not return a rellabeled copy of the input
#'   \code{obj}, instead the \link{binClst_instance} itself is relabelled. It is
#'   intended only for visualization pourposes, as it does not recompute the GMM
#'   parameters nor the binary delimiters of the \link{binClst_instance}.
#'
#' @export
#' @rdname rlbl
#'
#' @examples
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth,info=-1)
#' # -- manually merge clusters 1 and 2 --
#' rlbl(mybcp,1,2)
#' # -- reset to the original state --
#' rlbl(mybcp,reset=TRUE)

setGeneric("rlbl",function(obj,old=0,new=0,reset=FALSE){standardGeneric("rlbl")})

#' @rdname rlbl
setMethod("rlbl",signature(obj="binClst"),function(obj,old=0,new=0,reset=FALSE){
	bC <- obj
	if (!reset && !(all(c(old,new) %in% seq(bC@k))))
		cat('nothing relabelled!!\n')
	else{
		if (reset) bC@A <- getClusters(bC)
		else bC@A[which(bC@A==old)] <- new
		bC@P <- getTheta(bC,bC@stdv**2)
		assign(deparse(substitute(obj)),bC,envir=parent.frame())
		}
	})


# binClstPath Methods
# -------------------

#' @title Variables' profile plots
#'
#' @description \code{varp} easy plot of input, output and intermediate
#'   variables of a \link{binClstPath_instance}.
#'
#' @param obj Either a matrix or a  \link{binClstPath_instance}.
#'
#' @param lims A numeric vector with lower and upper bounds to limit the plot.
#'
#' @param ... Parameter \code{lims} is optional.
#'
#' @details
#'
#' If \code{obj} is a matrix, axes labels are automatically generated from the
#' \code{colnames()} of the matrix, hence they can be changed as desired.
#'
#' If \code{obj} is a \link{binClstPath_instance} it plots the values of the
#' intermediate computations saved in slots mybcp@@spn (span times), mybcp@@dst
#' (distances) and mybcp@@hdg (local heading directions).
#'
#' @export
#' @rdname varp
#'
#' @examples
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth,info=-1)
#' # -- plot clustering data points --
#' varp(mybcp@@X)
#' # -- plot data points' certainties --
#' varp(mybcp@@U)
#' # -- plot intermediate computations (span-times, distances and headings) in one figure --
#' varp(mybcp)
#' \dontrun{
#' # -- plot only span-times between locations a and b --
#' plot(seq(a,b),mybcp@@spn[a:b],col=4,type='l',xlab='loc',ylab='spanTime (s)')
#' }

setGeneric("varp",function(obj,...){standardGeneric("varp")})

#' @rdname varp
setMethod("varp",signature(obj="binClstPath"),function(obj,lims=NULL,...){
	parDef <- parSet(mtx=matrix(c(1,2,3),3,1,byrow=TRUE))
	if (is.numeric(lims) && (length(lims)==2)) {a <- max(1,lims[1]); b <- min(lims[2],obj@n)}
	else {a <- 1; b <- obj@n}
	plot(seq(a,b),obj@spn[a:b],col=4,type='l',xlab='loc',ylab='span-time (s)')
	plot(seq(a,b),obj@dst[a:b],col=4,type='l',xlab='loc',ylab='distance (m)')
	plot(seq(a,b),obj@hdg[a:b],col=4,type='l',xlab='loc',ylab='heading (rad)')
	par(parDef)
	})

#' @rdname varp
setMethod("varp",signature(obj="matrix"),function(obj,lims=NULL,...){
	parDef <- parSet(mtx=matrix(seq(dim(obj)[2]),dim(obj)[2],1,byrow=TRUE))
	if (is.numeric(lims) && (length(lims)==2)) {a <- max(1,lims[1]); b <- min(lims[2],nrow(obj))}
	else {a <- 1; b <- nrow(obj)}
	plts <- lapply(seq(ncol(obj)),function(i){
		plot(seq(a,b),obj[a:b,i],col=4,type='l',xlab='loc',ylab=colnames(obj)[i])})
	par(parDef)
	})


# Path visualization Methods
# --------------------------

#' @title Path fast view
#'
#' @description \code{view} provides a fast plot of a segmented trajectory or
#'   specific chunks of it.
#'
#' @param obj A \link{binClstPath_instance} or a data.frame with the format
#'   described for slot \code{binClstPath@@pth}.
#'
#' @param lbl A numeric vector with location labels. If \code{obj} is a
#'   \link{binClstPath_instance} the clustering labels are used by default.
#'
#' @param lims A numeric vector with lower and upper limit locations to show
#'   only a chunk of the trajectory.
#' @param bg A valid colour especification to be used as background colour. By default a light-grey colour is used to enhance data points visibility.
#'
#' @param ... Parameters \code{lbl} and \code{lims} are optional.
#'
#' @export
#' @rdname view
#'
#' @examples
#' # -- Fast view of the binClstPath instance included in the package --
#' view(expth)
#' # -- the same with reference labels --
#' view(expth,lbl=TRUE)

setGeneric("view",function(obj,...){standardGeneric("view")})

#' @rdname view
setMethod("view",signature(obj="binClstPath"),function(obj, lbl=NULL, lims=NULL, bg=NULL, ...){
	if (is.null(bg)) bg <- brewer.pal(8,'Greys')[4]
	parDef <- parSet(mtx=matrix(c(1,2),2,1), heights=c(0.2,1), bg=bg)
	if (is.numeric(lbl) && (length(lbl)>0)) pal <- getColors(length(table(lbl)))
	else {lbl <- obj@A; pal <- obj@C}
	if (is.numeric(lims) && (length(lims)==2)) {a <- max(1,lims[1]); b <- min(lims[2],obj@n)}
	else {a <- 1; b <- obj@n}
	par(mar=c(1,4,1,1))
	plot(1, 1, xlab='', ylab='', xlim=c(a,b), ylim=c(0,1), xaxt="n", yaxt="n", bty="n", type="n")
	rect(seq(a,b-1), 0, seq(a+1,b), 1, col=pal[lbl[a:b]], border=NA)
	par(mar=c(4,4,1,1))
	plims <- getPropLims(obj@pth, a, b)
	plot(1, 1, xlab='lon', ylab='lat', xlim=plims$x, ylim=plims$y, type="n", asp=1)
	points(obj@pth$lon[a:b], obj@pth$lat[a:b], col=pal[lbl[a:b]], pch=20)
	lines(obj@pth$lon[a:b], obj@pth$lat[a:b], col=pal[length(pal)])
#  	world(xlim=plims$x, ylim=plims$y)
 	par(parDef)
	})

#' @rdname view
setMethod("view",signature(obj="data.frame"),function(obj, lbl=NULL, lims=NULL, bg=NULL, ...){
	if (is.null(bg)) bg <- brewer.pal(8,'Greys')[4]
	parDef <- parSet(mtx=matrix(c(1,2),2,1), heights=c(0.2,1), bg=bg)
	colnames(obj)[c(1,2,3)] <- c('dTm','lon','lat')
	if (!(is.numeric(lbl) && (length(lbl)>0))){
		if ('lbl' %in% names(obj)) lbl <- obj$lbl
		else lbl <- rep(1,nrow(obj))
		}
	pal <- getColors(length(table(lbl)))
	if (is.numeric(lims) && (length(lims)==2)) {a <- max(1,lims[1]); b <- min(lims[2],obj@n)}
	else {a <- 1; b <- nrow(obj)}
	par(mar=c(1,4,1,1))
	plot(1, 1, xlab='', ylab='', xlim=c(a,b), ylim=c(0,1), xaxt="n", yaxt="n", bty="n", type="n")
	rect(seq(a,b-1), 0, seq(a+1,b), 1, col=pal[lbl[a:b]], border=NA)
	par(mar=c(4,4,1,1))
	plims <- getPropLims(obj, a, b)
	plot(1, 1, xlab='lon', ylab='lat', xlim=plims$x, ylim=plims$y, type="n", asp=1)
	points(obj$lon[a:b], obj$lat[a:b], col=pal[lbl[a:b]], pch=20)
	lines(obj$lon[a:b], obj$lat[a:b], col=pal[length(pal)])
 	par(parDef)
	})

#' @title Generate a pointwise .kml file of a \link{binClstPath_instance}
#'
#' @description \code{pkml} generates a pointwise KML file of a
#'   \link{binClstPath_instance}, which can be viewed using Google Earth or
#'   other GIS software.
#'
#' @param obj A \link{binClstPath_instance}.
#'
#' @param folder A character string indicating the name of the folder in which
#'   the .kml file will be saved. If the folder does not exist it is
#'   automatically created, (defaults to '~/embcDocs').
#'
#' @param markerRadius A numeric value indicating the radius of the markers to
#'   be plotted, (defaults to 5 pixels).
#'
#' @param display A boolean value (defaults to FALSE) to automatically launch
#'   Google-Earth from within R to visualize the generated .kml document.
#'   (Google Earth must already be installed on the system. In Windows, it must
#'   be associated with the .kml file type.)
#'
#' @param showClst A numeric vector indicating a subset of clusters to be shown.
#'
#' @param ... Parameters \code{folder}, \code{markerRadius}, \code{display} and
#'   \code{showClst} are optional.
#'
#' @return The path/name of the saved kml file.
#'
#' @seealso \code{\link{bkml}}
#'
#' @export
#' @rdname pkml
#'
#' @examples
#' \dontrun{
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth,info=-1)
#' # -- generate a pointwise .kml of the output --
#' pkml(mybcp)
#' # -- show only stopovers and automatically display the .kml document --
#' pkml(mybcp,showClst=c(1,2),display=TRUE)
#' }

setGeneric("pkml",
	function(obj,folder='embcDocs',markerRadius=15,display=FALSE,...)
	{standardGeneric("pkml")})

#' @rdname pkml
setMethod("pkml",signature(obj="binClstPath"),
	function(obj,folder,markerRadius,display,showClst=numeric(),...){
		if (length(showClst)>0) obj <- getSubPth(obj,showClst)
		kmlD <- new('kmlDoc',obj,bwise=FALSE,folder,markerRadius)
		if(display){
			if(.Platform$OS.type=='windows')
				shell.exec(kmlD@fName)
			else{
			  google <- c('/usr/bin/google-earth','/usr/bin/googleearth')
				system(paste(google[which(file.exists(google))],kmlD@fName,'&',sep=" "),ignore.stdout=TRUE,ignore.stderr=TRUE)
			  }
		}
	return(cat('saved file ',kmlD@fName,'\n'))
	})


#' @title Generate a burstwise .kml file of a \link{binClstPath_instance}.
#'
#' @description \code{bkml} generates a burstwise .kml file of a
#'   \link{binClstPath_instance}, which can be viewed using Google Earth or
#'   other GIS software. At first issue, this command can take some time because
#'   bursted segmentation has to be computed.
#'
#' @param obj A \link{binClstPath_instance}.
#'
#' @param folder A character string indicating the name of the folder in which
#'   the .kml file will be saved. If the folder does not exist it is
#'   automatically created, (defaults to '~/embcDocs').
#'
#' @param markerRadius A numeric value indicating the radius of the markers to
#'   be plotted, (defaults to 5 pixels).
#'
#' @param display A boolean value (defaults to FALSE) to automatically launch
#'   Google-Earth from within R to visualize the generated .kml document.
#'   (Google Earth must already be installed on the system. In Windows, it must
#'   be associated with the .kml file type.)
#'
#' @return The path/name of the saved kml file.
#'
#' @seealso \code{\link{pkml}}
#'
#' @export
#' @rdname bkml
#'
#' @examples
#' \dontrun{
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth,info=-1)
#' # -- generate a burstwise kml of the output --
#' bkml(mybcp)
#' }

setGeneric("bkml",
	function(obj,folder='embcDocs',markerRadius=15,display=FALSE)
	{standardGeneric("bkml")})

#' @rdname bkml
setMethod("bkml",signature(obj="binClstPath"),
	function(obj,folder,markerRadius,display){
		if (!obj@bursted) obj <- burstIt(obj,deparse(substitute(obj)))
		kmlD <- new('kmlDoc',obj,bwise=TRUE,folder,markerRadius)
		if(display){
			if(.Platform$OS.type=='windows')
				shell.exec(kmlD@fName)
		  else{
		    google <- c('/usr/bin/google-earth','/usr/bin/googleearth')
		    system(paste(google[which(file.exists(google))],kmlD@fName,'&',sep=" "),ignore.stdout=TRUE,ignore.stderr=TRUE)
		  }
		}
	return(cat('saved file ',kmlD@fName,'\n'))
	})

#' @title Generate an HTML pointwise webmap of a \link{binClstPath_instance}.
#'
#' @description \code{pmap} generates a pointwise .html file-map of a
#'   \link{binClstPath_instance} in HTML5, using Google Maps JavaScript API v3
#'   (\url{https://developers.google.com/maps/documentation/javascript/}). The
#'   resulting file can be viewed locally in most browsers (an internet
#'   connection is required for displaying the map tiles) or posted online.
#'
#' @param obj A \link{binClstPath_instance}.
#'
#' @param folder A character string indicating the name of the folder in which
#'   the .html file will be saved. If the folder does not exist it is
#'   automatically created, (defaults to '~/embcDocs').
#'
#' @param apiKey A character string specifying the API Key to be passed to the
#'   Google Maps server. No Key is needed for using Google Maps JavaScript API
#'   v3, but users may wish to specify a key in order to monitor web traffic if
#'   the document is being posted online.
#'
#' @param mapType A character string specifying the type of map to be used in
#'   the background. This value is passed directly to the Google Maps server,
#'   and currently can be set to ROADMAP, SATELLITE, HYBRID, or TERRAIN. (See
#'   the Google Maps API documentation for more information.)
#'
#' @param markerRadius A numeric value indicating the radius of the markers to
#'   be plotted, (defaults to 5 pixels).
#'
#' @param display A boolean value (defaults to FALSE) to automatically launch
#'   the system's default browser from witin R to visualize the generated .html
#'   document.
#'
#' @return The path/name of the saved html file.
#'
#' @export
#' @rdname pmap
#'
#' @examples
#' \dontrun{
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth,info=-1)
#' # -- generate a pointwise HTML of the output --
#' pmap(mybcp)
#' }

setGeneric("pmap",
           function(obj,folder='embcDocs',apiKey="",mapType='SATELLITE',markerRadius=15,display=FALSE)
           {standardGeneric("pmap")})

#' @rdname pmap
setMethod("pmap",signature(obj="binClstPath"),
          function(obj,folder,apiKey,mapType,markerRadius,display){
            wMap <- new('webMap',obj,bwise=FALSE,folder,apiKey,mapType,markerRadius)
            if(display) browseURL(wMap@fName)
            return(cat('saved file ',wMap@fName,'\n'))})

#' @title Generate an HTML burstwise webmap of a \link{binClstPath_instance}.
#'
#' @description \code{bmap} generates a burstwise .html file map of a
#'   \link{binClstPath_instance} in HTML5, using Google Maps JavaScript API v3
#'   (\url{https://developers.google.com/maps/documentation/javascript/}). The
#'   resulting file can be viewed locally in most browsers (an internet
#'   connection is required for displaying the map tiles) or posted online.
#'
#' @inheritParams pmap
#'
#' @return The path/name of the saved .html file.
#'
#' @export
#' @rdname bmap
#'
#' @examples
#' \dontrun{
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth,info=-1)
#' # -- generate a burstwise HTML of the output --
#' bmap(mybcp)
#' }

setGeneric("bmap",
           function(obj,folder='embcDocs',apiKey="",mapType='SATELLITE',markerRadius=15,display=FALSE)
           {standardGeneric("bmap")})

#' @rdname bmap
setMethod("bmap",signature(obj="binClstPath"),
          function(obj,folder,apiKey,mapType,markerRadius,display){
            if (!obj@bursted) obj <- burstIt(obj,deparse(substitute(obj)))
            wMap <- new('webMap',obj,bwise=TRUE,folder,apiKey,mapType,markerRadius)
            if(display) browseURL(wMap@fName)
            return(cat('saved file ',wMap@fName,'\n'))})
