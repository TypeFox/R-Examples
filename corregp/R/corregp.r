#' Functions and Methods for Correspondence Regression
#'
#' This package provides functions and methods for performing correspondence regression, i.e. the correspondence analysis of the
#'   crosstabulation of a categorical variable Y in function of another one X, where X can in turn be made up of the combination of various
#'   categorical variables.  
#'   Consequently, correspondence regression can be used to analyze the effects for a polytomous or multinomial outcome variable.  
#'   The central function in the package is \code{\link{corregp}}, which enables methods for printing, summarizing and plotting the output.  
#'   Additionally, there are functions for computing confidence intervals, ellipses or 3D ellipsoids (by means of bootstrapping).  
#' @section Contents:
#' This package consists of the following datasets, functions, generics and methods:  
#' \subsection{Datasets}{
#'   \itemize{
#'     \item{\code{\link{HairEye}} }{Hair and eye color of statistics students (data frame).}
#'     \item{\code{\link{COMURE}} }{The use of linguistic variants in translations vs. non-translations and in six different registers.}
#'     \item{\code{\link{TSS}}} {The use of inflected or uninflected determiners in vernacular Belgian Dutch.}
#'   }
#' }
#' \subsection{Functions}{
#'   \itemize{
#'     \item{\code{\link{ci}} }{A helper function to compute confidence intervals on the basis of a numeric vector.}
#'     \item{\code{\link{confGet}} }{An internal function to retrieve the \code{conf} components in a \code{corregp} object.}
#'     \item{\code{\link{corregp}} }{The basic function to perform correspondence regression. Typically, one starts here, and then one uses \code{print}, \code{summary}, \code{screeplot} or \code{plot} methods.}
#'   }
#' }
#' \subsection{Generics}{
#'   \itemize{
#'     \item{\code{\link{cint}} }{Compute confidence intervals.}
#'     \item{\code{\link{cell}} }{Compute confidence ellipses.}
#'     \item{\code{\link{cell3d}} }{Compute 3D confidence ellipsoids.}
#'     \item{\code{\link{ciplot}} }{Plot confidence intervals.}
#'     \item{\code{\link{pcplot}} }{Plot parallel coordinates.}
#'   }
#' }
#' \subsection{Methods}{
#'   \itemize{
#'     \item{\code{\link{print.corregp}} }{Print the output of a correspondence regression.}
#'     \item{\code{\link{summary.corregp}} }{Give a \code{summary} of a correspondence regression.}
#'     \item{\code{\link{print.summary.corregp}} }{Print the \code{summary} of a correspondence regression.}
#'     \item{\code{\link{screeplot.corregp}} }{Make a scree plot on the basis of the output of a correspondence regression.}
#'     \item{\code{\link{cint.corregp}} }{Compute confidence intervals on the basis of the output of a correspondence regression. Typically, this function is not so much used directly as it is called by a \code{ciplot.corregp} command.}
#'     \item{\code{\link{ciplot.corregp}} }{Plot confidence intervals on the basis of the output of a correspondence regression.}
#'     \item{\code{\link{pcplot.corregp}} }{Make a parallel coordinate plot on the basis of the output of a correspondence regression.}
#'     \item{\code{\link{cell.corregp}} }{Compute confidence ellipses on the basis of the output of a correspondence regression. Typically, this function is not so much used directly as it is called by a \code{plot.corregp} command.}
#'     \item{\code{\link{plot.corregp}} }{Plot the output (and the optional confidence ellipses) of a correspondence regression.}
#'     \item{\code{\link{cell3d.corregp}} }{Compute 3D confidence ellipsoids on the basis of a correspondence regression. Typically, this function is not so much used directly as it is called by a \code{plot3d.corregp} command.}
#'     \item{\code{\link{plot3d.corregp}} }{Plot the 3D output (and the optional confidence ellipsoids) of a correspondence regression.}
#'   }
#' }
#' 
#' @section Future prospects:
#' \itemize{
#'   \item Include a vignette.
#'   \item Allow for more than one random factor (e.g. enable the argument \code{parm} to be a character vector of more than one element).
#'   \item Specify methods for \code{anova}, \code{predict}, \code{residuals} etc.
#'   \item Enable scale transformations for all plots (and corresponding confidence regions).
#'   \item Provide the possibility for so-called "calibration lines".
#' }
#'
#' @section Author:
#' Koen Plevoets, \email{koen.plevoets@@ugent.be}  
#'
#' @section Acknowledgements:
#' This package has benefited greatly from the helpful comments of Isabelle Delaere and Gert De Sutter. Thanks to Kurt Hornik and Uwe Ligges for proofing this package.  
#' @docType package
#' @name corregp-package
NULL

#' @import graphics
NULL

#' @import rgl
NULL

#' @import stats
NULL

#' @importFrom ellipse ellipse
NULL

#' @importFrom gplots barplot2
NULL

#' @importFrom gplots plotCI
NULL

#' Hair and Eye Color of Statistics Students (Data Frame)
#'
#' The distribution of hair color, eye color and sex among 592 statistics students (from Snee 1974 and Friendly 1992).  
#' @format A data frame with 592 rows and 3 variables.
#' @source This is simply a data frame version of the in-built data set \code{\link[datasets]{HairEyeColor}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' haireye.crg
#' summary(haireye.crg, parm = "b", add_ci = TRUE)
#' screeplot(haireye.crg, add_ci = TRUE)
#' plot(haireye.crg, x_ell = TRUE, xsub = c("Hair", "Sex"))
#' @docType data
#' @name HairEye
NULL

#' The Use of Linguistic Variants in Translations vs. Non-translations and in Six Different Registers
#'
#' This data set was a case study in the COMURE project ("\strong{co}rpus-based, \strong{mu}ltivariate research of \strong{re}gister variation in translated and
#'   non-translated Belgian Dutch") which was conducted at the Department of Translation, Interpreting and Communication of Ghent University between 2010 and 2014.  
#' @format A data frame with 3762 rows and 4 variables.
#' \itemize{
#'   \item{\code{Variant} }{The linguistic variant used in a set of alternatives (27 levels).}
#'   \item{\code{Variable} }{The linguistic variable specifying a set of alternatives (13 levels).}
#'   \item{\code{Register} }{The register or "Text type" of the data (6 levels).}
#'   \item{\code{Language} }{The language (and source language) of the data (3 levels).}
#' }
#' @source
#' Delaere, I., G. De Sutter and K. Plevoets (2012) Is translated language more standardized than non-translated language? \emph{Target} \strong{24} (2), 203--224.
#' @examples
#' \donttest{
#' data(COMURE)
#' # The execution of corregp may be slow, due to bootstrapping:  
#' comure.crg <- corregp(Variant ~ Register * Language, data = COMURE, part = "Variable", b = 3000)
#' comure.crg
#' summary(comure.crg, parm = "b", add_ci = TRUE)
#' screeplot(comure.crg, add_ci = TRUE)
#' comure.col <- ifelse( xtabs(~ Variant + Variety, data = COMURE)[, "Standard"] > 0, "blue", "red")
#' plot(comure.crg, x_ell = TRUE, xsub = c("Register", "Language"), col_btm = comure.col, 
#'   col_top = "black")
#' }
#' @docType data
#' @name COMURE
NULL

#' The Use of Inflected or Uninflected Determiners in the Belgian Dutch Vernacular
#'
#' The distribution of the Belgian Dutch \emph{-e(n)}-suffix with 14 determiners in 14 registers and for several speaker characteristics.  
#' @format A data frame with 40778 rows and 12 variables.
#' \itemize{
#'   \item{\code{Variant} }{The linguistic variant used in a set of alternatives (35 levels).}
#'   \item{\code{Variable} }{The linguistic variable specifying a set of alternatives (14 levels).}
#'   \item{\code{Inflected} }{Numeric variable specifying whether the linguistic variant is inflected (\code{1}) or not (\code{0}).}
#'   \item{\code{Register} }{The register of the data in the Spoken Dutch Corpus (14 levels: see 
#'     \href{http://lands.let.ru.nl/cgn/doc_English/topics/version_1.0/overview.htm}{here} for their definition).}
#'   \item{\code{SpeakerID} }{The ID of the speaker in the Spoken Dutch Corpus (1144 levels).}
#'   \item{\code{Region} }{The region in which the speaker lived until the age of 18 (4 levels).}
#'   \item{\code{Sex} }{The sex of the speaker (2 levels).}
#'   \item{\code{BirthYear} }{The year in which the speaker was born (63 levels).}
#'   \item{\code{Decade} }{The decade in which the speaker was born (7 levels).}
#'   \item{\code{Generation} }{The generation cohort in which the speaker was born (5 levels).}
#'   \item{\code{Education} }{The level of education of the speaker (3 levels).}
#'   \item{\code{Occupation} }{The level of occupation of the speaker (10 levels: see 
#'     \href{http://lands.let.ru.nl/cgn/doc_English/topics/version_1.0/metadata/speakers.htm}{here} for their definition).}
#' }
#' @source
#' Plevoets, K. (2008) \emph{Tussen spreek- en standaardtaal}. Leuven, Doctoral dissertation. Available online \href{https://lirias.kuleuven.be/handle/1979/1760}{here} and
#'   \href{http://statbel.fgov.be/nl/modules/digibib/bevolking/1431_tussen_spreek-_en_standaardtaal.jsp}{here}.
#' @examples
#' \donttest{
#' data(TSS)
#' # The execution of corregp may be slow, due to bootstrapping:  
#' tss.crg <- corregp(Variant ~ Register * Region, data = TSS, part = "Variable", b = 3000)
#' tss.crg
#' summary(tss.crg, parm = "b", add_ci = TRUE)
#' screeplot(tss.crg, add_ci = TRUE)
#' tss.col <- ifelse( xtabs(~ Variant + Inflected, data = TSS)[, 1] > 0, "blue", "red")
#' plot(tss.crg, x_ell = TRUE, xsub = c("Register", "Region"), col_btm = tss.col, col_top = "black")
#' }
#' @docType data
#' @name TSS
NULL

#' Correspondence Regression
#'
#' This is the basic function for \emph{correspondence regression}, i.e. the correspondence analysis of a contingency table formed
#'   by the categorical variables Y and X, where X can be in turn made up of the combinations of various categorical variables.  
#' @param formula A \code{\link[stats]{formula}} specification of which factors to cross with each other. The left-hand (\code{y}) side must be a single factor.
#'   The right-hand side (\code{x}) can involve all the usual specifications of interactions and/or nested analyses.
#' @param data The data frame containing the variables specified in the \code{formula}.
#' @param part The name of a factor partitioning the levels of the left-hand side \code{y} into groups.
#'   This argument is relevant for analyses in which one wants to remove between-item variation.
#' @param b Number of the bootstrap replications. If \code{0} (i.e. the default), then the analysis is exploratory.
#' @param xep Logical specifying whether to output the separate terms in the right-hand side (\code{x}) as components in a list.
#'   If \code{FALSE}, then all \code{x} output is collected in a matrix.
#' @param std Logical specifying whether to output the standardized coordinates. Defaults to \code{FALSE}.
#' @param rel Logical specifying whether to divide the coordinates by the \code{sqrt} of their totals, so that one obtains coordinates for
#'   the relative frequencies (as is customary in correspondence analysis). Defaults to \code{TRUE}.
#' @param phi Logical specifying whether to compute the output on the scale of the \emph{Chi-squared} value of the contingency table or of the \emph{Phi-squared} value
#'   (which is \emph{Chi-squared} divided by \emph{N}). Reminiscent of \code{\link[MASS]{corresp}} in package \pkg{MASS}, defaults to \code{FALSE}.
#' @param chr Character specifying the separator string for constructing the interaction terms.
#' @details
#' Correspondence regression rests on the idea, described by Van der Heijden et al. (1989) and quoted in Greenacre (2007: 272), of using correspondence analysis to inspect
#'   the interactions in a log-linear analysis.  More specifically, as log-linear analysis or Poisson regression is sometimes used to model a polytomous or multinomial response
#'   variable (in a GLM), correspondence regression enables the analysis of a categorical factor (\code{Y}) in terms of other (possibly interacting) factors (\code{X}).  These
#'   are specified in the argument \code{formula}, which can be constructed in all the usual ways of specifying a model formula: e.g. \code{Y ~ X1 * X2} as a shorthand for
#'   \code{Y ~ X1 + X2 + X1 : X2}, or \code{Y ~ X1 * X2 - X1 : X2}, \code{Y ~ (X1 + X2 + X3) ^ 2}, etc.  Correspondence regression then crosstabulates the \code{Y} factor with all the
#'   combinations in \code{X}, thus producing a typical contingency table, on which a simple correspondence analysis is performed (see Greenacre 2007: 121-128 for the outline of
#'   this approach). The more general effects in \code{X} are obtained by aggregating the combinations.   
#'
#' Correspondence regression also allows for inferential validation of the effects, which is done by means of the bootstrap. Setting the argument \code{b} to a number \eqn{> 0}, \code{b}
#'   replicates of the contingency table are generated with multinomial sampling. From these, \code{b} new values are derived for the coordinates in both \code{Y} and
#'   \code{X} as well as for the eigenvalues (also called the "principal inertias").  On the basis of the replicate values, confidence intervals, ellipses or ellipsoids can
#'   be computed.  CAUTION: bootstrapping is computationally quite intensive, so it can take a while to reach results, especially with a large \code{b}.  
#'
#' The argument \code{parm} can be used when the levels of \code{Y} are grouped/partitioned/nested into clusters and one wants to exclude the heterogeneity between the clusters.  Thus,
#'   \code{parm} is equivalent to a \emph{random factor}, although \code{corregp} currently allows for only one such factor.  The use of \code{parm} can be relevant for so-called
#'   \emph{lectometric} analyses in linguistics.  
#' @return An object of class "corregp", i.e. a list with components:
#' \item{\code{eigen} }{A vector of eigenvalues of the correpondence regression.}
#' \item{\code{y} }{The coordinates (matrix) of the Y levels.}
#' \item{\code{x} }{The coordinates of the X levels. If \code{xep} is \code{TRUE}, then this is a list with a component for each term name.}
#' \item{\code{freq} }{A list of the frequencies of every Y and X level.}
#' \item{\code{conf} }{If \eqn{b>0}. A list of bootstrap replicates for the eigenvalues, the Y levels and the X levels.}
#' \item{\code{aux} }{A list of auxiliary information (such as the U and V matrices of the SVD, the specified values for all the arguments) to be passed to other functions and methods.}
#' @references
#' Greenacre, M. (2007) \emph{Correspondence analysis in practice, Second edition}. Boca Raton: Chapman and Hall/CRC.
#'
#' Van der Heijden, P.G.M., A. de Falguerolles and J. de Leeuw (1989) A combined approach to contingency table analysis using correspondence analysis and log-linear analysis.
#'   \emph{Applied Statistics} \strong{38} (2), 249--292.
#' @seealso \code{\link{print.corregp}}, \code{\link{summary.corregp}}, \code{\link{screeplot.corregp}}, \code{\link{plot.corregp}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' haireye.crg
#' @export
corregp <- function(formula,data,part=NULL,b=0,xep=TRUE,std=FALSE,rel=TRUE,phi=FALSE,chr=".") {
	n.tot <- nrow(data)
	f.var <- all.vars(formula,functions=FALSE,unique=TRUE)
	if (length(part)>1) {
		stop("argument part can currently only contain one element",call.=FALSE)
		}
	part <- ifelse(is.null(part),NA,part)
	if (is.na(part)) {
		data <- cbind(data.frame(lapply(data[,c(f.var[-1],f.var[1])],as.factor)),part=as.factor(rep("FULL",times=n.tot)))
		}
	else {
		data <- data.frame(lapply(data[,c(f.var[-1],f.var[1],part)],as.factor))
		}
	f.dim <- as.vector(sapply(data,nlevels))
	f.seq <- 1:length(f.dim)
	y.num <- match(f.var[1],table=names(data))
	f.num <- 1:y.num
	p.num <- lapply(strsplit(labels(terms(formula)),split=":"),match,table=names(data))
	p.lev <- lapply(p.num,function(x){levels(interaction(data[,x],sep=chr,drop=FALSE))})
	p.seq <- 1:length(p.num)
	p.lst <- cumsum(sapply(p.lev,length))
	p.fst <- c(0,p.lst[-length(p.lst)]) + 1
	x.len <- p.lst[length(p.lst)]
	y.len <- f.dim[y.num]
	n.rot <- ifelse(phi,1,sqrt(n.tot))
	r <- 0
	while (r<=b) {
		if (r>0) {
			r.tab <- array(rmultinom(1,size=n.tot,prob=f.req),dim=f.dim)
			}
		else {
			r.tab <- array(table(data),dim=f.dim)
			if (b>0) {
				f.req <- as.vector(r.tab)
				}
			}
		y.tot <- apply(r.tab,c(y.num,f.seq[-f.num]),sum)
		x.tot <- apply(r.tab,f.seq[-y.num],sum)
		r.dev <- apply(r.tab-sweep(sweep(aperm(array(x.tot,dim=c(f.dim[-y.num],f.dim[y.num])),perm=c(f.seq[-y.num],y.num)),c(y.num,f.seq[-f.num]),y.tot,"*"),f.seq[-f.num],apply(x.tot,y.num,sum),"/"),f.num,sum)
		r.den <- sweep(array(apply(x.tot,f.num[-y.num],sum),dim=f.dim[f.num]),y.num,y.tot <- apply(y.tot,1,sum),"*")
		r.tab <- array(n.rot * r.dev / sqrt(r.den),dim=c(prod(f.dim[f.num[-y.num]]),y.len))
		r.tab[is.nan(r.tab)] <- 0
		if (r==0) {
			r.svd <- svd(r.tab)
			u.svd <- r.svd$u
			v.svd <- r.svd$v
			s.svd <- r.svd$d
			r.ank <- length(s.svd)
			x.frq <- numeric(length=x.len)
			}
		x.loc <- matrix(nrow=x.len,ncol=r.ank)
		for (x1 in p.seq) {
			p.ind <- p.num[[x1]]
			p.dim <- c(p.ind,y.num)
			p.dev <- apply(r.dev,p.dim,sum)
			p.den <- sqrt(1 / apply(r.den,p.dim,sum))
			p.den[is.infinite(p.den)] <- 0
			p.tot <- as.vector(apply(x.tot,p.ind,sum))
			if (r==0) {
				x.frq[seq(p.fst[x1],p.lst[x1],by=1)] <- p.tot
				}
			if (rel) {
				p.wgt <- sqrt(1 / p.tot)
				}
			else {
				p.wgt <- rep(1,times=prod(f.dim[p.ind]))
				}
			p.wgt[p.tot==0] <- NA
			x.loc[seq(p.fst[x1],p.lst[x1],by=1),] <- n.rot * sweep( array(p.dev * p.den,dim=c(prod(f.dim[p.ind]),y.len)) %*% v.svd %*% diag(1/s.svd) ,1,p.wgt,"*")
			}
		if (rel) {
			y.wgt <- sqrt(1 / y.tot)
			}
		else {
			y.wgt <- rep(1,times=y.len)
			}
		y.wgt[y.tot==0] <- NA
		y.loc <- sweep ( t(r.tab) %*% u.svd %*% diag(1/s.svd) ,1,y.wgt,"*")
		if (!std) {
			y.loc <- y.loc %*% diag(s.svd)
			x.loc <- x.loc %*% diag(s.svd)
			}
		if (r>0) {
			r.con[r,] <- diag(t(u.svd) %*% r.tab %*% v.svd)^2
			for (yl in y.seq) {
				y.con[[yl]][r,] <- y.loc[yl,]
				}
			for (xl in x.seq) {
				x.con[[xl]][r,] <- x.loc[xl,]
				}
			}
		else {
			f.out <- list(eigen=s.svd^2,y=y.loc,x=x.loc,freq=list(y=y.tot,x=x.frq))
			if (b>0) {
				y.seq <- 1:y.len
				x.seq <- 1:x.len
				r.con <- matrix(0,nrow=b,ncol=r.ank)
				y.con <- rep(list(r.con),times=y.len)
				x.con <- rep(list(r.con),times=x.len)
				}
			}
		r <- r + 1
		}
	r.seq <- 1:r.ank
	dimnames(f.out$y) <- list(y.lab <- levels(data[,y.num]),r.seq)
	dimnames(f.out$x) <- list(x.lab <- unlist(p.lev),r.seq)
	names(f.out$freq$y) <- y.lab
	names(f.out$freq$x) <- x.lab
	if (b>0) {
		dimnames(r.con) <- list(1:b,r.seq)
		names(y.con) <- y.lab
		names(x.con) <- x.lab
		f.out <- c(f.out,list(conf=list(eigen=r.con,y=y.con,x=x.con)))
		}
	xep <- ifelse(length(f.var)==2,FALSE,xep)
	if (xep) {
		f.out$x <- lapply(p.seq,function(x2){f.out$x[seq(p.fst[x2],p.lst[x2],by=1),]})
		names(f.out$x) <- p.lab <- sapply(p.num,function(x3){paste(names(data)[x3],sep="",collapse=chr)})
		f.out$freq$x <- lapply(p.seq,function(x4){f.out$freq$x[seq(p.fst[x4],p.lst[x4],by=1)]})
		names(f.out$freq$x) <- p.lab
		if (b>0) {
			f.out$conf$x <- lapply(p.seq,function(x5){f.out$conf$x[seq(p.fst[x5],p.lst[x5],by=1)]})
			names(f.out$conf$x) <- p.lab
			}
		}
	colnames(u.svd) <- r.seq
	colnames(v.svd) <- r.seq
	f.out$aux <- list(U=u.svd,V=v.svd,formula=formula,data=as.list(match.call())$data,part=part,std=std,rel=rel,phi=phi,chr=chr)
	class(f.out) <- "corregp"
	return(f.out)
	}

#' Printing Correspondence Regression
#'
#' Method to print the output of \code{\link{corregp}}.  
#' @param x The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param nf The number of dimensions to print. Defaults to the first two dimensions.
#' @param ... Further arguments passed to or from other methods.
#' @return The output of a call to \code{\link{corregp}}.
#' @seealso \code{\link{corregp}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' haireye.crg
#' print(haireye.crg, nf = 3)
#' @export
print.corregp <- function(x,nf=2,...) {
	crg <- x
	nf <- ifelse(is.null(nf),length(crg$eigen),nf)
	nf <- ifelse(is.character(nf),match(nf,table=colnames(crg$y)),nf)
	f.var <- all.vars(crg$aux$formula,functions=FALSE,unique=TRUE)
	cat("Correspondence regression of ",format(crg$aux$formula)," in ",crg$aux$data,"\n",sep="",fill=FALSE)
	if (!is.na(crg$aux$part)) {
		cat("Partitioning factor: ",crg$aux$part,"\n",sep="",fill=FALSE)
		}
	cat("\nEigenvalues: ",crg$eigen[1:nf],"\n",sep=" ",fill=FALSE)
	cat("\nY (",f.var[1],"):\n",sep="",fill=FALSE)
	print(crg$y[,1:nf],...)
	if (is.list(crg$x)) {
		p.var <- names(crg$x)
		cat("\nX:\n",sep="",fill=FALSE)
		for (p in 1:length(crg$x)) {
			cat(p.var[p],":\n",sep="",fill=FALSE)
			print(crg$x[[p]][,1:nf],...)
			cat("\n",sep="",fill=FALSE)
			}
		}
	else {
		cat("\nX",ifelse(length(f.var)==2,paste(" (",f.var[2],")",sep="",collapse=NULL),""),":\n",sep="",fill=FALSE)
		print(crg$x[,1:nf],...)
		cat("\n",sep="",fill=FALSE)
		}
	invisible()
	}

#' Summarizing Correspondence Regression
#'
#' Method to produce a summary of a correspondence regression.  
#' @param object The outout of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param parm The parameter for which to compute the contributions \code{contrib}. Can be either \code{"y"} for the Y contributions, \code{"x"} for the X contributions, \code{"both"}
#'   which can be abbreviated to \code{"b"}, or a vector of term names in X. Defaults to \code{"b"}.
#' @param contrib The type of contributions to be computed: either \emph{from points to axes} (absolute contributions) or \emph{from axes to points} (squared correlations).
#'   The specification can be \code{"pnts_to_axes"} or \code{"axes_to_pnts"}, \code{"pts2axs"} or \code{"axs2pts"}, \code{"p_a"} or \code{"a_p"}, or any other reasonable abbreviation.
#' @param nf The number of dimensions to be retained in the reduced space. Defaults to all dimenions (no reduction).
#' @param add_ci Logical specifying whether to compute confidence intervals for the eigenvalues (and eigenvalues only). Defaults to \code{FALSE}.
#' @param cl The confidence level for the confidence intervals. Defaults to \code{0.95}.
#' @param nq Logical specifying whether to use a normal quantile (i.e. apply \code{\link[stats]{qnorm}}) in the computation of the confidence intervals.
#'   Defaults to \code{TRUE}. If \code{FALSE}, then the confidence intervals are computed directly with the \code{\link[stats]{quantile}} function.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class "summary.corregp", providing a summary of a correspondence regression, i.e. a list with components:
#' \item{\code{formula} }{The \code{\link[stats]{formula}} specified to the \code{formula} argument in the call to \code{corregp}.}
#' \item{\code{data} }{The name of the data frame specified to the \code{data} argument in the call to \code{corregp}.}
#' \item{\code{part} }{The name of the factor specified to the\code{part} argument in the call to \code{corregp}.}
#' \item{\code{chi_squared} }{The chi-squared value of the correspondence regression.}
#' \item{\code{phi_squared} }{The phi-squared value of the correspondence regression, i.e. the chi-squared value divided by \code{N}.}
#' \item{\code{N} }{The total number of observations.}
#' \item{\code{eigen} }{Depending on \code{add_ci}: if \code{FALSE}, a matrix of the actual eigenvalues, their percentages and cumulative percentages; if \code{TRUE}, a list of the actual eigenvalues, their percentages and cumulative percentages together with the lower and upper confidence limits for each.}
#' \item{\code{y} }{If \code{parm} is \code{"y"} or \code{"b"}. A list of components \code{p_a} for the absolute contributions and//or \code{a_p} for the squared correlations, depending on \code{contrib}.}
#' \item{\code{x} }{If \code{parm} is \code{"y"}, \code{"b"} or any of the term names in X. A list of components \code{p_a} for the absolute contributions and/or \code{a_p} for the squared correlations, depending in \code{contrib}.}
#' @seealso \code{\link{corregp}}, \code{\link{print.summary.corregp}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' summary(haireye.crg, add_ci = TRUE)
#' summary(haireye.crg, parm = "y", contrib = "pts_axs", nf = 2)
#' @export
summary.corregp <- function(object,parm=NULL,contrib=NULL,nf=NULL,add_ci=FALSE,cl=0.95,nq=TRUE,...) {
	crg <- object
	if (is.null(parm) && !is.null(contrib)) {
		parm <- "b"
		}
	if (!is.null(parm) && is.null(contrib)) {
		contrib <- "b"
		}
	parm <- ifelse(is.null(parm),NA,parm)
	contrib <- tolower(contrib)
	contrib <- ifelse(contrib %in% c("p_a","pts_axs","pts2axs","ptstoaxs","pts_to_axs","pnts_axes","pnts2axes","pntstoaxes","pnts_to_axes"),"p_a",contrib)
	contrib <- ifelse(contrib %in% c("a_p","axs_pts","axs2pts","axstopts","axs_to_pts","axes_pnts","axes2pnts","axestopnts","axes_to_pnts"),"a_p",contrib)
	nf <- ifelse(is.null(nf),length(crg$eigen),nf)
	nf <- ifelse(is.character(nf),match(nf,table=colnames(crg$y)),nf)
	o.eig <- crg$eigen
	o.tot <- sum(crg$freq$y)
	o.phi <- ifelse(crg$aux$phi,o.tot,1)
	s.out <- list(formula=crg$aux$formula,data=crg$aux$data,part=crg$aux$part,chi_squared=sum(o.eig)*o.phi,phi_squared=sum(o.eig)*o.phi/o.tot,N=o.tot)
	o.lab <- colnames(crg$y[,1:nf])
	if (add_ci && is.null(crg$conf)) {
		add_ci <- FALSE
		warning(paste("no bootstrapping was done in",as.list(match.call())$object,sep=" ",collapse=NULL),call.=FALSE)
		}
	if (add_ci) {
		s.out$eigen$value <- as.matrix(rbind(o.eig[1:nf],apply(crg$conf$eigen[,1:nf],2,ci,cl=cl,nq=nq)))
		s.out$eigen$"%" <- as.matrix(rbind(o.eig[1:nf]/sum(o.eig),apply(sweep(crg$conf$eigen[,1:nf],1,apply(crg$conf$eigen,1,sum),"/"),2,ci,cl=cl,nq=nq)))
		s.out$eigen$"cum_%" <- as.matrix(rbind(cumsum(o.eig[1:nf])/sum(o.eig),apply(sweep(t(apply(crg$conf$eigen[,1:nf],1,cumsum)),1,apply(crg$conf$eigen,1,sum),"/"),2,ci,cl=cl,nq=nq)))
		dimnames(s.out$eigen$value) <- list(c("","lower","upper"),o.lab)
		dimnames(s.out$eigen$"%") <- list(c("","lower","upper"),o.lab)
		dimnames(s.out$eigen$"cum_%") <- list(c("","lower","upper"),o.lab)
		}
	else {
		s.out$eigen <- as.matrix(rbind(o.eig[1:nf],o.eig[1:nf]/sum(o.eig),cumsum(o.eig[1:nf])/sum(o.eig)))
		dimnames(s.out$eigen) <- list(c("value","%","cum_%"),o.lab)
		}
	if (!is.na(parm)) {
		if (crg$aux$std) {
			o.std <- rep(1,length(o.eig))
			}
		else {
			o.std <- o.eig
			}
		}
	if (parm %in% c("y","b","both")) {
		if (contrib %in% c("p_a","b","both")) {
			if (crg$aux$rel) {
				y.rel <- crg$freq$y
				}
			else {
				y.rel <- rep(1,nrow(crg$y))
				}
			s.out$y$p_a <- sweep(sweep(crg$y[,1:nf]^2,1,y.rel,"*"),2,o.std[1:nf],"/")
			}
		if (contrib %in% c("a_p","b","both")) {
			s.out$y$a_p <- sweep(sweep(crg$y[,1:nf]^2,2,(o.eig/o.std)[1:nf],"*"),1,apply(sweep(crg$y^2,2,(o.eig/o.std),"*"),1,sum),"/")
			}
		}
	if (!is.list(crg$x) && parm %in% c("x","b","both")) {
		if (contrib %in% c("p_a","b","both")) {
			if (crg$aux$rel) {
				x.rel <- crg$freq$x
				}
			else {
				x.rel <- rep(1,nrow(crg$x))
				}
			s.out$x$p_a <- sweep(sweep(crg$x[,1:nf]^2,1,x.rel,"*"),2,o.std[1:nf],"/")
			}
		if (contrib %in% c("a_p","b","both")) {
			s.out$x$a_p <- sweep(sweep(crg$x[,1:nf]^2,2,(o.eig/o.std)[1:nf],"*"),1,apply(sweep(crg$x^2,2,(o.eig/o.std),"*"),1,sum),"/")
			}
		}
	if (is.list(crg$x) && all(parm %in% c("x","b","both",names(crg$x)))) {
		if (parm %in% c("x","b","both")) {
			parm <- names(crg$x)
			}
		if (contrib %in% c("p_a","b","both")) {
			if (crg$aux$rel) {
				x.rel <- crg$freq$x[parm]
				}
			else {
				x.rel <- lapply(crg$x[parm],function(x1){rep(1,nrow(x1))})
				}
			s.out$x$p_a <- lapply(parm,function(x2){sweep(sweep(crg$x[[x2]][,1:nf]^2,1,x.rel[[x2]],"*"),2,o.std[1:nf],"/")})
			names(s.out$x$p_a) <- parm
			}
		if (contrib %in% c("a_p","b","both")) {
			s.out$x$a_p <- lapply(parm,function(x3){sweep(sweep(crg$x[[x3]][,1:nf]^2,2,(o.eig/o.std)[1:nf],"*"),1,apply(sweep(crg$x[[x3]]^2,2,(o.eig/o.std),"*"),1,sum),"/")})
			names(s.out$x$a_p) <- parm
			}
		}
	class(s.out) <- "summary.corregp"
	return(s.out)
	}

#' Printing the Summary of Correspondence Regression
#'
#' Method to print the output of \code{\link{summary.corregp}}.  
#' @param x The output of a call to \code{summary} on a "corregp" object (i.e. an object of class "summary.corregp").
#' @param ... Further arguments passed to or from other methods.
#' @return The output of a call to \code{summary} on a "corregp" object. The eigenvalues and contributions are printed with \code{TOTAL}s.
#' @seealso \code{\link{summary.corregp}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' summary(haireye.crg, add_ci = TRUE)
#' summary(haireye.crg, parm = "y", contrib = "pts_axs", nf = 2)
#' @export
print.summary.corregp <- function(x,...) {
	crs <- x
	f.var <- all.vars(crs$formula,functions=FALSE,unique=TRUE)
	cat("Summary of correspondence regression of ",format(crs$formula)," in ",crs$data,"\n",sep="",fill=FALSE)
	if (!is.na(crs$part)) {
		cat("Partitioning factor: ",crs$part,"\n",sep=" ",fill=FALSE)
		}
	cat("\nChi-squared: ",crs$chi_squared,"\nPhi-squared: ",crs$phi_squared,"\nN: ",crs$N,"\n",sep="",fill=FALSE)
	cat("\n\nEigenvalues:\n",sep="",fill=FALSE)
	if (is.list(crs$eigen)) {
		cat(" Value:\n",sep="",fill=FALSE)
		print(cbind(crs$eigen$value,TOTAL=c(sum(crs$eigen$value[1,]),NA,NA)),na.print="",...)
		cat("\n Percentage (%):\n",sep="",fill=FALSE)
		print(cbind(crs$eigen$"%",TOTAL=c(sum(crs$eigen$"%"[1,]),NA,NA)),na.print="",...)
		cat("\n Cumulative percentage (cum_%):\n",sep="",fill=FALSE)
		print(cbind(crs$eigen$"cum_%",TOTAL=c(sum(crs$eigen$"%"[1,]),NA,NA)),na.print="",...)
		cat("\n",sep="",fill=FALSE)
		}
	else {
		print(cbind(crs$eigen,TOTAL=c(apply(crs$eigen[1:2,],1,sum),sum(crs$eigen[2,]))),...)
		cat("\n",sep="",fill=FALSE)
		}
	if ("y" %in% names(crs) || "x" %in% names(crs)) {
		cat("\nContributions:\n\n",sep="",fill=FALSE)
		}
	if ("y" %in% names(crs)) {
		cat("Y (",f.var[1],"):\n",sep="",fill=FALSE)
		if ("p_a" %in% names(crs$y)) {
			cat(" Points to axes (Absolute contributions):\n",sep="",fill=FALSE)
			print(rbind(crs$y$p_a,TOTAL=apply(crs$y$p_a,2,sum)),...)
			cat("\n",sep="",fill=FALSE)
			}
		if ("a_p" %in% names(crs$y)) {
			cat(" Axes to points (Squared correlations):\n",sep="",fill=FALSE)
			print(cbind(crs$y$a_p,TOTAL=apply(crs$y$a_p,1,sum)),...)
			cat("\n",sep="",fill=FALSE)
			}
		}
	if ("x" %in% names(crs)) {
		cat("X",ifelse(length(f.var)==2,paste(" (",f.var[2],")",sep="",collapse=NULL),""),":\n",sep="",fill=FALSE)
		if ("p_a" %in% names(crs$x)) {
			cat(" Points to axes (Absolute contributions):\n",sep="",fill=FALSE)
			if (is.list(crs$x$p_a)) {
				p.var <- names(crs$x$p_a)
				for (p1 in 1:length(p.var)) {
					cat("\n   ",p.var[p1],":\n",sep="",fill=FALSE)
					print(rbind(crs$x$p_a[[p1]],TOTAL=apply(crs$x$p_a[[p1]][complete.cases(crs$x$p_a[[p1]]),],2,sum)),...)
					}
				cat("\n",sep="",fill=FALSE)
				}
			else {
				print(rbind(crs$x$p_a,TOTAL=apply(crs$x$p_a[complete.cases(crs$x$p_a),],2,sum)),...)
				cat("\n",sep="",fill=FALSE)
				}
			}
		if ("a_p" %in% names(crs$x)) {
			cat(" Axes to points (Squared correlations):\n",sep="",fill=FALSE)
			if (is.list(crs$x$a_p)) {
				p.var <- names(crs$x$a_p)
				for (p2 in 1:length(p.var)) {
					cat("\n   ",p.var[p2],":\n",sep="",fill=FALSE)
					print(cbind(crs$x$a_p[[p2]],TOTAL=apply(crs$x$a_p[[p2]],1,sum)),...)
					}
				cat("\n",sep="",fill=FALSE)
				}
			else {
				print(cbind(crs$x$a_p,TOTAL=apply(crs$x$a_p,1,sum)),...)
				cat("\n",sep="",fill=FALSE)
				}
			}
		}
	invisible()
	}

#' Confidence Interval
#'
#' This is the basic function for computing a confidence interval on the basis of a sample of data values.  
#' @param x A numeric vector.
#' @param cl The confidence level for the confidence interval. Defaults to \code{0.95}.
#' @param nq Logical specifying whether to use a normal quantile (i.e. apply \code{\link[stats]{qnorm}}) in the computation of the confidence interval.
#'   Defaults to \code{TRUE}. If \code{FALSE}, then the confidence interval is computed directly with the \code{\link[stats]{quantile}} function.
#' @return A vector with two components \code{Lower} and \code{Upper} giving the lower and upper confidence limits respectively.
#' @seealso \code{\link{ciplot.corregp}}, \code{\link[stats]{confint}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' ci(haireye.crg$conf$eigen[, 1])
#' ci(haireye.crg$conf$eigen[, 2])
#' @export
ci <- function(x,cl=0.95,nq=TRUE) {
	x <- x[complete.cases(x)]
	if (nq) {
		avg <- mean(x)
		s.e <- sd(x)*abs(qnorm((1-cl)/2))
		out <- c(Lower=avg-s.e,Upper=avg+s.e)
		}
	else {
		out <- c(Lower=quantile(x,probs=(1-cl)/2,names=FALSE),Upper=quantile(x,probs=1-(1-cl)/2,names=FALSE))
		}
	return(ifelse(!is.na(out),out,c(Lower=NA,Upper=NA)))
	}

#' Scree Plotting
#'
#' Method to produce a \emph{scree plot}, i.e. a bar chart of the eigenvalues.  
#' @param x The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param type A character specification  of which type of values to plot: either \code{"value"} for the \emph{actual eigenvalues}, \code{"\%"} for \emph{percentages} or
#'   \code{"cum_\%"} for \emph{cumulative percentages}. Defaults to \code{"value"}.
#' @param add_ci Logical specifying whether to include the confidence intervals. Defaults to \code{FALSE}.
#' @param cl The confidence level for the confidence intervals. Defaults to \code{0.95}.
#' @param nq Logical specifying whether to use a normal quantile (i.e. apply \code{\link[stats]{qnorm}}) in the computation of the confidence intervals.
#'   Defaults to \code{TRUE}. If \code{FALSE}, then the confidence intervals are computed directly with the \code{\link[stats]{quantile}} function.
#' @param ... Further arguments passed to or from other methods.
#' @details \code{screeplot} (of a \code{corregp} output) makes use of \code{\link[gplots]{barplot2}} from the package \pkg{gplots}.  
#' @return A plot window containing the scree plot.
#' @seealso \code{\link{corregp}}, \code{\link{summary.corregp}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' screeplot(haireye.crg, add_ci = TRUE)
#' @export
screeplot.corregp <- function(x,type="value",add_ci=FALSE,cl=0.95,nq=TRUE,...) {
	crg <- x
	if (add_ci && is.null(crg$conf)) {
		add_ci <- FALSE
		warning(paste("no bootstrapping was done in",as.list(match.call())$x,sep=" ",collapse=NULL),call.=FALSE)
		}
	type <- tolower(type)
	if (add_ci) {
		r.int <- summary(object=crg,parm=NULL,contrib=NULL,nf=NULL,add_ci=TRUE,cl=cl,nq=nq)$eigen[[type]]
		r.val <- r.int[1,]
		}
	else {
		r.val <- summary(object=crg,parm=NULL,contrib=NULL,nf=NULL,add_ci=FALSE)$eigen[type,]
		}
	gplots::barplot2(r.val,plot.ci=add_ci,ci.l=r.int[2,],ci.u=r.int[3,],...)
	invisible()
	}

#' Getting \code{conf} Components from \code{corregp} Objects
#'
#' Internal function for retrieving the \code{conf} component(s) in a \code{corregp} object.  
#' @param crg The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param parm The parameter for which to retrieve the \code{conf} components. Can be either \code{"y"}, \code{"x"}, or any vector of term names in X, level names
#'   in X or level names in Y.
#' @details
#' \code{confGet} is an internal function to be called by \code{\link{cint.corregp}}, \code{\link{cell.corregp}} or \code{\link{cell3d.corregp}}, but not by users.  
#' @return A list of components selected with \code{parm}.
#' @export
confGet <- function(crg,parm) {
	stopifnot(class(crg)=="corregp")
	c.out <- NULL
	if (length(parm)==1 && parm=="x") {
			if(is.list(crg$x)) {
				c.out <- do.call(what=c,args=crg$conf$x)
				names(c.out) <- unlist(lapply(crg$x,rownames))
				}
			else {
				c.out <- crg$conf$x
				}
		}
	if (is.list(crg$x) && all(parm %in% names(crg$x))) {
		c.out <- do.call(what=c,args=crg$conf$x[parm])
		names(c.out) <- unlist(lapply(crg$x[parm],rownames))
		}
	if (is.list(crg$x) && all(parm %in% unlist(lapply(crg$x,rownames)))) {
		c.out <- do.call(what=c,args=crg$conf$x)
		names(c.out) <- unlist(lapply(crg$x,rownames))
		c.out <- c.out[parm]
		}
	if (all(parm %in% rownames(crg$x))) {
		c.out <- crg$conf$x[parm]
		}
	if (all(parm %in% rownames(crg$y))) {
		c.out <- crg$conf$y[parm]
		}
	if (length(parm)==1 && parm=="y") {
		c.out <- crg$conf$y
		}
	return(c.out)
	}

#' Confidence Intervals for Correspondence Regression
#'
#' Method to compute confidence intervals for coordinates in correspondence regression.  
#' @param object The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param parm The parameter for which to compute the confidence intervals. Can be either \code{"y"}, \code{"x"}, or any vector of term names in X, level names in X or
#'   level names in Y. Defaults to \code{"x"}.
#' @param axis The axis for which to compute the confidence intervals.
#' @param cl The confidence level for the confidence interval. Defaults to \code{0.95}.
#' @param nq Logical specifying whether to use a normal quantile (i.e. apply \code{\link[stats]{qnorm}}) in the computation of the confidence intervals.
#' @param ... Further arguments passed to or from other methods.
#'   Defaults to \code{TRUE}. If \code{FALSE}, then the confidence intervals are computed directly with the \code{\link[stats]{quantile}} function.
#' @details
#' \code{cint} (of a \code{corregp} output) makes use of \code{\link{ci}}.  
#'
#' Typically, \code{cint} is not so much used directly as it is called by a \code{\link{ciplot.corregp}} command.  
#' @return A matrix with \code{Lower} and \code{Upper} confidence limits for the coordinates of interest.
#' @seealso \code{\link{ci}}, \code{\link{ciplot.corregp}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' cint(haireye.crg, parm = "y", axis = 1)
#' cint(haireye.crg, parm = c("Hair", "Sex"), axis = 1)
#' @export
cint.corregp <- function(object,parm="x",axis,cl=0.95,nq=TRUE,...) {
	crg <- object
	if (is.null(crg$conf)) {
		stop(paste("no bootstrapping was done in",as.list(match.call())$object,sep=" ",collapse=NULL),call.=FALSE)
		}
	if (is.character(axis)) {
		axis <- match(axis,table=colnames(crg$y))
		}
	p.con <- confGet(crg=crg,parm=parm)
	return(do.call(what=rbind,args=lapply(p.con,function(p){ci(p[,axis],cl=cl,nq=nq)})))
	}

#' @rdname cint.corregp
#' @export
cint <- function(object,...) {
	UseMethod("cint")
	}

#' Plotting Confidence Intervals for Correspondence Regression
#'
#' Method to plot confidence intervals for coordinates in correspondence regression.  
#' @param x The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param parm The parameter for which to plot the confidence intervals. Can be either \code{"y"}, \code{"x"}, or any vector of term names in X, level names in X or
#'   level names in Y. Defaults to \code{"x"}.
#' @param axis The axis for which to plot the confidence intervals.
#' @param cl The confidence level for the confidence intervals. Defaults to \code{0.95}.
#' @param nq Logical specifying whether to use a normal quantile (i.e. apply \code{\link[stats]{qnorm}}) in the computation of the confidence intervals.
#'   Defaults to \code{TRUE}. If \code{FALSE}, then the confidence intervals are computed directly with the \code{\link[stats]{quantile}} function.
#' @param horiz Logical specifying whether the confidence intervals should be plotted horizontally or not. Defaults to \code{FALSE}.
#' @param na.rm Logical specifying whether to omit \code{NA} coordinates from the plot. Defaults to \code{FALSE}.
#' @param type The type of plot: see \code{\link[graphics]{plot.default}}. For correspondence regression, there is an additional option \code{"labs"} which
#'   plots the text labels at the centers of the confidence intervals. Defaults to \code{"p"}.
#' @param col Color of the text labels: either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param cex Character expansion factor: a number to specify the size of the text labels.
#' @param font Font of the text labels: \code{1} for plain, \code{2} for bold, \code{3} for italic, and \code{4} for bold italic. Defaults to \code{1}.
#' @param alim Vector of two values specifying the lower and upper limit between which to plot the axis.
#' @param adir Reading direction of the text labels on the (horizontal) axis: either a \code{numeric} value between \code{0} and \code{3} (see the \code{las} argument in
#'   the graphical parameters \code{\link[graphics]{par}}) or a \code{character} value matching either \code{"horizontal"} or \code{"vertical"}. Defaults to \code{1} (horizontal).
#' @param ecol Color of the error bars: either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param ewid Width of the error bars: a number to specify the line width.
#' @param etyp Line type of the error bars: \code{0} or \code{"blank"}, \code{1} or \code{"solid"}, \code{2} or \code{"dashed"}, \code{3} or \code{"dotted"},
#'   \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or \code{"twodash"}. Defaults to \code{1}.
#' @param psym The symbol (or "plotting character") to use for the centers of the confidence intervals.
#' @param pcol Color of the center symbol: either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param pcex Character expansion factor of the center symbol.
#' @param pbgc Background color of the center symbol: either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param lwd Width of all lines except for the error bars, e.g. the connecting lines: a number to specify the line width.
#' @param lty Line type of all lines except for the error bars, e.g. the connecting lines: \code{0} or \code{"blank"}, \code{1} or \code{"solid"}, \code{2} or \code{"dashed"},
#'   \code{3} or \code{"dotted"}, \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or \code{"twodash"}. Defaults to \code{1}.
#' @param sfrac Width of "crossbar" at the end of error bar as a fraction of the x plotting region. Defaults to 0.01.
#' @param gap Space left between the center of the error bar and the lines marking the error bar in units of the height (width) of the letter "O". Defaults to 1.0.
#' @param main The main title of the plot.
#' @param sub The subtitle of the plot.
#' @param ... Further arguments passed to or from other methods.
#' @details \code{ciplot} (of a \code{corregp} output) makes use of \code{\link[gplots]{plotCI}} from the package \pkg{gplots}.  
#' @return A plot window containing the confidence intervals.
#' @seealso \code{\link{ci}}, \code{\link[gplots]{plotCI}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' ciplot(haireye.crg, parm = "y", axis = 1)
#' ciplot(haireye.crg, parm = c("Hair", "Sex"), axis = 1)
#' @export
ciplot.corregp <- function(x,parm="x",axis,cl=0.95,nq=TRUE,horiz=FALSE,na.rm=FALSE,type="p",col="darkgrey",cex=par("cex"),font=par("font"),alim=NULL,adir=1,ecol="darkgrey",ewid=par("lwd"),etyp=par("lty"),psym=16,pcol=par("col"),pcex=cex,pbgc=par("bg"),lwd=ewid,lty=etyp,sfrac=0.01,gap=1,main=NULL,sub=NULL,...) {
	crg <- x
	if (is.null(crg$conf)) {
		stop(paste("no bootstrapping was done in",as.list(match.call())$x,sep=" ",collapse=NULL),call.=FALSE)
		}
	if (type=="labs" && horiz) {
		horiz <- FALSE
		warning("type='labs' will plot the confidence intervals vertically",call.=FALSE)
		}
	if (is.character(axis)) {
		axis <- match(axis,table=colnames(crg$y))
		}
	if (is.character(adir)) {
		adir <- pmatch(tolower(adir),table=c("horizontal","vertical"))
		}
	a.lab <- colnames(crg$y)[axis]
	p.loc <- NULL
	if (length(parm)==1 && parm=="x") {
			if(is.list(crg$x)) {
				p.loc <- do.call(what=rbind,args=crg$x)[,axis]
				}
			else {
				p.loc <- crg$x[,axis]
				}
		}
	if (is.list(crg$x) && all(parm %in% names(crg$x))) {
		p.loc <- do.call(what=rbind,args=crg$x[parm])[,axis]
		}
	if (is.list(crg$x) && all(parm %in% unlist(lapply(crg$x,rownames)))) {
		p.loc <- do.call(what=rbind,args=crg$x)[parm,axis]
		if (length(p.loc)==1) {
			names(p.loc) <- parm
			}
		}
	if (all(parm %in% rownames(crg$x))) {
		p.loc <- crg$x[parm,axis]
		if (length(p.loc)==1) {
			names(p.loc) <- parm
			}
		}
	if (all(parm %in% rownames(crg$y))) {
		p.loc <- crg$y[parm,axis]
		if (length(p.loc)==1) {
			names(p.loc) <- parm
			}
		}
	if (length(parm)==1 && parm=="y") {
		p.loc <- crg$y[,axis]
		}
	p.int <- cint(object=crg,parm=parm,axis=axis,cl=cl,nq=nq)
	if (na.rm) {
		p.int <- p.int[!is.na(p.loc),]
		p.loc <- p.loc[!is.na(p.loc)]
		}
	p.mai <- par("mai")
	n.mai <- p.mai
	if (horiz) {
		a.err <- "x"
		x.arg <- p.loc
		y.arg <- 1:length(p.loc)
		if (is.null(alim)) {
			x.lim <- signif(range(range(p.int[,1],na.rm=TRUE),range(p.int[,2],na.rm=TRUE)))
			}
		else {
			x.lim <- alim
			}
		y.lim <- c(1,length(p.loc))
		x.axt <- "s"
		y.axt <- "n"
		x.lab <- a.lab
		y.lab <- ""
		a.num <- 2
		a.las <- 1
		n.mai[2] <- max(c(p.mai[2] - 0.3,strwidth(names(p.loc),units="inches",cex=cex,font=font))) + 0.3
		}
	else {
		a.err <- "y"
		x.arg <- 1:length(p.loc)
		y.arg <- p.loc
		x.lim <- c(1,length(p.loc))
		if (is.null(alim)) {
			y.lim <- signif(range(range(p.int[,1],na.rm=TRUE),range(p.int[,2],na.rm=TRUE)))
			}
		else {
			y.lim <- alim
			}
		x.axt <- "n"
		y.axt <- "s"
		x.lab <- ""
		y.lab <- a.lab
		a.num <- 1
		a.las <- adir
		if (adir > 1) {
			n.mai[1] <- max(c(p.mai[1] - 0.3,strwidth(names(p.loc),units="inches",cex=cex,font=font))) + 0.3
			}
		}
	par(mai=n.mai)
	gplots::plotCI(x=x.arg,y=y.arg,uiw=NULL,ui=p.int[,2],li=p.int[,1],err=a.err,xlim=x.lim,ylim=y.lim,pch=NA,barcol=ecol,lwd=ewid,lty=etyp,sfrac=sfrac,gap=gap,xaxt=x.axt,yaxt=y.axt,labels=FALSE,add=FALSE,xlab=x.lab,ylab=y.lab,main=main,sub=sub,...)
	if (type=="labs") {
		text(x=x.arg,y=y.arg,labels=names(p.loc),col=col,cex=cex,font=font)
		}
	else {
		points(x=x.arg,y=y.arg,type=type,pch=psym,col=pcol,cex=pcex,bg=pbgc,lwd=lwd,lty=lty)
		axis(side=a.num,at=1:length(p.loc),labels=names(p.loc),col.axis=col,cex.axis=cex,font.axis=font,las=a.las)
		}
	par(mai=p.mai)
	invisible()
	}

#' @rdname ciplot.corregp
#' @export
ciplot <- function(x,...) {
	UseMethod("ciplot")
	}

#' Parallel Coordinate Plotting for Correspondence Regression
#'
#' Method to produce a \emph{parallel coordinate plot} of the output of a correspondence regression.  
#' @param x The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param parm The parameter for which to plot the coordinates. Can be either \code{"y"}, \code{"x"}, or any vector of term names in X, level names in X or
#'   level names in Y. Defaults to \code{"x"}.
#' @param axes The axes to plot.
#' @param add_ci Logical specifying whether to include the confidence intervals. Defaults to \code{FALSE}.
#' @param cl The confidence level for the confidence intervals. Defaults to \code{0.95}.
#' @param nq Logical specifying whether to use a normal quantile (i.e. apply \code{\link[stats]{qnorm}}) in the computation of the confidence intervals.
#'   Defaults to \code{TRUE}. If \code{FALSE}, then the confidence intervals are computed directly with the \code{\link[stats]{quantile}} function.
#' @param col Color of the text labels: either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param cex Character expansion factor: a number to specify the size of the text labels.
#' @param font Font of the text labels: \code{1} for plain, \code{2} for bold, \code{3} for italic, and \code{4} for bold italic. Defaults to \code{1}.
#' @param lwd Width of the connecting lines: a number to specify the line width.
#' @param lty Line type of the connecting lines: \code{0} or \code{"blank"}, \code{1} or \code{"solid"}, \code{2} or \code{"dashed"}, \code{3} or
#'   \code{"dotted"}, \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or \code{"twodash"}. Defaults to \code{1}.
#' @param lcol Color of the connecting lines: either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param psym The symbol (or "plotting character") for the values of the coordinates on the axes.
#' @param pcol Color of the symbol for the values on the axes: either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param pcex Character expansion factor of the symbol for the values on the axes.
#' @param ecol Color of the error lines (connecting the confidence intervals on each axis): either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param ewid Width of the error lines (connecting the confidence intervals on each axis): a number to specify the line width.
#' @param etyp Line type of the error lines (connecting the confidence intervals on each axis): \code{0} or \code{"blank"}, \code{1} or \code{"solid"},
#'   \code{2} or \code{"dashed"}, \code{3} or \code{"dotted"}, \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or
#'   \code{"twodash"}. Defaults to \code{2}.
#' @param acol Color of the parallel axes: either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param awid Width of the parallel axes: a number to specify the line width.
#' @param atyp Line type of the parallel axes: \code{0} or \code{"blank"}, \code{1} or \code{"solid"}, \code{2} or \code{"dashed"}, \code{3} or
#'   \code{"dotted"}, \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or \code{"twodash"}. Defaults to \code{1}.
#' @param acex Character expansion factor for the labels of the parallel axes.
#' @param afnt Font for the labels of the parallel axes: \code{1} for plain, \code{2} for bold, \code{3} for italic, and \code{4} for bold italic.
#' @param adir Reading direction of the labels on the parallel axes: either a \code{numeric} value between \code{0} and \code{3} (see the \code{las}
#'   argument in the graphical parameters \code{\link[graphics]{par}}) or a \code{character} value matching either \code{"horizontal"} or
#'   \code{"vertical"}. Defaults to \code{1} (horizontal).
#' @param add_scale Logical specifying whether to add a scale for the parallel axes (which are normalised).
#' @param main The main title of the plot.
#' @param sub The subtitle of the plot.
#' @param ... Further arguments passed to or from other methods.
#' @details Although adding lines for confidence intervals is possible, it is not recommended, as it typically leads to an unreadable plot.  
#' @return A parallel coordinate plot containing the output of a correspondence regression.
#' @seealso \code{\link{ciplot.corregp}}, \code{\link{plot.corregp}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' pcplot(haireye.crg, parm = "y", axes = 1:3)
#' pcplot(haireye.crg, parm = c("Hair", "Sex"), axes = 1:3)
#' @export
pcplot.corregp <- function(x,parm="x",axes,add_ci=FALSE,cl=0.95,nq=TRUE,col="darkgrey",cex=par("cex"),font=par("font"),lwd=par("lwd"),lty=par("lty"),lcol=col,psym=NULL,pcol=col,pcex=cex,ecol="red",ewid=1,etyp=2,acol="black",awid=1,atyp=1,acex=cex,afnt=font,adir=1,add_scale=FALSE,main=NULL,sub=NULL,...) {
	crg <- x
	if (add_ci && is.null(crg$conf)) {
		add_ci <- FALSE
		warning(paste("no bootstrapping was done in",as.list(match.call())$x,sep=" ",collapse=NULL),call.=FALSE)
		}
	if (length(axes)==1) {
		stop("pcplot is suited for more than 1 axis; use 'ciplot' instead",call.=FALSE)
		}
	else {
		a.len <- length(axes)
		}
	if (is.character(axes)) {
		axes <- match(axes,table=colnames(crg$y))
		}
	if (is.character(adir)) {
		adir <- pmatch(tolower(adir),table=c("horizontal","vertical"))
		}
	if (length(parm)==1 && parm=="x") {
			if(is.list(crg$x)) {
				p.loc <- do.call(what=rbind,args=crg$x)[,axes]
				}
			else {
				p.loc <- crg$x[,axes]
				}
		}
	if (is.list(crg$x) && all(parm %in% names(crg$x))) {
		p.loc <- do.call(what=rbind,args=crg$x[parm])[,axes]
		}
	if (is.list(crg$x) && all(parm %in% unlist(lapply(crg$x,rownames)))) {
		p.loc <- do.call(what=rbind,args=crg$x)[parm,axes]
		if (nrow(p.loc)==1) {
			rownames(p.loc) <- parm
			}
		}
	if (all(parm %in% rownames(crg$x))) {
		p.loc <- rbind(crg$x[parm,axes])
		if (nrow(p.loc)==1) {
			rownames(p.loc) <- parm
			}
		}
	if (all(parm %in% rownames(crg$y))) {
		p.loc <- rbind(crg$y[parm,axes])
		if (nrow(p.loc)==1) {
			rownames(p.loc) <- parm
			}
		}
	if (length(parm)==1 && parm=="y") {
		p.loc <- crg$y[,axes]
		}
	p.min <- apply(p.loc,2,min,na.rm=TRUE)
	p.max <- apply(p.loc,2,max,na.rm=TRUE)
	if (add_ci) {
		a.int <- lapply(axes,function(a){cint(object=crg,parm=parm,axis=a,cl=cl,nq=nq)})
		p.min <- pmin(p.min,sapply(a.int,function(a1){min(a1[,1],na.rm=TRUE)}))
		p.max <- pmax(p.max,sapply(a.int,function(a2){max(a2[,2],na.rm=TRUE)}))
		}
	a.loc <- sweep(sweep(p.loc,2,p.min,"-"),2,(p.max-p.min),"/")
	p.mai <- par("mai")
	n.mai <- p.mai
	n.mai[2] <- max(strwidth(rownames(a.loc),units="inches",cex=cex,font=font)) + 0.3
	if (adir > 1) {
		n.mai[1] <- max(c(p.mai[1] - 0.3,strwidth(colnames(a.loc),units="inches",cex=acex,font=afnt))) + 0.3
		}
	par(mai=n.mai)
	matplot(t(a.loc),type="n",xlim=c(1,a.len),ylim=c(0,1),ylab="",main=main,sub=sub,axes=FALSE,add=FALSE,...)
	axis(side=1,at=1:a.len,labels=colnames(a.loc),col.axis=acol,cex.axis=acex,font.axis=afnt,las=adir)
	if (add_scale) {
		axis(side=4,col.axis=acol,cex.axis=acex,font.axis=afnt,las=1)
		}
	abline(v=1:a.len,col=acol,lty=atyp,lwd=awid)
	matlines(t(a.loc),col=lcol,lty=lty,lwd=lwd)
	if (add_ci) {
		matlines(t(sweep(sweep(sapply(a.int,function(a3){a3[,1]}),2,p.min,"-"),2,(p.max-p.min),"/")),col=ecol,lty=etyp,lwd=ewid)
		matlines(t(sweep(sweep(sapply(a.int,function(a4){a4[,2]}),2,p.min,"-"),2,(p.max-p.min),"/")),col=ecol,lty=etyp,lwd=ewid)
		}
	psym <- ifelse(is.null(psym),NA,psym)
	if (!is.na(psym)) {
		matpoints(t(a.loc),pch=psym,col=pcol,cex=pcex)
		}
	mtext(text=rownames(a.loc),side=2,at=a.loc[,1],las=1,col=col,cex=cex,font=font)
	par(mai=p.mai)
	invisible()
	}

#' @rdname pcplot.corregp
#' @export
pcplot <- function(x,...) {
	UseMethod("pcplot")
	}

#' Confidence Ellipses for Correspondence Regression
#'
#' Method to compute confidence ellipses for coordinates in correspondence regression.  
#' @param object The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param parm The parameter for which to compute the confidence ellipses. Can be either \code{"y"}, \code{"x"}, or any vector of term names in X, level names in X or
#'   level names in Y. Defaults to \code{"x"}.
#' @param axes The axes for which to compute the confidence ellipses: a vector of two values. Defaults to the first two axes.
#' @param cl The confidence level for the confidence ellipses. Defaults to \code{0.95}.
#' @param np The number of points to represent the confidence ellipses. Defaults to \code{100}.
#' @param ... Further arguments passed to or from other methods.
#' @details
#' \code{cell} (of a \code{corregp} output) makes use of \code{\link[ellipse]{ellipse}} from the package \pkg{ellipse}.  
#'
#' Typically, \code{cell} is not so much used directly as it is called by a \code{\link{plot.corregp}} command.  
#' @return A list containing \code{np} points for each confidence ellipse of interest.
#' @seealso \code{\link{plot.corregp}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' cell(haireye.crg, parm = "y")
#' cell(haireye.crg, parm = c("Hair", "Sex"))
#' @export
cell.corregp <- function(object,parm="x",axes=1:2,cl=0.95,np=100,...) {
	crg <- object
	if (is.null(crg$conf)) {
		stop(paste("no bootstrapping was done in",as.list(match.call())$object,sep=" ",collapse=NULL),call.=FALSE)
		}
	if (is.character(axes)) {
		axes <- match(axes,table=colnames(crg$y))
		}
	p.con <- confGet(crg=crg,parm=parm)
	p.con <- p.con[unlist(lapply(p.con,function(p0){any(complete.cases(p0))}))]
	return(lapply(p.con,function(p1){ellipse::ellipse(cov(p1[complete.cases(p1),axes]),centre=apply(p1[complete.cases(p1),axes],2,mean),level=cl,npoints=np)}))
	}

#' @rdname cell.corregp
#' @export
cell <- function(object,...) {
	UseMethod("cell")
	}

#' Plotting Correspondence Regression
#'
#' Basic method to plot the output of a correspondence regression.  
#' @param x The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param axes The axes to plot: a vector of two values. Defaults to the first two axes.
#' @param y_btm Logical specifying whether the Y levels should be plotted first ("at the bottom") and then be overlaid by the X levels. Defaults to \code{TRUE}.
#' @param y_ell Logical specifying whether the confidence ellipses of the Y levels should be plotted. Defaults to \code{FALSE}.
#' @param x_ell Logical specifying whether the confidence ellipses of the X levels should be plotted. Defaults to \code{FALSE}.
#' @param ysub Vector of indices to select a subset of the Y levels.
#' @param xsub Vector of indices to select a subset of the X levels.
#' @param hlim Vector of two values specifying the lower and upper limit between which to plot the horizontal axis.
#' @param vlim Vector of two values specifying the lower and upper limit between which to plot the vertical axis.
#' @param expa_btm Expansion factor for the bottom coordinates: a number to rescale the axes.
#' @param expa_top Expansion factor for the top coordinates: a number to rescale the axes.
#' @param asp The aspect ratio for the whole plot. See \code{\link[graphics]{plot.window}}.
#' @param asp_btm The aspect ratio for the bottom coordinates. See \code{\link[graphics]{plot.window}}.
#' @param asp_top The aspect ratio for the top coordinates. See \code{\link[graphics]{plot.window}}.
#' @param col_btm Color of the bottom levels: either \code{numeric} or see \code{\link[grDevices]{colors}}. Defaults to \code{"darkgrey"}.
#' @param col_top Color of the top levels: either \code{numeric} or see \code{\link[grDevices]{colors}}. Defaults to \code{"red"}.
#' @param cex_btm Character expansion factor of the bottom levels: a number to specify the size of the text labels.
#' @param cex_top Character expansion factor of the top levels: a number to specify the size of the text labels.
#' @param font_btm Font of the bottom levels: \code{1} for plain, \code{2} for bold, \code{3} for italic, and \code{4} for bold italic. Defaults to \code{1}.
#' @param font_top Font of the top levels: \code{1} for plain, \code{2} for bold, \code{3} for italic, and \code{4} for bold italic. Defaults to \code{1}.
#' @param col_ell Color of the confidence ellipses: either a number or see \code{\link[grDevices]{colors}}.
#' @param lwd_ell Width of the confidence ellipses: a number to specify the line width.
#' @param lty_ell Line type of the confidence ellipses: \code{0} or \code{"blank"}, \code{1} or \code{"solid"}, \code{2} or \code{"dashed"}, \code{3} or \code{"dotted"},
#'   \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or \code{"twodash"}. Defaults to \code{1}.
#' @param col_ori Color of the lines through the origin: either \code{numeric} or see \code{\link[grDevices]{colors}}.
#' @param lwd_ori Width of the lines through the origin: a number to specify the line width.
#' @param lty_ori Line type of the lines through the origin: \code{0} or \code{"blank"}, \code{1} or \code{"solid"}, \code{2} or \code{"dashed"},
#'   \code{3} or \code{"dotted"}, \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or \code{"twodash"}. Defaults to \code{1}.
#' @param main The main title of the plot.
#' @param sub The subtitle of the plot.
#' @param hlab The title of the horizontal axis.
#' @param vlab The title of the vertical axis.
#' @param cl The confidence level for the confidence ellipses. Defaults to \code{0.95}.
#' @param np The number of points to represent the confidence ellipses. Defaults to \code{100}.
#' @param add_ori Logical specifying whether to add lines through the origin. Defaults to \code{TRUE}.
#' @param ... Further arguments passed to or from other methods.
#' @details The plot of a correspondence regression is by definition a \code{\link[stats]{biplot}}.  
#' @return A plot window containing the output of a correspondence regression.
#' @references
#' Gower, J., S. Lubbe and N. Le Roux (2011) \emph{Understanding biplots}. Chichester: Wiley.
#'
#' Greenacre, M. (2010) \emph{Biplots in practice}. Bilbao: Fundacion BBVA.
#' @seealso \code{\link{corregp}}, \code{\link{summary.corregp}}, \code{\link{screeplot.corregp}}, \code{\link[stats]{biplot}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' plot(haireye.crg, x_ell = TRUE, xsub = c("Hair", "Sex"))
#' @export
plot.corregp <- function(x,axes=1:2,y_btm=TRUE,y_ell=FALSE,x_ell=FALSE,ysub=NULL,xsub=NULL,hlim=NULL,vlim=NULL,expa_btm=1,expa_top=1,asp=1,asp_btm=asp,asp_top=asp,col_btm="darkgrey",col_top="red",cex_btm=par("cex"),cex_top=cex_btm,font_btm=par("font"),font_top=font_btm,col_ell=par("col"),lwd_ell=par("lwd"),lty_ell=par("lty"),col_ori=par("col"),lwd_ori=par("lwd"),lty_ori=1,main=NULL,sub=NULL,hlab=NULL,vlab=NULL,cl=0.95,np=100,add_ori=TRUE,...) {
	crg <- x
	if (is.character(axes)) {
		axes <- match(axes,table=colnames(crg$y))
		}
	y.mat <- crg$y[,axes]
	if (is.list(crg$x)) {
		x.mat <- do.call(what=rbind,args=crg$x)[,axes]
		if (is.numeric(xsub) || all(xsub %in% names(crg$x))) {
			xsub <- unlist(lapply(crg$x[xsub],rownames))
			}
		}
	else {
		x.mat <- crg$x[,axes]
		}
	if (is.character(ysub)) {
		ysub <- match(ysub,table=rownames(y.mat))
		}
	if (is.null(ysub)) {
		ysub <- 1:nrow(y.mat)
		}
	if (is.character(xsub)) {
		xsub <- match(xsub,table=rownames(x.mat))
		}
	if (is.null(xsub)) {
		xsub <- 1:nrow(x.mat)
		}
	ell_btm <- NULL
	ell_top <- NULL
	if ((y_ell || x_ell) && is.null(crg$conf)) {
		y_ell <- x_ell <- FALSE
		warning(paste("no bootstrapping was done in",as.list(match.call())$x,sep=" ",collapse=NULL),call.=FALSE)
		}
	if (y_btm) {
		loc_btm <- y.mat
		loc_top <- x.mat
		sub_btm <- ysub
		sub_top <- xsub
		if (y_ell) {
			ell_btm <- cell(object=crg,parm=rownames(y.mat)[ysub],axes=axes,cl=cl,np=np)
			}
		if (x_ell) {
			ell_top <- cell(object=crg,parm=rownames(x.mat)[xsub],axes=axes,cl=cl,np=np)
			}
		}
	else {
		loc_btm <- x.mat
		loc_top <- y.mat
		sub_btm <- xsub
		sub_top <- ysub
		if (y_ell) {
			ell_top <- cell(object=crg,parm=rownames(y.mat)[ysub],axes=axes,cl=cl,np=np)
			}
		if (x_ell) {
			ell_btm <- cell(object=crg,parm=rownames(x.mat)[xsub],axes=axes,cl=cl,np=np)
			}
		}
	if (is.null(hlim)) {
		hlim <- signif(range(range(y.mat[,1],na.rm=TRUE),range(x.mat[,1],na.rm=TRUE)))
		}
	if (is.null(vlim)) {
		vlim <- signif(range(range(y.mat[,2],na.rm=TRUE),range(x.mat[,2],na.rm=TRUE)))
		}
	if (is.null(hlab)) {
		hlab <- colnames(y.mat)[1]
		}
	if (is.null(vlab)) {
		vlab <- colnames(y.mat)[2]
		}
	par(pty="s")
	plot(x=loc_btm[sub_btm,1],y=loc_btm[sub_btm,2],type="n",xlim=(hlim/expa_btm),ylim=(vlim/expa_btm),main=main,sub=sub,xlab=hlab,ylab=vlab,asp=asp_btm,...)
	if (add_ori) {
		abline(h=0,v=0,col=col_ori,lwd=lwd_ori,lty=lty_ori)
		}
	if (!is.null(ell_btm)) {
		lapply(ell_btm,function(p1){lines(x=p1[,1],y=p1[,2],col=col_ell,lwd=lwd_ell,lty=lty_ell)})
		}
	text(x=loc_btm[sub_btm,1],y=loc_btm[sub_btm,2],labels=rownames(loc_btm)[sub_btm],col=col_btm,cex=cex_btm,font=font_btm)
	par(new=TRUE)
	plot(x=loc_top[sub_top,1],y=loc_top[sub_top,2],type="n",xlim=(hlim/expa_top),ylim=(vlim/expa_top),main=NA,sub=NA,xlab=NA,ylab=NA,asp=asp_top,axes=FALSE)
	if (expa_btm != expa_top) {
		axis(side=3,col.ticks=col_top,...)
		axis(side=4,col.ticks=col_top,...)
		}
	if (!is.null(ell_top)) {
		lapply(ell_top,function(p2){lines(x=p2[,1],y=p2[,2],col=col_ell,lwd=lwd_ell,lty=lty_ell)})
		}
	text(x=loc_top[sub_top,1],y=loc_top[sub_top,2],labels=rownames(loc_top)[sub_top],col=col_top,cex=cex_top,font=font_top)
	invisible()
	}

#' 3D Confidence Ellipsoids for Correspondence Regression
#'
#' Method to compute 3D confidence ellipsoids for coordinates in correspondence regression.  
#' @param object The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param parm The parameter for which to compute the confidence ellipsoids. Can be either \code{"y"}, \code{"x"}, or any vector of term names in X, level names in X or
#'   level names in Y. Defaults to \code{"x"}.
#' @param axes The axes for which to compute the confidence ellipsoids: a vector of three values. Defaults to the first three axes.
#' @param cl The confidence level for the confidence ellipsoids. Defaults to \code{0.95}.
#' @param ... Further arguments passed to or from other methods.
#' @details
#' \code{cell3d} (of a \code{corregp} output) makes use of \code{\link[rgl]{ellipse3d}} from the package \pkg{rgl}.  
#' 
#' Typically, \code{cell3d} is not so much used directly as it is called by a \code{\link{plot3d.corregp}} command.  
#' @return A list containing coordinate points for each confidence ellipsoid of interest.
#' @seealso \code{\link{plot3d.corregp}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' cell3d(haireye.crg, parm = "y")
#' cell3d(haireye.crg, parm = c("Hair", "Sex"))
#' @export
cell3d.corregp <- function(object,parm="x",axes=1:3,cl=0.95,...) {
	crg <- object
	if (is.null(crg$conf)) {
		stop(paste("no bootstrapping was done in",as.list(match.call())$object,sep=" ",collapse=NULL),call.=FALSE)
		}
	if (is.character(axes)) {
		axes <- match(axes,table=colnames(crg$y))
		}
	p.con <- confGet(crg=crg,parm=parm)
	p.con <- p.con[unlist(lapply(p.con,function(p0){any(complete.cases(p0))}))]
	return(lapply(p.con,function(p1){rgl::ellipse3d(cov(p1[complete.cases(p1),axes]),centre=apply(p1[complete.cases(p1),axes],2,mean),level=cl)}))
	}

#' @rdname cell3d.corregp
#' @export
cell3d <- function(object,...) {
	UseMethod("cell3d")
	}

#' 3D Plotting for Correspondence Regression
#'
#' Method to produce a 3D plot for a correspondence regression.  
#' @param x The output of a call to \code{\link{corregp}} (i.e. an object of class "corregp").
#' @param axes The axes to plot: a vector of three values. Defaults to the first three axes.
#' @param y_btm Logical specifying whether the Y levels should be plotted first ("at the bottom") and then be overlaid by the X levels. Defaults to \code{TRUE}.
#' @param y_ell Logical specifying whether the confidence ellipsoids of the Y levels should be plotted. Defaults to \code{FALSE}.
#' @param x_ell Logical specifying whether the confidence ellipsoids of the X levels should be plotted. Defaults to \code{FALSE}.
#' @param ysub Vector of indices to select a subset of the Y levels.
#' @param xsub Vector of indices to select a subset of the X levels.
#' @param hlim Vector of two values specifying the lower and upper limit between which to plot the horizontal axis.
#' @param vlim Vector of two values specifying the lower and upper limit between which to plot the vertical axis.
#' @param dlim Vector of two values specifying the lower and upper limit between which to plot the "depth" axis.
#' @param asp The aspect ratio for the whole plot. See \code{\link[rgl]{aspect3d}}.
#' @param col_btm Color of the bottom levels: either \code{numeric} or see \code{\link[grDevices]{colors}}. Defaults to \code{"darkgrey"}.
#' @param col_top Color of the top levels: either \code{numeric} or see \code{\link[grDevices]{colors}}. Defaults to \code{"red"}.
#' @param cex_btm Character expansion factor of the bottom levels: a number to specify the size of the text labels.
#' @param cex_top Character expansion factor of the top levels: a number to specify the size of the text labels.
#' @param font_btm Font of the bottom levels: \code{1} for plain, \code{2} for bold, \code{3} for italic, and \code{4} for bold italic.
#' @param font_top Font of the top levels: \code{1} for plain, \code{2} for bold, \code{3} for italic, and \code{4} for bold italic.
#' @param col_ell Color of the confidence ellipsoids: either a number or see \code{\link[grDevices]{colors}}. Defaults to \code{"black"}.
#' @param lwd_ell Width of the confidence ellipsoids: a number to specify the line width.
#' @param lty_ell Line type of the confidence ellipsoids: either \code{"shade"}, \code{"wire"}, or \code{"dots"}. Defaults to \code{"shade"}.
#' @param opa_ell Opaqueness of the confidence ellipsoids: a number between \code{0} for fully transparant and \code{1} for fully opaque. Defaults to \code{0.2}.
#' @param col_ori Color of the lines through the origin: either a number or see \code{\link[grDevices]{colors}}. Defaults to \code{"grey"}.
#' @param lwd_ori Width of the lines through the origin: a number to specify the line width. Defaults to \code{1}.
#' @param main The main title of the plot.
#' @param sub The subtitle of the plot.
#' @param hlab The title of the horizontal axis.
#' @param vlab The title of the vertical axis.
#' @param dlab The title of the "depth" axis.
#' @param cl The confidence level for the confidence ellipsoids. Defaults to \code{0.95}.
#' @param add_ori Logical specifying whether to add lines through the origin. Defaults to \code{TRUE}.
#' @param ... Further arguments passed to or from other methods.
#' @details \code{plot3d} (of a \code{corregp} output) makes use of \code{\link[rgl]{plot3d}} (and \code{\link[rgl]{text3d}} and \code{\link[rgl]{abclines3d}}) from the package \pkg{rgl}.  
#' @return A 3D plot window containing the output of a correspondence regression.
#' @seealso \code{\link{corregp}}, \code{\link[rgl]{plot3d}}.
#' @examples
#' data(HairEye)
#' haireye.crg <- corregp(Eye ~ Hair * Sex, data = HairEye, b = 3000)
#' \donttest{plot3d(haireye.crg, x_ell = TRUE, xsub = c("Hair", "Sex"))}
#' @export
plot3d.corregp <- function(x,axes=1:3,y_btm=TRUE,y_ell=FALSE,x_ell=FALSE,ysub=NULL,xsub=NULL,hlim=NULL,vlim=NULL,dlim=NULL,asp=par3d("scale"),col_btm="darkgrey",col_top="red",cex_btm=par3d("cex"),cex_top=cex_btm,font_btm=par3d("font"),font_top=font_btm,col_ell="black",lwd_ell=1,lty_ell="shade",opa_ell=0.2,col_ori="grey",lwd_ori=1,main=NULL,sub=NULL,hlab=NULL,vlab=NULL,dlab=NULL,cl=0.95,add_ori=TRUE,...) {
	crg <- x
	if (is.character(axes)) {
		axes <- match(axes,table=colnames(crg$y))
		}
	axes <- c(axes[1],axes[3],axes[2])
	y.mat <- crg$y[,axes]
	if (is.list(crg$x)) {
		x.mat <- do.call(what=rbind,args=crg$x)[,axes]
		if (is.numeric(xsub) || all(xsub %in% names(crg$x))) {
			xsub <- unlist(lapply(crg$x[xsub],rownames))
			}
		}
	else {
		x.mat <- crg$x[,axes]
		}
	if (is.character(ysub)) {
		ysub <- match(ysub,table=rownames(y.mat))
		}
	if (is.null(ysub)) {
		ysub <- 1:nrow(y.mat)
		}
	if (is.character(xsub)) {
		xsub <- match(xsub,table=rownames(x.mat))
		}
	if (is.null(xsub)) {
		xsub <- 1:nrow(x.mat)
		}
	ell_btm <- NULL
	ell_top <- NULL
	if ((y_ell || x_ell) && is.null(crg$conf)) {
		y_ell <- x_ell <- FALSE
		warning(paste("no bootstrapping was done in",as.list(match.call())$x,sep=" ",collapse=NULL),call.=FALSE)
		}
	if (y_btm) {
		loc_btm <- y.mat
		loc_top <- x.mat
		sub_btm <- ysub
		sub_top <- xsub
		if (y_ell) {
			ell_btm <- cell3d(object=crg,parm=rownames(y.mat)[ysub],axes=axes,cl=cl)
			}
		if (x_ell) {
			ell_top <- cell3d(object=crg,parm=rownames(x.mat)[xsub],axes=axes,cl=cl)
			}
		}
	else {
		loc_btm <- x.mat
		loc_top <- y.mat
		sub_btm <- xsub
		sub_top <- ysub
		if (y_ell) {
			ell_top <- cell3d(object=crg,parm=rownames(y.mat)[ysub],axes=axes,cl=cl)
			}
		if (x_ell) {
			ell_btm <- cell3d(object=crg,parm=rownames(x.mat)[xsub],axes=axes,cl=cl)
			}
		}
	if (is.null(hlim)) {
		hlim <- signif(range(range(y.mat[,1],na.rm=TRUE),range(x.mat[,1],na.rm=TRUE)))
		}
	if (is.null(vlim)) {
		vlim <- signif(range(range(y.mat[,2],na.rm=TRUE),range(x.mat[,2],na.rm=TRUE)))
		}
	if (is.null(dlim)) {
		dlim <- signif(range(range(y.mat[,3],na.rm=TRUE),range(x.mat[,3],na.rm=TRUE)))
		}
	if (is.null(hlab)) {
		hlab <- colnames(y.mat)[1]
		}
	if (is.null(vlab)) {
		vlab <- colnames(y.mat)[2]
		}
	if (is.null(dlab)) {
		dlab <- colnames(y.mat)[3]
		}
	rgl::plot3d(x=loc_btm[sub_btm,1],y=loc_btm[sub_btm,2],z=loc_btm[sub_btm,3],type="n",xlim=hlim,ylim=vlim,zlim=dlim,main=main,sub=sub,xlab=hlab,ylab=vlab,zlab=dlab,aspect=asp,...)
	if (add_ori) {
		rgl::abclines3d(x=0,y=0,z=0,a=1,b=0,c=0,col=col_ori,lwd=lwd_ori)
		rgl::abclines3d(x=0,y=0,z=0,a=0,b=1,c=0,col=col_ori,lwd=lwd_ori)
		rgl::abclines3d(x=0,y=0,z=0,a=0,b=0,c=1,col=col_ori,lwd=lwd_ori)
		}
	if (!is.null(ell_btm)) {
		lapply(ell_btm,rgl::plot3d,add=TRUE,xlim=hlim,ylim=vlim,zlim=dlim,lit=FALSE,col=col_ell,lwd=lwd_ell,type=lty_ell,alpha=opa_ell)
		}
	rgl::text3d(x=loc_btm[sub_btm,1],y=loc_btm[sub_btm,2],z=loc_btm[sub_btm,3],texts=rownames(loc_btm)[sub_btm],col=col_btm,cex=cex_btm,font=font_btm)
	rgl::plot3d(x=loc_top[sub_top,1],y=loc_top[sub_top,2],z=loc_top[sub_top,3],type="n",xlim=hlim,ylim=vlim,zlim=dlim,main="",sub="",xlab="",ylab="",zlab="",add=TRUE)
	if (!is.null(ell_top)) {
		lapply(ell_top,rgl::plot3d,add=TRUE,xlim=hlim,ylim=vlim,zlim=dlim,lit=FALSE,col=col_ell,lwd=lwd_ell,type=lty_ell,alpha=opa_ell)
		}
	rgl::text3d(x=loc_top[sub_top,1],y=loc_top[sub_top,2],z=loc_top[sub_top,3],texts=rownames(loc_top)[sub_top],col=col_top,cex=cex_top,font=font_top)
	invisible()
	}
