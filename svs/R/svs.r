#' @import graphics
NULL

#' @import stats
NULL

#' @import utils
NULL

#' @importFrom gtools rdirichlet
NULL

#' Tools for Semantic Vector Spaces
#'
#' This package offers various tools for semantic vector spaces. There are techniques for correspondence analysis (simple,
#'   multiple and discriminant), latent semantic analysis, probabilistic latent semantic analysis, non-negative matrix
#'   factorization, latent class analysis and EM clustering. Furthermore, the package has specialized distance measures and
#'   plotting functions as well as some helper functions.  
#' @section Contents:
#' This package contains the following raw data files (in the folder \emph{extdata}):
#' \itemize{
#'    \item{\code{\link{SndT_Fra.txt}} }{Seventeen Dutch source words and their French translations.}
#'    \item{\code{\link{SndT_Eng.txt}} }{Seventeen Dutch source words and their English translations.}
#'    \item{\code{\link{InvT_Fra.txt}} }{Seventeen Dutch target words and their French source words.}
#'    \item{\code{\link{InvT_Eng.txt}} }{Seventeen Dutch target words and their English source words.}
#'    \item{\code{\link{Ctxt_Dut.txt}} }{Context words for seventeen Dutch words.}
#'    \item{\code{\link{Ctxt_Fra.txt}} }{Context words for seventeen Dutch words translated from French.}
#'    \item{\code{\link{Ctxt_Eng.txt}} }{Context words for seventeen Dutch words translated from English.}
#' }
#' The (fast procedures for the) techniques in this package are:
#' \itemize{
#'   \item{\code{\link{fast_sca}} }{Simple correspondence analysis.}
#'   \item{\code{\link{fast_mca}} }{Multiple correspondence analysis.}
#'   \item{\code{\link{fast_dca}} }{Discriminant correspondence analysis.}
#'   \item{\code{\link{fast_lsa}} }{Latent semantic analysis.}
#'   \item{\code{\link{fast_psa}} }{Probabilistic latent semantic analysis.}
#'   \item{\code{\link{fast_nmf}} }{Non-negative matrix factorization.}
#'   \item{\code{\link{fast_lca}} }{Latent class analysis.}
#'   \item{\code{\link{fast_E_M}} }{EM clustering.}
#' }
#'
#' The complete overview of local and global weighting functions in this package can be found on \code{\link{weighting_functions}}.
#'
#' The specialized distance measures are:
#' \itemize{
#'   \item{\code{\link{dist_chisquare}} }{Chi-square distance.}
#'   \item{\code{\link{dist_cosine}} }{Cosine distance.}
#'   \item{\code{\link{dist_wrt}} }{Distance with respect to a certain point.}
#'   \item{\code{\link{dist_wrt_centers}} }{Distance with respect to cluster centers.}
#' }
#'
#' The specialized plotting functions are:
#' \itemize{
#'   \item{\code{\link{cd_plot}} }{Cumulative distribution plot.}
#'   \item{\code{\link{pc_plot}} }{Parallel coordinate plot.}
#' }
#'
#' There are two helper functions for correspondence analysis:
#' \itemize{
#'   \item{\code{\link{freq_ca}} }{Compute level frequencies (for a factor).}
#'   \item{\code{\link{centers_ca}} }{Compute coordinates for cluster centers.}
#' }
#'
#' There is one helper function for \pkg{pvclust}:
#' \itemize{
#'   \item{\code{\link{complete_pvpick}} }{Complete the output of \code{pvpick}.}
#' }
#'
#' The remaining helper functions in this package are:
#' \itemize{
#'   \item{\code{\link{vec2ind}} }{Transform a vector into an indicator matrix.}
#'   \item{\code{\link{tab2dat}} }{Transform a table into a data frame.}
#'   \item{\code{\link{tab2ind}} }{Transform a table into an indicator matrix.}
#'   \item{\code{\link{outerec}} }{Recursive application of the outer product.}
#'   \item{\code{\link{pmi}} }{Pointwise mutual information.}
#'   \item{\code{\link{MI}} }{Mutual information.}
#' }
#' @section Further reference:
#' \itemize{
#'   \item{Many packages contain correspondence analysis: \pkg{ca}, \pkg{FactoMineR}, \pkg{MASS} and others.}
#'   \item{For latent semantic analysis there is also the package \pkg{lsa}.}
#'   \item{The package \pkg{NMF} provides more flexibility for non-negative matrix factorization.}
#'   \item{For topic models there are the packages \pkg{lda} and \pkg{topicmodels}.}
#'   \item{Latent class analysis can also be run in the package \pkg{poLCA}.}
#' }
#' @section Author:
#' Koen Plevoets, \email{koen.plevoets@@ugent.be}  
#' @section Acknowledgements:
#' This package has benefited greatly from the helpful comments of Lore Vandevoorde, Pauline De Baets and Gert De Sutter. Thanks to
#'   Kurt Hornik, Uwe Ligges and Brian Ripley for their valuable recommendations when proofing this package.  
#' @docType package
#' @name svs-package
NULL

#' Seventeen Dutch Source Words and their French Translations
#'
#' The occurrences of seventeen Dutch synonyms of \emph{beginnen} ("to begin") and their French translations (from the Dutch Parallel Corpus).  
#' @format A data frame with 1487 rows and 2 variables.
#' \itemize{
#'   \item{\code{source_Dut} }{The Dutch source word.}
#'   \item{\code{target_Fra} }{The French target word.}
#' }
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt",package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.SndT_Fra <- fast_sca(SndT_Fra)
#' sca.SndT_Fra
#' lsa.SndT_Fra <- fast_lsa(SndT_Fra)
#' lsa.SndT_Fra
#' @docType data
#' @name SndT_Fra.txt
NULL

#' Seventeen Dutch Source Words and their English Translations
#'
#' The occurrences of seventeen Dutch synonyms of \emph{beginnen} ("to begin") and their English translations (from the Dutch Parallel Corpus).  
#' @format A data frame with 1117 rows and 2 variables.
#' \itemize{
#'   \item{\code{source_Dut} }{The Dutch source word.}
#'   \item{\code{target_Eng} }{The English target word.}
#' }
#' @examples
#' SndT_Eng <- read.table(system.file("extdata", "SndT_Eng.txt",package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.SndT_Eng <- fast_sca(SndT_Eng)
#' sca.SndT_Eng
#' lsa.SndT_Eng <- fast_lsa(SndT_Eng)
#' lsa.SndT_Eng
#' @docType data
#' @name SndT_Eng.txt
NULL

#' Seventeen Dutch Target Words and their French Source Words
#'
#' The occurrences of seventeen Dutch synonyms of \emph{beginnen} ("to begin") and their French source words (from the Dutch Parallel Corpus).  
#' @format A data frame with 856 rows and 2 variables.
#' \itemize{
#'   \item{\code{source_Fra} }{The French source word.}
#'   \item{\code{target_Dut} }{The Dutch target word.}
#' }
#' @examples
#' InvT_Fra <- read.table(system.file("extdata", "InvT_Fra.txt",package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.InvT_Fra <- fast_sca(InvT_Fra)
#' sca.InvT_Fra
#' lsa.InvT_Fra <- fast_lsa(InvT_Fra)
#' lsa.InvT_Fra
#' @docType data
#' @name InvT_Fra.txt
NULL

#' Seventeen Dutch Target Words and their English Source Words
#'
#' The occurrences of seventeen Dutch synonyms of \emph{beginnen} ("to begin") and their English source words (from the Dutch Parallel Corpus).  
#' @format A data frame with 782 rows and 2 variables.
#' \itemize{
#'   \item{\code{source_Eng} }{The English source word.}
#'   \item{\code{target_Dut} }{The Dutch target word.}
#' }
#' @examples
#' InvT_Eng <- read.table(system.file("extdata", "InvT_Eng.txt",package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.InvT_Eng <- fast_sca(InvT_Eng)
#' sca.InvT_Eng
#' lsa.InvT_Eng <- fast_lsa(InvT_Eng)
#' lsa.InvT_Eng
#' @docType data
#' @name InvT_Eng.txt
NULL

#' Context Words for seventeen Dutch Words
#'
#' The frequency table of seventeen Dutch synonyms of \emph{beginnen} ("to begin") and their context words (from the Dutch Parallel Corpus).  
#' @format A table with 17 rows and 1404 columns.
#' @examples
#' Ctxt_Dut <- read.table(system.file("extdata", "Ctxt_Dut.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.Ctxt_Dut <- fast_sca(data.matrix(Ctxt_Dut))
#' sca.Ctxt_Dut
#' lsa.Ctxt_Dut <- fast_lsa(data.matrix(Ctxt_Dut))
#' lsa.Ctxt_Dut
#' @docType data
#' @name Ctxt_Dut.txt
NULL

#' Context Words for seventeen Dutch Words Translated from French
#'
#' The frequency table of seventeen Dutch synonyms of \emph{beginnen} ("to begin") and their context words in texts translated from French
#'   (from the Dutch Parallel Corpus).  
#' @format A table with 17 rows and 612 columns.
#' @examples
#' Ctxt_Fra <- read.table(system.file("extdata", "Ctxt_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.Ctxt_Fra <- fast_sca(data.matrix(Ctxt_Fra))
#' sca.Ctxt_Fra
#' lsa.Ctxt_Fra <- fast_lsa(data.matrix(Ctxt_Fra))
#' lsa.Ctxt_Fra
#' @docType data
#' @name Ctxt_Fra.txt
NULL

#' Context Words for seventeen Dutch Words Translated from French
#'
#' The frequency table of seventeen Dutch synonyms of \emph{beginnen} ("to begin") and their context words in texts translated from English
#'   (from the Dutch Parallel Corpus).  
#' @format A table with 17 rows and 609 columns.
#' @examples
#' Ctxt_Eng <- read.table(system.file("extdata", "Ctxt_Eng.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.Ctxt_Eng <- fast_sca(data.matrix(Ctxt_Eng))
#' sca.Ctxt_Eng
#' lsa.Ctxt_Eng <- fast_lsa(data.matrix(Ctxt_Eng))
#' lsa.Ctxt_Eng
#' @docType data
#' @name Ctxt_Eng.txt
NULL

#' Transform a Vector into an Indicator Matrix
#'
#' A helper function for transforming a vector into an indicator matrix.  
#' @param x A vector (which will be converted to a factor).
#' @param add_names Logical specifying whether to add dimnames to the resulting indicator matrix.
#' @details
#' This is essentially the function \code{class.ind} from the package \pkg{MASS}.
#' @return An indicator matrix.
#' @export
vec2ind <- function(x,add_names=TRUE) {
	fac <- as.factor(x)
	lev <- levels(fac)
	len <- length(fac)
	out <- matrix(0,nrow=len,ncol=length(lev))
	out[(1:len)+len*(unclass(fac)-1)] <- 1
	if (add_names) {
		dimnames(out) <- list(names(fac),lev)
		}
	out
	}

#' Transform a Table into a Data Frame
#'
#' A helper function for transforming a table into a data frame.  
#' @param tab A table (i.e. an object which can be converted to an array).
#' @return A data frame.
#' @export
tab2dat <- function(tab) {
	tab <- as.array(tab)
	data.frame(lapply(expand.grid(dimnames(tab)),FUN=rep,times=as.vector(tab)))
	}

#' Transform a Table into an Indicator Matrix
#'
#' A helper function for transforming a table into an indicator matrix.  
#' @param tab A table (i.e. an object which can be converted to an array).
#' @param sep Character specifying the separator string for joining the levels.
#' @param add_names Logical specifying whether to add dimnames to the resulting indicator matrix.
#' @return An indicator matrix.
#' @export
tab2ind <- function(tab,sep="_",add_names=TRUE) {
	vec2ind(interaction(tab2dat(tab),drop=FALSE,sep=sep,lex.order=FALSE),add_names=add_names)
	}

#' Recursive Application of the Outer Product
#'
#' A helper function for computing the outer product of two or more arrays.  
#' @param ... The specification of two or more arrays (separated by comma's or contained in a list).
#' @return An array with the outer product of all the arrays specified in \code{...}.
#' @export
outerec <- function(...) {
	arg <- list(...)
	if (length(arg)==1 && is.list(arg[[1]])) {
		arg <- arg[[1]]
		}
	out <- outer(X=arg[[1]],Y=arg[[2]],FUN="*")
	if (length(arg)>2) {
		for (j in 3:length(arg)) {
			out <- outer(X=out,Y=arg[[j]],FUN="*")
			}
		}
	out
	}

#' Simple Correspondence Analysis
#'
#' A fast procedure for computing simple correspondence analysis.  
#' @param dat Input data: can be a table or a data frame (but the data frame must have only two columns).
#' @return A list with components:
#' \item{\code{val} }{The eigenvalues or principal inertias, indicating how much each latent axis explains.}
#' \item{\code{pos1} }{The coordinates of the first set of levels (\emph{viz.} the row levels of a frequency table).}
#' \item{\code{pos2} }{The coordinates of the second set of levels (\emph{viz.} the column levels of a frequency table).}
#' @references
#' Greenacre, M. (2007) \emph{Correspondence analysis in practice, Second edition}. Boca Raton: Chapman and Hall/CRC.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.SndT_Fra <- fast_sca(SndT_Fra)
#' sca.SndT_Fra
#' @export
fast_sca <- function(dat) {
	if (is.data.frame(dat)) {
		dat <- table(dat)
		}
	stopifnot(is.matrix(dat))
	exf <- outer(apply(dat,1,sum),apply(dat,2,sum))/sum(dat)
	dec <- svd((dat-exf)/sqrt(exf))
	rwc <- sweep(dec$u%*%diag(dec$d),1,sqrt(apply(dat,1,sum)),"/")
	dimnames(rwc) <- list(rownames(dat),paste("Dim",1:length(dec$d),sep=""))
	clc <- sweep(dec$v%*%diag(dec$d),1,sqrt(apply(dat,2,sum)),"/")
	dimnames(clc) <- list(colnames(dat),paste("Dim",1:length(dec$d),sep=""))
	list(val=(dec$d^2)/sum(dat),pos1=rwc,pos2=clc)
	}

#' Multiple Correspondence Analysis
#'
#' A fast procedure for computing multiple correspondence analysis.  
#' @param dat Input data: has to be a data frame (with any number of columns).
#' @param nfac Logical indicating whether the number of factors (i.e. the number of columns in \code{dat}) is a divisor for
#'   the eigenvalues (principal inertias) and the coordinates.
#' @return A list with components:
#' \item{\code{val} }{The eigenvalues or principal inertias, indicating how much each latent axis explains.}
#' \item{\code{pos} }{The coordinates of all levels.}
#' @references
#' Greenacre, M. (2007) \emph{Correspondence analysis in practice, Second edition}. Boca Raton: Chapman and Hall/CRC.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' mca.SndT_Fra <- fast_mca(SndT_Fra)
#' mca.SndT_Fra
#' @export
fast_mca <- function(dat,nfac=FALSE) {
	stopifnot(is.data.frame(dat))
	ind <- do.call(what=cbind,args=lapply(dat,FUN=vec2ind,add_names=TRUE))
	exf <- outer(apply(ind,2,sum),apply(ind,2,sum))/nrow(ind)
	dec <- eigen(((t(ind)%*%ind)-exf)/sqrt(exf))
	pos <- sweep(dec$vectors%*%diag(dec$values),1,sqrt(apply(ind,2,sum)*ifelse(nfac,ncol(dat),1)),"/")
	dimnames(pos) <- list(colnames(ind),paste("Dim",1:length(dec$values),sep=""))
	list(val=(dec$values^2)/(nrow(ind)*ifelse(nfac,ncol(dat)^2,1)),pos=pos)
	}

#' Discriminant Correspondence Analysis
#'
#' A fast procedure for computing discriminant correspondence analysis.  
#' @param dat Input data: can be a table or a data frame (but the data frame must have only two columns).
#' @param clusters1 A clustering of the first set of levels: either a list or the output of \code{kmeans}.
#' @param clusters2 A clustering of the second set of levels: either a list or the output of \code{kmeans}.
#' @param members Logical indicating whether the (supplementary) coordinates for the individual levels should also be computed.
#' @return A list with components:
#' \item{\code{val} }{The eigenvalues or principal inertias, indicating how much each latent axis explains.}
#' \item{\code{cen1} }{The coordinates of the cluster centers for the first set of levels.}
#' \item{\code{cen2} }{The coordinates of the cluster centers for the second set of levels.}
#' \item{\code{mem1} }{If \code{members = TRUE}: The coordinates of the first set of individual levels.}
#' \item{\code{mem2} }{If \code{members = TRUE}: The coordinates of the second set of individual levels.}
#' @references
#' Abdi, H. (2007) Discriminant correspondence analysis. In: N. Salkind (ed.) \emph{Encyclopedia of measurement and statistics}. Thousand Oaks: SAGE.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.SndT_Fra <- fast_sca(SndT_Fra)
#' kcl.SndT_Fra <- kmeans(sca.SndT_Fra$pos1, centers = 7)
#' dca.SndT_Fra <- fast_dca(SndT_Fra, clusters1 = kcl.SndT_Fra)
#' dca.SndT_Fra
#' @export
fast_dca <- function(dat,clusters1=NULL,clusters2=NULL,members=FALSE) {
	if (is.data.frame(dat)) {
		dat <- table(dat)
		}
	stopifnot(is.matrix(dat))
	if (class(clusters1)=="kmeans") {
		clusters1 <- lapply(sort(unique(clusters1$cluster)),function(x1){names(clusters1$cluster[clusters1$cluster==x1])})
		}
	if (class(clusters2)=="kmeans") {
		clusters2 <- lapply(sort(unique(clusters2$cluster)),function(x2){names(clusters2$cluster[clusters2$cluster==x2])})
		}
	if (is.null(clusters1)) {
		cl1 <- as.list(rownames(dat))
		dn1 <- rownames(dat)
		}
	else {
		cl1 <- clusters1
		dn1 <- 1:length(clusters1)
		}
	if (is.null(clusters2)) {
		cl2 <- as.list(colnames(dat))
		dn2 <- colnames(dat)
		}
	else {
		cl2 <- clusters2
		dn2 <- 1:length(clusters2)
		}
	grp <- do.call(what=rbind,args=lapply(cl1,function(y1){apply(rbind(dat[y1,]),2,sum)}))
	grp <- do.call(what=cbind,args=lapply(cl2,function(y2){apply(cbind(grp[,y2]),1,sum)}))
	exf <- outer(apply(grp,1,sum),apply(grp,2,sum))/sum(grp)
	dec <- svd((grp-exf)/sqrt(exf))
	rwc <- sweep(dec$u%*%diag(dec$d),1,sqrt(apply(grp,1,sum)),"/")
	dimnames(rwc) <- list(dn1,paste("Dim",1:length(dec$d),sep=""))
	clc <- sweep(dec$v%*%diag(dec$d),1,sqrt(apply(grp,2,sum)),"/")
	dimnames(clc) <- list(dn2,paste("Dim",1:length(dec$d),sep=""))
	out <- list(val=(dec$d^2)/sum(grp),cen1=rwc,cen2=clc)
	if (!is.null(clusters1) && members) {
		tmp <- do.call(what=cbind,args=lapply(cl2,function(z1){apply(cbind(dat[,z1]),1,sum)}))
		txf <- outer(apply(tmp,1,sum),apply(tmp,2,sum))/sum(tmp)
		rws <- sweep(((tmp-txf)/sqrt(txf))%*%dec$v,1,sqrt(apply(tmp,1,sum)),"/")
		dimnames(rws) <- list(rownames(dat),paste("Dim",1:length(dec$d),sep=""))
		out <- c(out,list(mem1=rws))
		}
	if (!is.null(clusters2) && members) {
		tmp <- do.call(what=rbind,args=lapply(cl1,function(z2){apply(rbind(dat[z2,]),2,sum)}))
		txf <- outer(apply(tmp,1,sum),apply(tmp,2,sum))/sum(tmp)
		cls <- sweep(t((tmp-txf)/sqrt(txf))%*%dec$u,1,sqrt(apply(tmp,2,sum)),"/")
		dimnames(cls) <- list(colnames(dat),paste("Dim",1:length(dec$d),sep=""))
		out <- c(out,list(mem2=cls))
		}
	out
	}

#' Latent Semantic Analysis
#'
#' A fast procedure for computing latent semantic analysis.  
#' @param dat Input data: can be a table or a data frame (but the data frame must have only two columns).
#' @param local_weights Character specification of the local weighting function (without a prefix): see
#'   \link[=weighting_functions]{Weighting functions}.
#' @param global_weights Character specification of the global weighting function (without a prefix): see
#'   \link[=weighting_functions]{Weighting functions}.
#' @return A list with components:
#' \item{\code{val} }{The singular values, indicating how much each latent axis explains.}
#' \item{\code{pos1} }{The coordinates of the first set of levels (\emph{viz.} the row levels of a frequency table).}
#' \item{\code{pos2} }{The coordinates of the second set of levels (\emph{viz.} the column levels of a frequency table).}
#' @references
#' Deerwester, S., S. T. Dumais, G. W. Furnas, Th. K. Landauer and R. Harshman (1990) Indexing by latent semantic analysis.
#'   \emph{Journal of the American society for information science} \strong{41} (6), 391--407.
#'
#' Landauer, Th. K. and S. T. Dumais (1997) A solution to Plato's problem: the latent semantic analysis theory of the
#'   acquisition, induction, and representation of knowledge. \emph{Psychological review} \strong{104}, 211--240.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' lsa.SndT_Fra <- fast_lsa(SndT_Fra)
#' lsa.SndT_Fra
#' @export
fast_lsa <- function(dat,local_weights="log",global_weights="idf") {
	if (is.data.frame(dat)) {
		dat <- table(dat)
		}
	stopifnot(is.matrix(dat))
	lwf <- match.fun(paste("lw",tolower(local_weights),sep="_",collapse=NULL))
	gwf <- match.fun(paste("gw",tolower(global_weights),sep="_",collapse=NULL))
	out <- svd(sweep(lwf(dat),1,gwf(dat),"*"))
	dimnames(out$u) <- list(rownames(dat),paste("Dim",1:length(out$d),sep="",collapse=NULL))
	dimnames(out$v) <- list(colnames(dat),paste("Dim",1:length(out$d),sep="",collapse=NULL))
	names(out) <- c("val","pos1","pos2")
	out
	}

#' @rdname fast_lsa
#' @export
fast_lsi <- function(dat,local_weights="log",global_weights="idf") {
	fast_lsa(dat=dat,local_weights=local_weights,global_weights=global_weights)
	}

#' Non-negative Matrix Factorization
#'
#' A fast procedure for non-negative matrix factorization.  
#' @param dat Input data: can be a table or a data frame (but the data frame must have only two columns).
#' @param k Numeric specification of the number of latent axes to compute.
#' @param type Character specification of the type of optimization: can in the current implementation be either \code{"KL"} for
#'   the Kullback-Leibler divergence, \code{"Frobenius"} or \code{"euclidean"} (or abbreviations thereof) for the euclidean distance, or
#'   \code{"ALS"} for alternating least squares.
#' @param tol Numeric specification of the convergence criterion.
#' @return A list with components:
#' \item{\code{pos1} }{The coordinates of the first set of levels (\emph{viz.} the row levels of a frequency table).}
#' \item{\code{pos2} }{The coordinates of the second set of levels (\emph{viz.} the column levels of a frequency table).}
#' @references
#' Lee, D. D. and H. S. Seung (1999) Learning the parts of objects by non-negative matrix factorization. \emph{Nature} \strong{401}, 788--791.
#'
#' Lee, D. D. and H. S. Seung (2001) Algorithms for non-negative matrix factorization. \emph{Advances in neural information processing systems}
#'   \strong{13}, 556--562.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' nmf.SndT_Fra <- fast_nmf(SndT_Fra, k = 7)
#' nmf.SndT_Fra
#' @export
fast_nmf <- function(dat,k,type="KL",tol=1e-8) {
	typ <- pmatch(x=toupper(type),table=c("KL","FROBENIUS","EUCLIDEAN","ALS"))
	if (typ==1) {
		out <- fast_nmf_KL(dat=dat,k=k,tol=tol)
		}
	if (typ==2||typ==3) {
		out <- fast_nmf_Fr(dat=dat,k=k,tol=tol)
		}
	if (typ==4) {
		out <- fast_nmf_Al(dat=dat,k=k,tol=tol)
		}
	out
	}

#' @rdname fast_nmf
#' @export
fast_nmf_KL <- function(dat,k,tol=1e-8) {
	if (is.data.frame(dat)) {
		dat <- table(dat)
		}
	stopifnot(is.matrix(dat))
	wmt <- sqrt(sum(dat)/k)*rmultinom(n=k,size=10^(ceiling(log10(nrow(dat)+1))+3),prob=apply(dat,1,sum)/sum(dat))/(10^(ceiling(log10(nrow(dat)+1))+3))
	hmt <- sqrt(sum(dat)/k)*t(rmultinom(n=k,size=10^(ceiling(log10(ncol(dat)+1))+3),prob=apply(dat,2,sum)/sum(dat))/(10^(ceiling(log10(ncol(dat)+1))+3)))
	exf <- wmt%*%hmt
	tmp <- dat*log(dat/exf)-dat+exf
	tmp[!is.finite(tmp)] <- 0
	cst <- sum(tmp)
	div <- tol*10
	itr <- FALSE
	while (div>tol) {
		wup <- prop.table(t(t(wmt)*(diag(1/apply(hmt,1,sum))%*%hmt%*%t(dat/exf))),2)
		hup <- hmt*(diag(1/apply(wmt,2,sum))%*%t(wmt)%*%(dat/exf))
		wmt <- wup
		hmt <- hup
		exf <- wmt%*%hmt
		tmp <- dat*log(dat/exf)-dat+exf
		tmp[!is.finite(tmp)] <- 0
		csn <- sum(tmp)
		div <- cst-ifelse(itr,csn,cst-tol*(10^ceiling(log10(prod(dim(dat))+1))))
		if (div>0) {
			wop <- wmt
			hop <- hmt
			}
		cst <- csn
		itr <- TRUE
		}
	wop <- wop[,order(apply(hop,1,sum),decreasing=TRUE)]
	hop <- hop[order(apply(hop,1,sum),decreasing=TRUE),]
	dimnames(wop) <- list(rownames(dat),1:k)
	dimnames(hop) <- list(1:k,colnames(dat))
	list(pos1=wop,pos2=hop)
	}

#' @rdname fast_nmf
#' @export
fast_nmf_Fr <- function(dat,k,tol=1e-8) {
	if (is.data.frame(dat)) {
		dat <- table(dat)
		}
	stopifnot(is.matrix(dat))
	wmt <- sqrt(sum(dat)/k)*rmultinom(n=k,size=10^(ceiling(log10(nrow(dat)+1))+3),prob=apply(dat,1,sum)/sum(dat))/(10^(ceiling(log10(nrow(dat)+1))+3))
	hmt <- sqrt(sum(dat)/k)*t(rmultinom(n=k,size=10^(ceiling(log10(ncol(dat)+1))+3),prob=apply(dat,2,sum)/sum(dat))/(10^(ceiling(log10(ncol(dat)+1))+3)))
	exf <- wmt%*%hmt
	cst <- sum((dat-exf)^2)/2
	div <- tol*10
	itr <- FALSE
	while (div>tol) {
		wup <- prop.table(wmt*(dat%*%t(hmt))/(wmt%*%hmt%*%t(hmt)),2)
		hup <- hmt*(t(wmt)%*%dat)/(t(wmt)%*%wmt%*%hmt)
		wmt <- wup
		hmt <- hup
		exf <- wmt%*%hmt
		csn <- sum((dat-exf)^2)/2
		div <- cst-ifelse(itr,csn,cst-tol*(10^ceiling(log10(prod(dim(dat))+1))))
		if (div>0) {
			wop <- wmt
			hop <- hmt
			}
		cst <- csn
		itr <- TRUE
		}
	wop <- wop[,order(apply(hop,1,sum),decreasing=TRUE)]
	hop <- hop[order(apply(hop,1,sum),decreasing=TRUE),]
	dimnames(wop) <- list(rownames(dat),1:k)
	dimnames(hop) <- list(1:k,colnames(dat))
	list(pos1=wop,pos2=hop)
	}

#' @rdname fast_nmf
#' @export
fast_nmf_Al <- function(dat,k,tol=1e-8) {
	if (is.data.frame(dat)) {
		dat <- table(dat)
		}
	stopifnot(is.matrix(dat))
	wmt <- sqrt(sum(dat)/k)*rmultinom(n=k,size=10^(ceiling(log10(nrow(dat)+1))+3),prob=apply(dat,1,sum)/sum(dat))/(10^(ceiling(log10(nrow(dat)+1))+3))
	hmt <- solve(t(wmt)%*%wmt,t(wmt)%*%dat)
	hmt[hmt<0] <- 0
	exf <- wmt%*%hmt
	cst <- sum((dat-exf)^2)/2
	div <- tol*10
	itr <- FALSE
	while (div>tol) {
		wup <- t(solve(hmt%*%t(hmt),hmt%*%t(dat)))
		wup[wup<0] <- 0
		wup <- prop.table(wup,2)
		hup <- solve(t(wmt)%*%wmt,t(wmt)%*%dat)
		hup[hup<0] <- 0
		wmt <- wup
		hmt <- hup
		exf <- wmt%*%hmt
		csn <- sum((dat-exf)^2)/2
		div <- cst-ifelse(itr,csn,cst-tol*(10^ceiling(log10(prod(dim(dat))+1))))
		if (div>0) {
			wop <- wmt
			hop <- hmt
			}
		cst <- csn
		itr <- TRUE
		}
	wop <- wop[,order(apply(hop,1,sum),decreasing=TRUE)]
	hop <- hop[order(apply(hop,1,sum),decreasing=TRUE),]
	dimnames(wop) <- list(rownames(dat),1:k)
	dimnames(hop) <- list(1:k,colnames(dat))
	list(pos1=wop,pos2=hop)
	}

#' Latent Class Analysis
#'
#' A fast procedure for computing latent class analysis.  
#' @param dat Input data: can be a table or a data frame.
#' @param k Numeric specification of the number of latent classes to compute.
#' @param tol Numeric specification of the convergence criterion.
#' @param posterior Logical indicating whether the posterior probabilities of the individual observations should also be returned.
#' @param sep Character specifying the separator string for joining the levels (if \code{posterior = TRUE}).
#' @details
#' The prior probabilities of the latent classes are initialized with a Dirichlet distribution (by means of \code{rdirichlet} from
#'   the package \pkg{gtools}) with \code{alpha =} the total frequency counts of every level.  
#' @return A list with components:
#' \item{\code{prob0} }{The probabilities of the latent classes.}
#' \item{\code{prob1}-\code{prob...} }{The probabilities for each set of levels. The columns of each \code{prob...} sum to 1.}
#' \item{\code{posterior} }{If \code{posterior = TRUE}: An indicator matrix with the posterior probabilities of each observation.}
#' @references
#' Agresti, A. (2013) \emph{Categorical data analysis}. Hoboken: John Wiley and Sons, 535--542.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' lca.SndT_Fra <- fast_lca(SndT_Fra, k = 7)
#' lca.SndT_Fra
#' @export
fast_lca <- function(dat,k,tol=1e-8,posterior=FALSE,sep="_") {
	if (is.data.frame(dat)) {
		dat <- table(dat)
		}
	stopifnot(is.array(dat))
	way <- length(dim(dat))
	pri <- rep(1/k,k)
	tht <- lapply(1:way,function(j){gtools::rdirichlet(n=k,alpha=apply(dat,j,sum))})
	llk <- sum(as.vector(dat)*log(apply(sapply(1:k,function(z){outerec(lapply(1:way,function(j){tht[[j]][z,]}))*pri[z]}),1,sum)))
	div <- tol*10
	while (div>tol) {
		pst <- prop.table(sapply(1:k,function(z){outerec(lapply(1:way,function(j){tht[[j]][z,]}))*pri[z]}),1)
		tmp <- array(as.vector(dat)*pst,dim=c(dim(dat),k))
		tht <- lapply(1:way,function(j){matrix(prop.table(apply(tmp,c(way+1,j),sum),1),nrow=k)})
		pri <- prop.table(apply(tmp,way+1,sum))
		tmp <- 1
		lln <- sum(as.vector(dat)*log(apply(sapply(1:k,function(z){outerec(lapply(1:way,function(j){tht[[j]][z,]}))*pri[z]}),1,sum)))
		div <- lln-llk
		if (is.na(div)) {
			break
			}
		if (div>0) {
			pso <- pst
			tho <- tht
			pro <- pri
			}
		llk <- lln
		}
	tho <- lapply(1:way,function(j){t(matrix(tho[[j]][order(pro,decreasing=TRUE),],nrow=k))})
	names(tho) <- paste("prob",1:way,sep="",collapse=NULL)
	for (j in 1:way) {
		dimnames(tho[[j]]) <- list(dimnames(dat)[[j]],1:k)
		}
	out <- c(list(prob0=sort(pro,decreasing=TRUE)),tho)
	if (posterior) {
		pso <- pso[,order(pro,decreasing=TRUE)]
		ind <- tab2ind(dat,sep=sep,add_names=TRUE)
		tmp <- ind%*%pso
		colnames(tmp) <- 1:k
		out <- c(out,list(posterior=cbind(ind,tmp)))
		}
	out
	}

#' Probabilistic Latent Semantic Analysis
#'
#' A fast procedure for computing probabilistic latent semantic analysis.  
#' @param dat Input data: can be a table or a data frame (but the data frame must have only two columns).
#' @param k Numeric specification of the number of latent classes to compute.
#' @param symmetric Logical indicating whether to compute the symmetric or the asymmetric solution.
#' @param tol Numeric specification of the convergence criterion.
#' @details
#' From version 1.1.0 of the \pkg{svs} package on, probabilistic latent semantic analysis is a special case of latent class analysis.  
#' @return A list with components:
#' \item{\code{prob0} }{The probabilities of the latent classes.}
#' \item{\code{prob1} }{The probabilities for the first set of levels (\emph{viz.} the row levels of a frequency table). The rows of \code{prob1} sum to 1 if \code{symmetric = FALSE}, the columns sum to 1 if \code{symmetric = TRUE}.}
#' \item{\code{prob2} }{The probabilities for the second set of levels (\emph{viz.} the column levels of a frequency table). The columns of \code{prob2} sum to 1.}
#' @references
#' Hofmann, Th. (1999). Probabilistic latent semantic indexing.
#'   \emph{SIGIR'99: Proceedings of the 22nd annual international SIGIR conference on research and development in information retrieval}, 50--57.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' psa.SndT_Fra <- fast_psa(SndT_Fra, k = 7)
#' psa.SndT_Fra
#' @export
fast_psa <- function(dat,k,symmetric=FALSE,tol=1e-8) {
	if (is.data.frame(dat)) {
		dat <- table(dat)
		}
	stopifnot(is.matrix(dat))
	out <- fast_lca(dat=dat,k=k,tol=tol,posterior=FALSE)
	if (!symmetric) {
		out$prob1 <- prop.table(sweep(out$prob1,2,out$prob0,"*"),1)
		}
	out
	}

#' @rdname fast_psa
#' @export
fast_psi <- function(dat,k,symmetric=FALSE,tol=1e-8) {
	fast_psa(dat=dat,k=k,symmetric=symmetric,tol=tol)
	}

#' @rdname fast_psa
#' @export
fast_plsa <- function(dat,k,symmetric=FALSE,tol=1e-8) {
	fast_psa(dat=dat,k=k,symmetric=symmetric,tol=tol)
	}

#' @rdname fast_psa
#' @export
fast_plsi <- function(dat,k,symmetric=FALSE,tol=1e-8) {
	fast_psa(dat=dat,k=k,symmetric=symmetric,tol=tol)
	}

#' EM clustering
#'
#' A fast procedure for Expectation-Maximization clustering.  
#' @param dat Input data: can be a table or a data frame (but the data frame must have only two columns).
#' @param k Numeric specification of the number of latent classes to compute.
#' @param tol Numeric specification of the convergence criterion.
#' @details
#' This function assumes that the rows of a frequency table come from a multinomial distribution. The prior probabilities of
#'   the latent classes are initialized with a Dirichlet distribution (by means of \code{rdirichlet} from the package \pkg{gtools}) with
#'   \code{alpha =} the total frequency counts of every level.  
#' @return A list with components:
#' \item{\code{prob0} }{The probabilities of the latent classes.}
#' \item{\code{prob1} }{The probabilities for the first set of levels (\emph{viz.} the row levels of a frequency table). The rows of \code{prob1} sum to 1.}
#' \item{\code{prob2} }{The probabilities for the second set of levels (\emph{viz.} the column levels of a frequency table). The rows of \code{prob2} sum to 1.}
#' @references
#' Dempster, A. P., N. M. Laird and D. B. Rubin (1977) Maximum likelihood from incomplete data via the EM algorithm.
#'   \emph{Journal of the royal statistical society, series B} \strong{39} (1), 1--38.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' E_M.SndT_Fra <- fast_E_M(SndT_Fra, k = 7)
#' E_M.SndT_Fra
#' @export
fast_E_M <- function(dat,k,tol=1e-8) {
	if (is.data.frame(dat)) {
		dat <- table(dat)
		}
	stopifnot(is.matrix(dat))
	len <- nrow(dat)
	pri <- rep(1/k,k)
	tht <- gtools::rdirichlet(n=k,alpha=apply(dat,2,sum))
	llk <- sum(log(sapply(1:len,function(i){sum(sapply(1:k,function(z){max(dmultinom(x=dat[i,],size=sum(dat[i,]),prob=tht[z,]),1e-300)*pri[z]}))})))
	div <- tol*10
	while (div>tol) {
		pst <- prop.table(sapply(1:len,function(i){sapply(1:k,function(z){max(dmultinom(x=dat[i,],size=sum(dat[i,]),prob=tht[z,]),1e-300)*pri[z]})}),2)
		tht <- prop.table(pst%*%dat,1)
		pri <- apply(pst,1,sum)/len
		lln <- sum(log(sapply(1:len,function(i){sum(sapply(1:k,function(z){max(dmultinom(x=dat[i,],size=sum(dat[i,]),prob=tht[z,]),1e-300)*pri[z]}))})))
		div <- lln-llk
		if (is.na(div)) {
			break
			}
		if (div>0) {
			pso <- pst
			tho <- tht
			pro <- pri
			}
		llk <- lln
		}
	pso <- pso[order(pro,decreasing=TRUE),]
	tho <- tho[order(pro,decreasing=TRUE),]
	dimnames(pso) <- list(1:k,rownames(dat))
	dimnames(tho) <- list(1:k,colnames(dat))
	list(prob0=sort(pro,decreasing=TRUE),prob1=t(pso),prob2=tho)
	}

#' @rdname fast_E_M
#' @export
fast_EM <- function(dat,k,tol=1e-8) {
	fast_E_M(dat=dat,k=k,tol=tol)
	}

#' Weighting Functions
#'
#' Local and global weighting functions.  
#' @param x A numeric matrix.
#' @details
#' There are many local and global weighting functions. In this package, local weighting functions are prefixed with \code{lw_} and
#'   global weighting functions with \code{gw_}, so users can define their own weighting functions.  
#' 
#' Local weighting functions (i.e. weighting every cell in the matrix):
#' \itemize{
#'   \item{\code{lw_tf} }{Term frequency: \emph{f(x) = x}.}
#'   \item{\code{lw_raw} }{Raw frequency, which is the same as the term frequency: \emph{f(x) = x}.}
#'   \item{\code{lw_log} }{Logarithm: \emph{f(x) = log(x + 1)}.}
#'   \item{\code{lw_bin} }{Binary: \emph{f(x) = 1} if \emph{x > 0} and \emph{0} otherwise.}
#' }
#' 
#' Global weighting functions, weighting the columns of the matrix (hence, these weighting functions work according to expectation for
#'   a document-term matrix, i.e. with the documents as the rows and the terms as the columns):
#' \itemize{
#'   \item{\code{gw_idf} }{Inverse document frequency: \emph{f(x) = log( nrow(x) / n + 1)} where \emph{n =} the number of rows in which the column \emph{>0}.}
#'   \item{\code{gw_idf_alt} }{Alternative definition of the inverse document frequency: \emph{f(x) = log( nrow(x) / n) + 1} where \emph{n =} the number of rows in which the column \emph{>0}.}
#'   \item{\code{gw_gfidf} }{Global frequency multiplied by inverse document frequency: \emph{f(x) = colSums(x) / n} where \emph{n =} the number of rows in which the column \emph{>0}.}
#'   \item{\code{gw_nor} }{Normal(ized) frequency: \emph{f(x) = x / colSums(x^2)}.}
#'   \item{\code{gw_ent} }{Entropy: \emph{f(x) = 1 +} the relative Shannon entropy.}
#'   \item{\code{gw_bin} }{Binary: \emph{f(x) = 1}.}
#'   \item{\code{gw_raw} }{Raw, which is the same as binary: \emph{f(x) = 1}.}
#' }
#' @return A numeric matrix.
#' @seealso \code{\link{fast_lsa}}.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' tab.SndT_Fra <- table(SndT_Fra)
#' lw_log(tab.SndT_Fra)
#' gw_idf(tab.SndT_Fra)
#' @name weighting_functions
NULL

#' @rdname weighting_functions
#' @export
lw_tf <- function(x) {
	x
	}

#' @rdname weighting_functions
#' @export
lw_raw <- function(x) {
	x
	}

#' @rdname weighting_functions
#' @export
lw_log <- function(x) {
	log(x+1)
	 }

#' @rdname weighting_functions
#' @export
lw_bin <- function(x) {
	(x>0)*1
	}

#' @rdname weighting_functions
#' @export
gw_idf <- function(x) {
	log2(nrow(x)/(apply(x>0,2,sum)+1))
	}

#' @rdname weighting_functions
#' @export
gw_idf_alt <- function(x) {
	log2(nrow(x)/apply(x>0,2,sum))+1
	}

#' @rdname weighting_functions
#' @export
gw_gfidf <- function(x) {
	apply(x,2,sum)/apply(x>0,2,sum)
	}

#' @rdname weighting_functions
#' @export
gw_nor <- function(x) {
	1/sqrt(apply(x^2,2,sum))
	}

#' @rdname weighting_functions
#' @export
gw_ent <- function(x) {
	tmp <- prop.table(x,2)
	tmp <- tmp*log(tmp)
	tmp[!is.finite(tmp)] <- 0
	1-apply(tmp/log(nrow(x)),2,sum)
	}

#' @rdname weighting_functions
#' @export
gw_bin <- function(x) {
	rep(1,ncol(x))
	}

#' @rdname weighting_functions
#' @export
gw_raw <- function(x) {
	rep(1,ncol(x))
	}

#' Pointwise Mutual Information
#'
#' A function for computing the pointwise mutual information of every entry in a table.  
#' @param x A table (i.e. an object which can be converted to an array).
#' @param normalize Logical indicating whether to normalize the pointwise mutual information.
#' @param base Numeric specification of the base with respect to which logarithms are computed.
#' @return An array with the pointwise mutual information of every entry.
#' @seealso \code{\link{MI}}.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' tab.SndT_Fra <- table(SndT_Fra)
#' pmi(tab.SndT_Fra)
#' @export
pmi <- function(x,normalize=FALSE,base=2) {
	stopifnot(is.array(x))
	stopifnot(all(x>=0))
	x <- x/sum(x)
	out <- log(x/outerec(lapply(1:length(dim(x)),function(j){apply(x,j,sum)})),base=base)
	if (normalize) {
		out <- -out/log(x,base=base)
		}
	out[!is.finite(out)] <- 0
	out
	}

#' @rdname pmi
#' @export
PMI <- function(x,normalize=FALSE,base=2) {
	pmi(x=x,normalize=normalize,base=base)
	}

#' Mutual Information
#'
#' A function for computing the mutual information.  
#' @param x A table (i.e. an object which can be converted to an array).
#' @param base Numeric specification of the base with respect to which logarithms are computed.
#' @return A numeric value containing the mutual information.
#' @seealso \code{\link{pmi}}.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' tab.SndT_Fra <- table(SndT_Fra)
#' MI(tab.SndT_Fra)
#' @export
MI <- function(x,base=2) {
	stopifnot(is.array(x))
	stopifnot(all(x>=0))
	x <- x/sum(x)
	cel <- x*pmi(x,normalize=FALSE,base=base)
	cel[!is.finite(cel)] <- 0
	sum(cel)
	}

#' @rdname MI
#' @export
mi <- function(x,base=2) {
	MI(x=x,base=base)
	}

#' Compute Level Frequencies (for a Factor)
#'
#' A helper function for computing the frequency of each factor level (typically used in correspondence analysis).  
#' @param dat A factor or a data frame.
#' @param nfac Logical indicating whether the number of factors (i.e. the number of columns in \code{dat}) is a divisor for
#'   the level frequencies.
#' @return A vector containing the frequency counts of every level.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' freq_ca(SndT_Fra)
#' @export
freq_ca <- function(dat,nfac=FALSE) {
	if (is.factor(dat)) {
		dat <- as.data.frame(dat)
		}
	stopifnot(is.data.frame(dat))
	out <- lapply(dat,function(fac){summary(fac,maxsum=nlevels(fac))})
	names(out) <- NULL
	unlist(out)*ifelse(nfac,ncol(dat),1)
	}

#' Compute Coordinates for Cluster Centers
#'
#' A helper function for computing the coordinates of cluster centers (typically used in correspondence analysis).  
#' @param x A numeric matrix.
#' @param clusters A clustering of the row levels of \code{x}: either a list or the output of \code{kmeans}.
#' @param freq An optional vector of frequency counts for the row levels of \code{x}.
#' @return A matrix containing the coordinates of the cluster centers.
#' @seealso \code{\link{freq_ca}}.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.SndT_Fra <- fast_sca(SndT_Fra)
#' kcl.SndT_Fra <- kmeans(sca.SndT_Fra$pos1, centers = 7)
#' centers_ca(sca.SndT_Fra$pos1, clusters = kcl.SndT_Fra, freq = freq_ca(SndT_Fra[, 1]))
#' @export
centers_ca <- function(x,clusters,freq) {
	stopifnot(is.matrix(x))
	stopifnot(class(clusters) %in% c("list","kmeans"))
	if (class(clusters)=="kmeans") {
		clusters <- lapply(sort(unique(clusters$cluster)),function(y){names(clusters$cluster[clusters$cluster==y])})
		}
	if (is.null(names(clusters))) {
		names(clusters) <- 1:length(clusters)
		}
	if (is.null(names(freq))) {
		names(freq) <- rownames(x)
		}
	do.call(what=rbind,args=lapply(clusters,function(z){apply(sweep(rbind(x[z,]),1,freq[z],"*"),2,sum)/sum(freq[z])}))
	}

#' Compute Chi-square Distances
#'
#' A function for computing chi-square distances.  
#' @param x A numeric matrix (containing coordinates).
#' @param diag Logical specifying whether the diagonal of the resulting distance matrix should be printed.
#' @param upper Logical specifying whether the upper triangle of the resulting distance matrix should be printed.
#' @return A distance matrix.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' tab.SndT_Fra <- table(SndT_Fra)
#' dist_chisquare(tab.SndT_Fra)
#' @export
dist_chisquare <- function(x,diag=FALSE,upper=FALSE) {
	stopifnot(is.matrix(x))
	stopifnot(all(x>=0))
	dist(sweep(sweep(x,1,apply(x,1,sum),"/"),2,sqrt(apply(x,2,sum)/sum(x)),"/"),method="euclidean",diag=diag,upper=upper)
	}

#' @rdname dist_chisquare
#' @export
dist_chisq <- function(x,diag=FALSE,upper=FALSE) {
	dist_chisquare(x=x,diag=diag,upper=upper)
	}

#' Compute Cosine Distances
#'
#' A function for computing cosine distances.  
#' @param x A numeric matrix (containing coordinates).
#' @param diag Logical specifying whether the diagonal of the resulting distance matrix should be printed.
#' @param upper Logical specifying whether the upper triangle of the resulting distance matrix should be printed.
#' @details
#' The cosine distance equals 1 - the cosine similarity.  
#' @return A distance matrix.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' lsa.SndT_Fra <- fast_lsa(SndT_Fra)
#' dist_cosine(lsa.SndT_Fra$pos1[, 1:7])
#' @export
dist_cosine <- function(x,diag=FALSE,upper=FALSE) {
	stopifnot(is.matrix(x))
	as.dist(1-x%*%t(x)/sqrt(apply(x^2,1,sum)%o%apply(x^2,1,sum)),diag=diag,upper=upper)
	}

#' @rdname dist_cosine
#' @export
dist_cos <- function(x,diag=FALSE,upper=FALSE) {
	dist_cosine(x=x,diag=diag,upper=upper)
	}

#' Compute Distances with respect to a certain Point
#'
#' A function for computing (euclidean) distances with respect to a certain specified point.  
#' @param x A numeric matrix (containing coordinates).
#' @param wrt A specification of the point with respect to which to compute all distances: can be either a vector or
#'   the character label of one of the row levels in \code{x}. If \code{NULL} or \code{NA}, then the origin (i.e.
#'   the point \code{c(0 , 0, 0,... )}) is taken as the value.
#' @return A matrix (containing distances between the rows of \code{x} and \code{wrt}).
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.SndT_Fra <- fast_sca(SndT_Fra)
#' dist_wrt(sca.SndT_Fra$pos1, wrt = "beginnen")
#' @export
dist_wrt <- function(x,wrt=NULL) {
	stopifnot(is.matrix(x))
	if (is.null(wrt) || is.na(wrt)) {
		wrt <- rep(0,times=ncol(x))
		}
	if (is.character(wrt)) {
		stopifnot(length(wrt)==1)
		wrt <- x[wrt,]
		}
	stopifnot(length(wrt)==ncol(x))
	sqrt(apply(sweep(x,2,wrt,"-")^2,1,sum))
	}

#' Compute Distances with respect to Cluster Centers
#'
#' A function for computing (euclidean) distances with respect to specified cluster centers.  
#' @param x A numeric matrix (containing coordinates).
#' @param clusters A clustering of the row levels of \code{x}: either a list or the output of \code{kmeans}.
#' @param freq An optional vector of frequency counts for the row levels of \code{x}.
#' @param members_only Logical specifying whether the distances from the cluster centers should only be computed for
#'   the cluster members.
#' @return A list with a matrix of distances for every cluster.
#' @seealso \code{\link{centers_ca}}, \code{\link{freq_ca}}.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.SndT_Fra <- fast_sca(SndT_Fra)
#' kcl.SndT_Fra <- kmeans(sca.SndT_Fra$pos1, centers = 7)
#' dist_wrt_centers(sca.SndT_Fra$pos1, clusters = kcl.SndT_Fra, freq = freq_ca(SndT_Fra[, 1]))
#' @export
dist_wrt_centers <- function(x,clusters,freq=NULL,members_only=TRUE) {
	stopifnot(is.matrix(x))
	stopifnot(class(clusters) %in% c("list","kmeans"))
	if (class(clusters)=="kmeans") {
		clusters <- lapply(sort(unique(clusters$cluster)),function(y){names(clusters$cluster[clusters$cluster==y])})
		}
	if (all(is.null(freq))) {
		cen <- lapply(clusters,function(z){apply(rbind(x[z,]),2,mean)})
		}
	else {
		if (is.null(names(freq))) {
			names(freq) <- rownames(x)
			}
		cen <- lapply(clusters,function(z){apply(sweep(rbind(x[z,]),1,freq[z],"*"),2,sum)/sum(freq[z])})
		}
	if (members_only) {
		out <- lapply(1:length(clusters),function(k){dist_wrt(matrix(x[clusters[[k]],],ncol=ncol(x),dimnames=list(clusters[[k]])),cen[[k]])})
		}
	else {
		out <- lapply(1:length(clusters),function(k){dist_wrt(x,cen[[k]])})
		}
	out
	}

#' Plotting a Cumulative Distribution
#'
#' A function for plotting a cumulative distribution.  
#' @param x A numeric vector.
#' @param inc The (numeric) increment for constructing the sequence from 0 to \code{ceiling(max(x))}, plotted on the horizontal axis. 
#' @param col The color of the line and the text labels: see \code{\link[grDevices]{colors}}.
#' @param cex The character expansion factor: a numberic value to specify the size of the text labels.
#' @param font The font of the text labels: \code{1} for plain, \code{2} for bold, \code{3} for italic, and \code{4} for bold italic.
#' @param family The font family of the text labels: \code{"serif"}, \code{"sans"}, \code{"mono"}, or one of the \code{\link[grDevices]{Hershey}} fonts.
#' @param srt The rotation angle (in degrees) of the text labels.
#' @param pch The plotting character for displaying points: see \code{\link[graphics]{points}}.
#' @param pcol The color of the plotting character: see \code{\link[grDevices]{colors}}.
#' @param pbg The background color of the plotting character: see \code{\link[grDevices]{colors}}.
#' @param pcex The character expansion factor of the plotting character: a numeric value to specify the size of the plotting character.
#' @param lcol The color of the line: see \code{\link[grDevices]{colors}}.
#' @param lwd The line width of the line: a numeric value to specify the width of the line.
#' @param lty The line type of the line: \code{0} or \code{"blank"}, \code{1} or \code{"solid"}, \code{2} or \code{"dashed"}, \code{3} or
#'   \code{"dotted"}, \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or \code{"twodash"}.
#' @param xlim A vector of two numeric values specifying the lower and upper limit between which to plot the horizontal axis.
#' @param ylim A vector of two numeric values specifying the lower and upper limit between which to plot the vertical axis.
#' @param xlab A character string for labelling of the horizontal axis.
#' @param ylab A character string for labelling of the vertical axis.
#' @param main A character string for the main title of the plot.
#' @param sub A character string for the subtitle of the plot.
#' @return A cumulative distribution plot.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.SndT_Fra <- fast_sca(SndT_Fra)
#' dis.SndT_Fra <- dist_wrt(sca.SndT_Fra$pos1)
#' cd_plot(dis.SndT_Fra)
#' @export
cd_plot <- function(x,inc=0.01,col="darkgrey",cex=1,font=1,family="",srt=-45,pch=20,pcol="black",pbg="white",pcex=cex,lcol=col,lwd=1,lty=1,xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,main=NULL,sub=NULL) {
	stopifnot(is.numeric(x))
	s <- sort(x)
	r <- seq(0,ceiling(max(s)),by=inc)
	if (is.null(ylim) || is.na(ylim)) {
		ylim <- c(0,length(s)+1)
		}
	xlab <- ifelse(is.null(xlab),NA,xlab)
	ylab <- ifelse(is.null(ylab),NA,ylab)
	plot(x=r,y=sapply(r,function(z){sum(s<z)}),type="l",col=lcol,lwd=lwd,lty=lty,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,sub=sub)
	points(x=s,y=1:(length(s)),pch=pch,col=pcol,bg=pbg,cex=pcex)
	text(x=s,y=1:(length(s)),labels=names(s),xpd=TRUE,adj=c(1.1,0.5),col=col,cex=cex,font=font,family=family,srt=srt)
	invisible()
	}

#' Plotting Parallel Coordinates
#'
#' A function for plotting parallel coordinates.  
#' @param x A numeric matrix.
#' @param col The color of the text labels, points and connecting lines: see \code{\link[grDevices]{colors}}.
#' @param cex The character expansion factor: A numeric value to specify the size of the text labels and the points.
#' @param font The font of the text labels: \code{1} for plain, \code{2} for bold, \code{3} for italic, and \code{4} for bold italic.
#' @param family The font family of the text labels: \code{"serif"}, \code{"sans"}, \code{"mono"}, or one of the \code{\link[grDevices]{Hershey}} fonts.
#' @param pch The plotting character for displaying points: see \code{\link[graphics]{points}}.
#' @param pcol The color of the plotting character: see \code{\link[grDevices]{colors}}.
#' @param pcex The character expansion factor of the plotting character: a numeric value to specify the size of the plotting character.
#' @param lcol The color of the connecting lines: see \code{\link[grDevices]{colors}}.
#' @param lwd The line width of the connecting lines: a numeric value to specify the width of the connecting lines.
#' @param lty The line type of the connecting lines: \code{0} or \code{"blank"}, \code{1} or \code{"solid"}, \code{2} or \code{"dashed"}, \code{3} or
#'   \code{"dotted"}, \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or \code{"twodash"}.
#' @param acol The color of the parallel axes: see \code{\link[grDevices]{colors}}.
#' @param alwd The line width of the parallel axes: a numeric value to specify the width of the parallel axes.
#' @param alty The line type of the parallel axes: \code{0} or \code{"blank"}, \code{1} or \code{"solid"}, \code{2} or \code{"dashed"}, \code{3} or
#'   \code{"dotted"}, \code{4} or \code{"dotdash"}, \code{5} or \code{"longdash"}, \code{6} or \code{"twodash"}.
#' @param las The reading direction of the labels on the axes ("label axis style"): either a numeric value between \code{0} and \code{3} (see \code{las} in
#'   \code{\link[graphics]{par}}), or a character value matching either \code{"horizontal"} or \code{"vertical"}.
#' @param add_scale Logical specifying whether to add a scale for the parallel axes (which are normalized).
#' @param main A character string for the main title of the plot.
#' @param sub A character string for the subtitle of the plot.
#' @return A parallel coordinate plot.
#' @examples
#' SndT_Fra <- read.table(system.file("extdata", "SndT_Fra.txt", package = "svs"),
#'    header = TRUE, sep = "\t", quote = "\"", encoding = "UTF-8")
#' sca.SndT_Fra <- fast_sca(SndT_Fra)
#' pc_plot(sca.SndT_Fra$pos1, las = "vertical")
#' @export
pc_plot <- function(x,col="darkgrey",cex=1,font=1,family="",pch=20,pcol=col,pcex=cex,lcol=col,lwd=1,lty=1,acol="black",alwd=1,alty=1,las=1,add_scale=FALSE,main=NULL,sub=NULL) {
	stopifnot(is.matrix(x))
	if (is.character(las)) {
		las <- pmatch(tolower(las),table=c("horizontal","vertical"))
		}
	len <- ncol(x)
	ran <- apply(x,2,range,na.rm=TRUE)
	loc <- sweep(sweep(x,2,ran[1,],"-"),2,(ran[2,]-ran[1,]),"/")
	mai <- par("mai")
	nai <- mai
	nai[2] <- max(strwidth(rownames(x),units="inches",cex=cex,font=font))+0.3
	if (las>1) {
		nai[1] <- max(c(mai[1]-0.3,strwidth(colnames(x),units="inches",cex=cex,font=font)))+0.3
		hdj <- 1.2-ifelse(cex>1,cex/10,0)
		pdj <- 0.5
		}
	else {
		hdj <- 0.5
		pdj <- 1.1
		}
	par(mai=nai)
	matplot(t(loc),type="n",xlim=c(1,len),ylim=c(0,1),xlab="",ylab="",main=main,sub=sub,axes=FALSE,add=FALSE)
	axis(side=1,at=1:len,labels=FALSE)
	mtext(text=colnames(x),side=1,at=1:len,las=las,adj=hdj,padj=pdj,col=acol,cex=cex,font=font,family=family)
	mtext(text=rownames(x),side=2,at=loc[,1],las=1,adj=NA,padj=NA,col=col,cex=cex,font=font,family=family)
	abline(v=1:len,col=acol,lty=alty,lwd=alwd)
	matlines(t(loc),col=lcol,lty=lty,lwd=lwd)
	pch <- ifelse(is.null(pch),NA,pch)
	if (!is.na(pch)) {
		matpoints(t(loc),pch=pch,col=pcol,cex=pcex)
		}
	if (add_scale) {
		axis(side=4,col.axis=acol,cex.axis=cex,font.axis=font,las=1)
		}
	par(mai=mai)
	invisible()
	}

#' Complete the Output of \code{pvpick}
#'
#' A helper function to add the missing singleton clusters in the output of \code{pvpick} (from the package \pkg{pvclust}).  
#' @param clusters A clustering by a call to \code{pvpick}.
#' @param labels A character vector containing the exhaustive set of levels.
#' @return A list with the singleton clusters inserted at the end (so that the set of clusters is exhaustive).
#' @export
complete_pvpick <- function(clusters,labels) {
	if (!is.null(names(clusters)) && names(clusters)[1]=="clusters") {
		clusters <- clusters$clusters
		}
	c(clusters,as.list(labels[!labels %in% unlist(clusters)]))
	}
