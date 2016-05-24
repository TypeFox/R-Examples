#' Computes (average) Idenity-by-State for a set of people and markers
#' 
#' Given a set of SNPs, computes a matrix of average IBS for a group 
#' of individuals. 
#' This function facilitates quality control of genomic data. 
#' E.g. people with exteremly high (close to 1) IBS may indicate 
#' duplicated samples (or twins), simply high values of IBS may 
#' indicate relatives. 
#' 
#' When weight "freq" is used, IBS for a pair of people i and j is computed as
#' 
#' \deqn{
#' f_{i,j} = \frac{1}{N} \Sigma_k \frac{(x_{i,k} - p_k) * (x_{j,k} - p_k)}{(p_k * (1 - p_k))}
#' }
#' 
#' where k changes from 1 to N = number of SNPs GW, \eqn{x_{i,k}} is 
#' a genotype of ith person at the kth SNP, coded as 0, 1/2, 1 and 
#' \eqn{p_k} is the frequency 
#' of the "+" allele. This apparently provides an unbiased estimate of the 
#' kinship coefficient.
#' 
#' With "eVar" option above formula changes by using ( 2 * empirical variance 
#' of the genotype ) in the denominator. The empirical variance is computed 
#' according to the formula 
#' 
#' \deqn{
#' Var(g_k) = \frac{1}{M} \Sigma_i g_{ik}^2 - E[g_k]^2
#' }
#' 
#' where M is the number of people
#' 
#' Only with "freq" option monomorphic SNPs are regarded as non-informative.
#' 
#' ibs() operation may be very lengthy for a large number of people.
#'  
#' @return 	A (Npeople X Npeople) matrix giving average IBS (kinship) values
#' between a pair below the diagonal and number of SNP genotype  
#' measured for both members of the pair above the diagonal. 
#' 
#' On the diagonal, homozygosity 0.5*(1+inbreeding) is provided with option 
#' 'freq';  with option 'eVar' the diagonal is set to 0.5; the diagonal is set 
#' to homozygosity with option 'no'.  
#' 
#' attr(computedobject,"Var") returns variance (replacing the 
#' diagonal when the object is used by \code{\link{egscore}}
#' 
#' @param data object of \code{\link{snp.data-class}}
#' @param snpsubset Index, character or logical vector with subset of SNPs to run analysis on. 
#' If missing, all SNPs from \code{data} are used for analysis.
#' @param idsubset IDs of people to be analysed. 
#' If missing, all people from \code{data} are used for analysis.
#' @param cross.idsubset Parameter allowing parallel implementation. Not to be used normally.
#' If supplied together with idsubset, the ibs/kinship for all pairs between 
#' idsubset and cross.idsubset computed. 
#' @param weight "no" for direct IBS computations, "freq" to weight by allelic frequency 
#' asuming HWE and "eVar" for empirical variance to be used
#' @param snpfreq when option weight="freq" used, you can provide 
#' fixed allele frequencies
#' 
#' @author Yurii Aulchenko
#' 
#' @seealso \code{\link{check.marker}}, \code{\link{summary.snp.data}},
#' \code{\link{snp.data-class}}
#' 
#' @examples 
#' data(ge03d2c)
#' set.seed(7)
#' # compute IBS based on a random sample of 1000 autosomal marker
#' selectedSnps <- sample(autosomal(ge03d2c),1000,replace=FALSE)
#' a <- ibs(ge03d2c,snps=selectedSnps)
#' a[1:5,1:5]
#' mds <- cmdscale(as.dist(1-a))
#' plot(mds)
#' # identify smaller cluster of outliers
#' km <- kmeans(mds,centers=2,nstart=1000)
#' cl1 <- names(which(km$cluster==1))
#' cl2 <- names(which(km$cluster==2))
#' if (length(cl1) > length(cl2)) cl1 <- cl2;
#' cl1
#' # PAINT THE OUTLIERS IN RED
#' points(mds[cl1,],pch=19,col="red")
#' # compute genomic kinship matrix to be used with e.g. polygenic, mmscore, etc
#' a <- ibs(ge03d2c,snps=selectedSnps,weight="freq")
#' a[1:5,1:5]
#' # now replace diagonal with EIGENSTRAT-type of diaganal to be used for egscore
#' diag(a) <- hom(ge03d2c[,autosomal(ge03d2c)])$Var
#' a[1:5,1:5]
#' ##############################
#' # compare 'freq' with 'eVar'
#' ##############################
#' ibsFreq <- ibs(ge03d2c,snps=selectedSnps, weight="freq") 
#' ibsEvar <- ibs(ge03d2c,snps=selectedSnps, weight="eVar")
#' mdsEvar <- cmdscale( as.dist( 0.5 - ibsEvar ) )
#' plot(mdsEvar)
#' outliers <- (mdsEvar[,1]>0.1)
#' ibsFreq[upper.tri(ibsFreq,diag=TRUE)] <- NA
#' ibsEvar[upper.tri(ibsEvar,diag=TRUE)] <- NA
#' plot(ibsEvar,ibsFreq)
#' points(ibsEvar[outliers,outliers],ibsFreq[outliers,outliers],col="red")
#' 
#' @keywords htest
#' 
#' '

"ibs" <- 
		function (data,snpsubset,idsubset=NULL,cross.idsubset=NULL,weight="no",snpfreq=NULL) {
# idsubset, cross.idsubset: should be real names, not indexes!
	if (is(data,"gwaa.data")) data <- gtdata(data)
	if (!is(data,"snp.data")) stop("The data argument must be of snp.data-class or gwaa.data-class")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (is.null(idsubset) && !is.null(cross.idsubset)) stop("cross.idsubset arg cannot be used (idsubset missing)",immediate. = TRUE)
	if (!is.null(snpfreq)) {
		if (length(snpfreq) != data@nsnps) stop("snpfreq argument not equal in length to the number of SNPs in data")
		if (any(snpfreq<0.) || any(snpfreq>1.)) stop("snpfreq argument: frequencies out of [0,1]")
		if (!is(snpfreq,"numeric")) stop("snpfreq argument: non-numeric class")
	} else {
		snpfreq <- summary(data)[,"Q.2"]
	}
	
	wargs <- c("no","freq","homo","eVar")
	if (!(match(weight,wargs,nomatch=0)>0)) {
		out <- paste("weight argument should be one of",wargs,"\n")
		stop(out)
	}
#	if (npairs > 2500000) stop("Too many pairs to evaluate... Stopped")
	if (weight == "no") {
		homodiag <- hom(data)[,"Hom"]
		option = 0
	} else if (weight == "freq") {
		homodiag <- 0.5*(1+hom(data,snpfreq=snpfreq)[,"F"])
		option = 1
	} else if (weight == "homo") {
		smr <- summary(data)
		if (any(smr[,"P.12"] != 0 )) warning("'weight=homo' used, but data contain heterozygotes")
		rm(smr);gc();
		homodiag <- rep(1.0,nids(data))
		option = 2
	} else if (weight == "eVar") {
		homodiag <- rep(0.5,nids(data))
		option = 3
	} else {
		stop("'weight' argument value not recognised")
	}
	varidiag <- hom(data)[,"Var"]
	ibs.C.option <- 0
	if (!is.null(idsubset) && !(is.numeric(idsubset) || is.logical(idsubset) || is.character(idsubset))) stop("idsubset must be numeric, logical, or character")
	if (!is.null(cross.idsubset) && !(is.numeric(cross.idsubset) || is.logical(cross.idsubset) || is.character(cross.idsubset))) stop("cross.idsubset must be numeric, logical, or character")
	if (!is.null(idsubset) && (is.numeric(idsubset) || is.logical(idsubset))) idsubset <- data@idnames[idsubset]
	if (!is.null(cross.idsubset) && (is.numeric(cross.idsubset) || is.logical(cross.idsubset))) cross.idsubset <- data@idnames[cross.idsubset]
	if (!is.null(idsubset) && !is.null(cross.idsubset)) {
		idset1 <- idsubset
		idset2 <- cross.idsubset
		if (any(idset1 %in% idset2)) stop("idsubset and cross.idsubset should not overlap!")
		idsorder <- c(idset1,idset2)
		homodiag <- homodiag[match(idsorder,data@idnames)]
		varidiag <- varidiag[match(idsorder,data@idnames)]
		if (length(idsorder) != data@nids) data <- data[idsorder,]
		if (any(idsorder!=data@idnames)) data <- data[idsorder,]
		ibs.C.option <- 1
	} else if (!is.null(idsubset) & is.null(cross.idsubset)) {
		idset1 <- idsubset
		idset2 <- idsubset
		idsorder <- idsubset
		homodiag <- homodiag[match(idsorder,data@idnames)]
		varidiag <- varidiag[match(idsorder,data@idnames)]
		data <- data[idsorder,]
	} else if (is.null(idsubset) & is.null(cross.idsubset)) {
		idset1 <- data@idnames
		idset2 <- data@idnames
	} else {
		stop("can not be: impossible combination of idsubset and cross.idsubset")
	}
	gc()
	idset1.num <- which(data@idnames %in% idset1)
	idset2.num <- which(data@idnames %in% idset2)
	
	if (ibs.C.option==1) {
		if (option<=1) {
			sout <- .C("ibspar", as.raw(data@gtps), as.integer(data@nids), as.integer(data@nsnps), 
					as.integer(length(idset1.num)), as.integer(idset1.num-1), 
					as.integer(length(idset2.num)), as.integer(idset2.num-1), as.double(snpfreq), 
					as.integer(option), sout = double(2*length(idset1.num)*length(idset2.num)), 
					PACKAGE="GenABEL")$sout
		} else if (option>1) {
			stop("parallel execution facilities with options >1 only available in optimized version")
		} else {
			stop("unknown 'weight' option")
		}
		out <- list()
		out$ibs <- sout[1:(length(idset1.num)*length(idset2.num))]
		out$num <- sout[(length(idset1.num)*length(idset2.num)+1):(length(idset1.num)*length(idset2.num)*2)]
		dim(out$ibs) <- c(length(idset2.num),length(idset1.num))
		dim(out$num) <- c(length(idset1.num),length(idset2.num))
		rownames(out$ibs) <- idset2
		colnames(out$ibs) <- idset1
		rownames(out$num) <- idset1
		colnames(out$num) <- idset2
	} else if (ibs.C.option==0) {
		out <- .C("ibsnew", as.raw(data@gtps), as.integer(data@nids), as.integer(data@nsnps), 
				as.double(snpfreq), as.integer(option), sout = double(data@nids*data@nids), PACKAGE="GenABEL")$sout
		dim(out) <- c(length(idset2.num),length(idset1.num))
		out <- t(out)
		diag(out) <- homodiag
		rownames(out) <- idset1
		colnames(out) <- idset2
	} else {
		stop("can not be: incorrect ibs.C.option")
	}
	attributes(out) <- c(attributes(out),list(Var=varidiag))
	out
}
