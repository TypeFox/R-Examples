#' Miscellaneous Functions for Community Ecology
#'
#' This is just a small collection of miscellaneous functions
#' that may be useful, primarily for community ecology analyses,
#' particularly for paleoecological data. They are here mainly for
#' pedagogical reasons (i.e. for students) as they don't appear
#' to be available in other ecology-focused packages.

#' @param x The community abundance matrix. Must be a matrix with two dimensions for
#' \code{pairwiseSpearmanRho}; if a vector is supplied to \code{HurlbertPIE}, then
#' it is treated as if it was matrix with a single row and number of columns equal to
#' its length. Taxonomic units are assumed to be the columns and sites (samples)
#'  are assumed to be the rows, for both functions.

#' @param dropAbsent Should absent taxa be dropped? Must be one of either 'bothAbsent' (drop taxa
#' absent in both sites for a given pairwise comparison),'eitherAbsent' (drop taxa absent in either
#' site), or 'noDrop' (drop none of the taxa). The default 'bothAbsent' is recommended, see examples.

#' @param asDistance Should the rho coefficients be rescaled on a scale similar to
#' dissimilarity metrics, i.e. bounded 0 to 1, with 1 representing maximum dissimilarity (i.e.
#' a Spearman rho correlation of -1)? ( dissimilarity = (1 - rho) / 2 )

#' @param diag Should the diagonal of the output distance matrix be included?

#' @param upper Should the upper triangle of the output distance matrix be included?

#' @param na.rm Should taxa listed with NA values be dropped from a pair-wise site comparison?
#' If FALSE, the returned value for that site pair will be NA if NAs are present.

#' @param nAnalyze Allows users to select that PIE be calculated only on \code{nAnalyze} most
#' abundant taxa in a site sample. \code{nAnalyze} must be a vector of length 1, consisting of
#' whole-number value at least equal to 2. By default, \code{nAnalyze = Inf} so all taxa are
#' accepted. Note that if there are less taxa in a sample than nAnalyze, the number present will
#' be used.

#' @details
#' \code{pairwiseSpearmanRho} returns Spearman rho correlation coefficients
#' based on the rank abundances of taxa (columns) within sites (rows) from
#' the input matrix, by internally wrapping the function \code{cor.test}.
#' It allows for various options that ultimatically allow
#' for dropping taxa not shared between two sites (the default), as well as
#' several other options. This allows the rho coefficient to behave like the
#' Bray-Curtis distance, in that it is not affected by the number of taxa absent
#' in both sites.
#' 
#' \code{pairwiseSpearmanRho} can also rescale the rho coefficients with (1-rho)/2
#' to provide a measure similar to a dissimilarity metric, bounded between 0 and 1.
#' This function was written so several arguments would be in a similar format to
#' the \code{vegan} library function \code{vegdist}. If used to obtain rho
#' rescaled as a dissimilarity, the default output will be the lower triangle of
#' a distance matrix object, just as is returned by default by \code{vegdist}.
#' This behavior can be modified via the arguments for including the diagonal 
#' and upper triangle of the matrix. Otherwise, a full matrix is returned (by default)
#' if the \code{asDistance} argument is not enabled. 
#'
#' \code{HurlbertPIE} provides the Probability of Interspecific Encounter metric for
#' relative community abundance data, a commonly used metric for evenness of community
#' abundance data based on derivations in Hurlbert (1971). An optional argument allows
#' users to apply Hurlbert's PIE to only a subselection of the most abundant taxa.

#' @return
#' \code{pairwiseSpearmanRho} will return either a full matrix (the default) or (if
#' \code{asDistance} is true, a distance matrix, with only the lower triangle
#' shown (by default). See details.
#'
#' \code{HurlbertPIE} returns a named vector of PIE values for the input data.

#' @aliases communityEcology pairwiseSpearmanRho HurlbertPIE

#' @seealso
#' Example dataset: \code{\link{kanto}}

#' @name communityEcology

#' @author David W. Bapst

#' @references
#' Hurlbert, S. H. 1971. The nonconcept of species diversity:
#' a critique and alternative parameters. \emph{Ecology} 52(4):577-586.

#' @examples
#'
#' # let's load some example data:
#' # a classic dataset collected by Satoshi and Okido from the Kanto region
#' 
#' data(kanto)
#' 
#' rhoBothAbsent<-pairwiseSpearmanRho(kanto,dropAbsent="bothAbsent")
#' 
#' #other dropping options
#' rhoEitherAbsent<-pairwiseSpearmanRho(kanto,dropAbsent="eitherAbsent")
#' rhoNoDrop<-pairwiseSpearmanRho(kanto,dropAbsent="noDrop")
#' 
#' #compare
#' layout(1:3)
#' lim<-c(-1,1)
#' plot(rhoBothAbsent, rhoEitherAbsent, xlim=lim, ylim=lim)
#' 	abline(0,1)
#' plot(rhoBothAbsent, rhoNoDrop, xlim=lim, ylim=lim)
#' 	abline(0,1)
#' plot(rhoEitherAbsent, rhoNoDrop, xlim=lim, ylim=lim)
#' 	abline(0,1)
#' layout(1)
#' 
#' #using dropAbsent="eitherAbsent" reduces the number of taxa so much that
#' 	# the number of taxa present drops too low to be useful
#' #dropping none of the taxa restricts the rho measures to high coefficients
#' 	# due to the many shared 0s for absent taxa
#' 
#' #############
#' 
#' # Try the rho coefficients as a rescaled dissimilarity
#' rhoDist<-pairwiseSpearmanRho(kanto,asDistance=TRUE,dropAbsent="bothAbsent")
#' 
#' # What happens if we use these in typical distance matrix based analyses?
#' 
#' # Cluster analysis
#' clustRes<-hclust(rhoDist)
#' plot(clustRes)
#' 
#' # Principle Coordinates Analysis
#' pcoRes <- pcoa(rhoDist,correction="lingoes")
#' scores <- pcoRes$vectors
#' #plot the PCO
#' plot(scores,type="n")
#' text(labels=rownames(kanto),scores[,1],scores[,2],cex=0.5)
#' 
#' ##################################
#' 
#' # measuring evenness with Hurlbert's PIE
#' 
#' kantoPIE<-HurlbertPIE(kanto)
#' 
#' #histogram
#' hist(kantoPIE)
#' #evenness of the kanto data is fairly high
#' 
#' #barplot
#' parX<-par(mar=c(7,5,3,3))
#' barplot(kantoPIE,las=3,cex.names=0.7,
#' 	ylab="Hurlbert's PIE",ylim=c(0.5,1),xpd=FALSE)
#' par(parX)
#' 
#' #and we can see that the Tower has extremely low unevenness
#' 	#...overly high abundance of ghosts?
#' 
#' #let's look at evenness of 5 most abundant taxa
#' kantoPIE_5<-HurlbertPIE(kanto,nAnalyze=5)
#'
#' #barplot
#' parX<-par(mar=c(7,5,3,3))
#' barplot(kantoPIE_5,las=3,cex.names=0.7,
#' 	ylab="Hurlbert's PIE for 5 most abundant taxa",ylim=c(0.5,1),xpd=FALSE)
#' par(parX)

#' @rdname communityEcology
#' @export
pairwiseSpearmanRho<-function(x, dropAbsent="bothAbsent", asDistance=FALSE,
	diag=NULL, upper=NULL, na.rm=FALSE){
	#for ecology: SPECIES ARE COLUMNS, SAMPLES ARE ROWS
	if(length(dim(x))!=2){
		stop('x must be a matrix, with samples as rows')
		}
	dropPar<-c('bothAbsent','eitherAbsent','noDrop')
	dropAbsent<-dropPar[pmatch(dropAbsent,dropPar)]
	if(asDistance){
		if(is.null(diag)){diag<-FALSE}
		if(is.null(upper)){upper<-FALSE}
	}else{
		if(is.null(diag)){diag<-TRUE}
		if(is.null(upper)){upper<-TRUE}	
		}
	if(is.na(dropAbsent)){
		stop(paste0('dropAbsent must be one of ',
			paste(dropPar,collapse=', ')))}
	rhos<-matrix(,nrow(x),nrow(x))
	for(i in 1:nrow(x)){
		for(j in 1:nrow(x)){
			if(i>=j){
				#if no dropping NA and there are any NAs
				if(!na.rm & any(is.na(x[c(i,j),]))){
					#just return NA			
					rhos[i,j]<-rhos[j,i]<-NA
				}else{
					if(na.rm){
						selector<-!is.na(x[i,]) & !is.na(x[j,])
					}else{
						selector<-rep(TRUE,ncol(x))
						}
					#now, if dropping absent, add to selector
					if(dropAbsent=='bothAbsent'){
						selector<-selector & (x[i,]!=0 | x[j,]!=0)
						}
					if(dropAbsent=='eitherAbsent'){
						selector<-selector & (x[i,]!=0 & x[j,]!=0)
						}
					#check selector
					if(sum(selector)<2){
						#then these two samples aren't comparable, return NA						
						rhos[i,j]<-rhos[j,i]<-NA
					}else{
						samp1<-x[i,selector]
						samp2<-x[j,selector]
						#print(sum(selector))
						rho<-suppressWarnings(
							cor.test(samp1,samp2,method='spearman')$estimate
							)
						rhos[i,j]<-rhos[j,i]<-rho
						}
					}
				}
			}
		}
	colnames(rhos)<-rownames(rhos)<-rownames(x)
	if(asDistance){
		rhos<-(1-rhos)/2
		result<-as.dist(rhos)
		attr(result, 'Diag') <- diag
		attr(result, 'Upper') <- upper
	}else{
		result<-rhos
		}
	return(result)
	}

#' @rdname communityEcology
#' @export
HurlbertPIE<-function(x,nAnalyze=Inf){
	if(is.vector(x)){ 
		x<-matrix(x,1,length(x))
		}
	if(length(nAnalyze)!=1 | !is.numeric(nAnalyze) | nAnalyze<2){
		stop("nAnalyze must be a numeric vector of length 1, with a value at least equal to 2")
		}
	if(!is.infinite(nAnalyze)){if(!is.integer(nAnalyze)){
		nAnalyze2<-as.integer(nAnalyze)
		if(nAnalyze2!=nAnalyze){
			stop("nAnalyze must be a whole number")
			}
		nAnalyze<-nAnalyze2
		}}
	PIE<-numeric()
	for(i in 1:nrow(x)){
		#first need to test there is actually more than one species
		samp<-x[i,]
		#remove all but n most abundant
		if(!is.infinite(nAnalyze)){
			samp<-samp[rank(-samp)<(nAnalyze+1)]
			}
		#drop zeroes
		samp<-samp[samp>0]
		#first need to test there is actually more than one species
		diversity<-length(samp)
		if(diversity>1){
			PIE[i]<-diversity/(diversity-1) * (1-sum((samp/sum(samp))^2))
		}else{
			PIE[i]<-0
			}
		}
	names(PIE)<-rownames(x)
	return(PIE)
	}

