#' Estimate Sampling Rate from Sampling Horizon Data (Solow and Smith, 1997)
#'
#' This function implements the exact maximum likelihood estimator for the
#' instantaneous sampling rate from Solow and Smith (1997, Paleobiology),
#' which is based on the relationship between the number of collections for a
#' set of taxa and their durations (known precisely in continuous time).

#' @details 
#' Given a dataset of taxa with a vector \eqn{N}, representing the number of
#' collections for each taxon, and a vector \eqn{D}, giving the precise duration
#' for each taxon, we can use the following maximum likelihood estimator from
#' Solow and Smith (1997) to obtain the instantaneous sampling rate:
#'
#' \eqn{samplingRate = (sum(N-1)^2)/(sum(D)*sum(N))}
#'
#' This method is exclusively for datasets with very precisely dated horizons,
#' such as microfossils from deep sea cores with very precise age models. The
#' first and last appearance must be known very precisely to provide an equally
#' precise estimate of the duration. Most datasets are not precise enough
#' for this method, due to chronostratigraphic uncertainty. However, note that the age
#' of individual collections other than the first and last appearance dates
#' do not need to be known: its only the number of collections that matters.

#' @param sampOcc A list with the number of elements equal to the number of taxa,
#' and each element of the list being a numerical vector with the length equal
#' to the number of collections for each taxon, and each value equal to the
#' precise date of that fossil's time of collection. These dates do not need
#' to be ordered. If not supplied, the elements \code{durations} and \code{nCollections} must
#' be supplied.

#' @param durations A vector of precise durations in continuous time, with the
#' length equal to the number of taxa. If not supplied, this is calculated from
#' \code{SampOcc}, which must be supplied.

#' @param nCollections A vector of integers representing the number of
#' collections for each taxon in the input durations. If not supplied
#' this is calculated from \code{SampOcc}, which must be supplied.

#' @return
#' Returns the instantaneous sampling (in per lineage*time-units) as a
#' single numerical value. Note that this is the instantaneous sampling
#' rate and not the probability of sampling a taxon per interval.

#' @seealso
#' Duration frequency methods (Foote and Raup, 1996; Foote, 1997) use
#' ranges alone to estimate sampling parameters, implemented in
#' \code{\link{durationFreq}}.
#'
#' Also see the conversion functions for sampling parameters at
#' \code{\link{SamplingConv}}.

#' @references
#' Solow, A. R., and W. Smith. 1997. On Fossil Preservation and the
#' Stratigraphic Ranges of Taxa. \emph{Paleobiology} 23(3):271-277.

#' @examples
#' #can simulate this type of data with sampleRanges
#'     # just set ranges.only = FALSE
#' #let's try a simulation example:
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' sampledOccurrences <- sampleRanges(taxa,r = 0.5,ranges.only = FALSE)
#' 
#' # now try with horizonSampRate
#' horizonSampRate(sampOcc = sampledOccurrences)
#' 
#' # but we could also try with the *other* inputs
#'    # useful because some datasets we may only have durations
#'    # and number of sampling events for
#' filtered <- sampledOccurrences[!is.na(sampledOccurrences)] 
#' dur <- sapply(filtered,max) - sapply(filtered,min)
#' nCol <- sapply(filtered,length)
#' # supply as durations and nCollections
#' horizonSampRate(durations = dur, nCollections=nCol)

#' @export
horizonSampRate<-function(sampOcc=NULL,durations=NULL,nCollections=NULL){
	#for estimating sampling rate from continuous time data
		#taken from Solow & Smith, 1997, Paleobiology
	#input is a list, with each element a taxon, consisting of n sampling occurrences
		#just like sampleRanges(data,ranges.only=FALSE)
	if(is.null(durations) & is.null(nCollections)){
		if(is.null(sampOcc)){
			stop("No input data supplied!")}
		if(!is.list(sampOcc)){
			stop("sampOcc is supplied but in't a list of species occurrences")}
		sampOcc<-sampOcc[!is.na(sampOcc)]
		nCollections<-sapply(sampOcc,length)
		durations<-sapply(sampOcc,max)-sapply(sampOcc,min)
		names(durations)<-names(nCollections)<-names(sampOcc)
		}
	#check lengths
	if(length(durations)!=length(nCollections)){
		stop("durations and nCollections are not the same length")
		}
	#check names
	if(is.null(names(durations)) | is.null(names(nCollections))){
		message("Input data lacks names, assuming ")
	}else{
		# check that names are the same, if not re-sort and re-check
		if(!identical(names(durations),names(nCollections))){
			message("Attempting to reorder durations and nCollections so names match")
			nCollections<-sapply(names(durations), USE.NAMES = FALSE, function(x){ 
				matches<-names(nCollections)==x
				if(sum(matches)==1){
					nCollections[matches]
				}else{
					stop("Name matches between durations and nCollections are not one-to-one")
					}
				})
			if(is.list(nCollections)){
				stop("Multiple matches of names on durations to names on nCollections")
				}
			if(length(durations)!=length(nCollections)){
				stop("durations and nCollections do not contain the same set of taxa, based on input names")
				}
			if(!identical(names(durations),names(nCollections))){
				stop("names still not matching after resorting")
				}
			}
		}
	if(is.null(durations) | is.null(nCollections)){
		stop("durations and nCollections have to be both supplied, if one is given")}
	sampRate<-(sum(nCollections-1)^2)/(sum(durations)*sum(nCollections))
	names(sampRate)<-"sampRate"
	return(sampRate)
	}