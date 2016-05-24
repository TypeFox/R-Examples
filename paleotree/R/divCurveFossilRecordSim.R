#' Diversity-Curve Plotting for Simulations of Diversification and Sampling In the Fossil Record
#'
#' An extremely simple plotting function, which plots the original taxonomic diversity
#' versus the sampled taxonomic diversity, for use with output from the function \code{simFossilRecord}.
#' If sampling processes were not included in the model, then it plots simply the
#' single diversity curve.

#' @details
#' This function is essentially a wrapper for \code{paleotree} function \code{multiDiv}.

#' @param fossilRecord A list object output by \code{simFossilRecord}, often composed
#' of multiple elements, each of which is data for 'one taxon', with the first
#' element being a distinctive six-element vector composed of numbers, corresponding
#' to the six variable tables by \code{fossilRecord2fossilTaxa} after simulating with
#' \code{simFossilRecord} (originally produced by deprecated function \code{simFossilTaxa}).

#' @param merge.cryptic If \code{TRUE}, cryptic taxon-units (i.e.
#' those in the same cryptic complex) will be merged into single taxa for the
#' sake of being counted in the diversity curves presented by this function.

#' @param plotLegend A logical. Should a legend be plotted? Only applies if sampling
#' processes were modeled.

#' @param legendPosition Where should the legend be plotted? See help for \code{legend}
#' for details. Only applies if sampling processes were modeled.

#' @param curveColors A vector of length two indicating what colors the original and sampled
#' diversity curves should be displayed in. Only applies if sampling processes were modeled.	

#' @param curveLineTypes A vector of length two indicating what colors the original and sampled
#' diversity curves should be displayed in. Only applies if sampling processes were modeled.

#' @return
#' This function returns nothing: it just creates a plot.

#' @seealso
#' \code{\link{simFossilRecord}}

#' @author David W. Bapst

#' @examples
#'
#' set.seed(44)
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1,
#' 	nTotalTaxa=c(20,30) ,nExtant=0, plot=FALSE)
#' 
#' # now let's plot it
#' divCurveFossilRecordSim(record)
#'

#' @name divCurveFossilRecordSim
#' @rdname divCurveFossilRecordSim
#' @export		
divCurveFossilRecordSim<-function(fossilRecord,merge.cryptic=TRUE,plotLegend=TRUE,
	legendPosition="topleft",curveColors=c("black","red"),curveLineTypes=c(1,2)){
	#
	taxaConvert<-fossilRecord2fossilTaxa(fossilRecord=fossilRecord)
	#taxicDivCont(taxaConvert,int.length=0.2)
	#are any sampled AT NOT THE PRESENT?
	areSampled<-whichSampledInPast(fossilRecord)
	if(length(areSampled)>0){
		fossilRanges<-fossilRecord2fossilRanges(fossilRecord=fossilRecord, 
			merge.cryptic=merge.cryptic, ranges.only = TRUE)
		curveList<-list(taxaConvert,fossilRanges)
		#if(i==5 & nruns==1000){browser()}
		multiDiv(curveList,plotMultCurves=TRUE,drop.cryptic=merge.cryptic,
			divPalette=curveColors,divLineType=curveLineTypes,main="")
		if(plotLegend){
			legend(legendPosition,legend=c("True Richness", "Sampled Richness"),
				col=curveColors,lty=curveLineTypes)
			}
	}else{ #none sampled
		taxicDivCont(taxaConvert,int.length=0.2)
		}
	}	
	
whichSampledInPast<-function(taxa){
	#which taxa are sampled prior to the modern?
	areSamp<-sapply(taxa,function(x) length(x[[2]])>0)
	#are there any extant taxa
	areLive<-which(sapply(taxa,function(x) x[[1]][5]==1))
	nExtant<-length(areLive)
	if(sum(areSamp)>0 & nExtant>0){
		#identify modern time
			# assuming that latest sampling time of extant taxa is modernTime
		possModernTime<-min(unlist(sapply(taxa[areLive],function(x) x[[2]])))
		sampOnce<-sapply(taxa,function(x) length(x[[2]])==1)
		sampModernOnly<-sapply(taxa,function(x)
			if(length(x[[2]]==1)){
				any(x[[2]]==possModernTime)
			}else{
				FALSE
				})
		res<-which(areSamp & !sampModernOnly)
	}else{
		res<-which(areSamp)
		}
	return(res)
	}
