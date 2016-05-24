#' Plotting Occurrence Data Across Taxa
#'
#' \code{plotOccData} takes occurrence data which has been sorted into a by-taxon list,
#' such as that output by \code{taxonSortPBDBocc} or may be output by simulations using
#' \code{sampleRanges} and produces a plot showing the age uncertainty associated with
#' individual occurrences, with occurrences of the same taxon grouped by color.

#' @details
#' This function was originally conceived of in the following blog post: 
#' \href{http://nemagraptus.blogspot.com/2015/02/how-do-we-treat-fossil-age-data-dates.html}{Link}
#'

#' @param occList A list where every element is a table of occurrence data for a different taxon,
#' such as that returned by \code{\link{taxonSortPBDBocc}}. The occurrence data can be either a 
#' two-column matrix composed of the lower and upper age bounds on each taxon occurrence, or has
#' two named variables which match any of the field names given by the PBDB API under either
#' the 'pbdb' vocab or 'com' (compact) vocab for early and late age bounds.

#' @param groupLabel A character vector with a single string giving the name for
#' the occurrence dataset used, such as the taxonomic name of the group examined.
#' If not given (the default) a generic plot title is appended.

#' @param occColors A vector of numbers or characters indicating colors on a color
#' palette for use with basic plot. Must be the same length as occList. If empty, the
#' default, the colors used are sampled randomly from the \code{rainbow()} function.

#' @param lineWidth A numeric value giving the length to be used for the width of lines
#' plotted in \code{plotOccData}. If not given (the default), this is calculated using
#' an algorithm that selects an optimal line width for plotting.

#' @param xlims A two element vector controlling the width of the horizontal time-scale
#' the occurrence bars are plotted against. By default, this is not given and calculated
#' internally.

#' @return
#' This function will invisibly return a list, with each per-taxon element containing the
#' two-column matrix of age bounds for occurrences.

#' @seealso
#' \code{\link{taxonSortPBDBocc}}, \code{\link{occData2timeList}} 
#' and the example graptolite dataset at \code{\link{graptPBDB}}

#' @author David W. Bapst

#' @examples
#' #load example graptolite PBDB occ dataset
#' data(graptPBDB)
#' 
#' #get formal genera
#' occSpecies<-taxonSortPBDBocc(graptOccPBDB, rank="species")
#' 
#' #plot it!
#' plotOccData(occSpecies)
#' 
#' #this isn't too many occurrences, because there are so few
#'     #formal grapt species in the PBDB
#'
#' #genera is messier...
#' 
#' #get formal genera
#' occGenus<-taxonSortPBDBocc(graptOccPBDB, rank="genus")
#' 
#' #plot it!
#' plotOccData(occGenus)
#' 
#' #some of those genera have occurrences with very large
#'    #age uncertainties on them!
#' 

#' @name plotOccData
#' @rdname plotOccData
#' @export
plotOccData<-function(occList,groupLabel=NULL,occColors=NULL,lineWidth=NULL,xlims=NULL){
	#check groupLabel
	if(!is.null(groupLabel)){if(length(groupLabel)!=1){
		stop("groupLabel must be length=1")}}
	#check occColors
	if(!is.null(occColors)){if(length(occColors)!=length(occList)){
		stop("occColors must be the same length as occList")}}
	
	#use pullOccListData to get occurrence data
	occList<-pullOccListData(occList)
	#
	#number of occurrences
	sumOcc<-sum(sapply(occList,nrow))
	#
	#order taxa by earliest occurrence first, then later occurrence
	occList<-occList[order(sapply(occList,max))]
	#order the occurrences within taxa
	occList<-lapply(occList,function(x) x[order(x[,1],x[,1]),, drop=FALSE])
	#
	#set xlims
	if(is.null(xlims)){
		xlims<-c(max(sapply(occList,max)), min(sapply(occList,min))) 
		xlimMod<-(xlims[1]-xlims[2])*0.01
		xlims<-c(xlims[1]+xlimMod,xlims[2]-xlimMod)
		}
	#plot title
	plotMainTitle<-ifelse(is.null(groupLabel),"Age Uncertainty of Occurrence Data",
             paste("Age Uncertainty of Occurrence Data for",groupLabel))
	#initiate the plot with a modifiable main title
	origPar<-par(no.readonly=TRUE)	#save original parameters	
	par(yaxt="n")
	plot(0, 0, type="n", xlim=xlims,
		ylim=c(0,sum(sapply(occList,nrow))+(2*length(occList))),
		ylab = "", xlab = "Time (Ma)",
		main=plotMainTitle)
	#set colors
	if(is.null(occColors)){
		occColors<-sample(rainbow(length(occList)))  #scramble colors
		}
	#set line width
	if(is.null(lineWidth)){
		lineWidth<-(-0.02*nrow(sumOcc))+3.2
		lineWidth<-ifelse(lineWidth>0,lineWidth,0.01)
		}
	#now plot the occurrences as lines
	count<-1
	for(i in 1:length(occList)){
		for(j in 1:nrow(occList[[i]])){
			lines(occList[[i]][j,],c(count,count),lwd=lineWidth,col=occColors[i])
			count<-count+1
			}
		count<-count+2
		}
	#reset graph par
	par(origPar)
	#return the occList as an invisible data object
	return(invisible(occList))
	}