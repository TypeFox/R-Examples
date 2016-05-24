#' Converting Datasets of Taxon Ranges in Intervals Between timeList format and fourDate format
#'
#' Functions for manipulating data where the first and last appearances of taxa
#' are known from bounded intervals of time. The two main functions listed here
#' are for converting between (1)
#' a data structure consisting of a single 'flat' table where each taxon is listed as a
#' set of four dates (a \code{fourDate} data type), and (2) a list format where each taxon
#' is listed as its first and last intervals, with an associated table of age bounds for
#' the intervals referred to in the first table (referred to as a \code{timeList} data
#' structure by many \code{paleotree} functions).

#' @details
#' \code{timeList2fourDate} is for converting from a \code{timeList} format to
#' a \code{fourDate} format. \code{fourDate2timeList} is for converting from
#' a \code{fourDate} format to a \code{timeList} format.

#' @param timeList A list composed of two matrices with two columns each, the first
#' giving interval start and end date bounds, and the second giving taxon first and
#' last interval appearances in reference to the intervals listed in the first matrix.

#' @param fourDate A four column matrix where each row is a different taxon, the first two columns
#' are the lower and upper bounds on the time of first appearance for that taxon and the third
#' and fourth columns are respectively the lower and upper bounds on the time of last appearance
#' for that taxon, all in time before present.

#' @return
#' A converted data object, respective to the function applied.

#' @aliases fourDateFunctions

#' @seealso
#' \code{\link{bin_timePaleoPhy}} and \code{\link{taxicDivDisc}} for common applications;
#' \code{\link{binTimeData}} for a simulation function for such data objects

#' @author David W. Bapst

#' @references
#' See my recent blog post on temporal datasets in paleontology for some details:
#' 
#' http://nemagraptus.blogspot.com/2015/02/how-do-we-treat-fossil-age-data-dates.html

#' @examples
#' #timeList object from the retiolinae dataset
#' data(retiolitinae)
#'
#' str(retioRanges)
#'
#' taxicDivDisc(retioRanges)
#' 
#' fourDateRet<-timeList2fourDate(retioRanges)
#' 
#' # total uncertainty in retio first and last appearances?
#' sum((fourDateRet[,1]-fourDateRet[,2])+(fourDateRet[,3]-fourDateRet[,4]))
#'
#' #convert back
#' newTimeList<-fourDate2timeList(fourDateRet)
#' taxicDivDisc(retioRanges)
#'

#' @name timeList2fourDate
#' @rdname timeList2fourDate
#' @export
timeList2fourDate<-function(timeList){
	if(!is(timeList,"list")){stop("timeList is not a list?")}
	if(!(length(timeList)==2)){stop("timeList is not length=2")}
	timeList<-lapply(timeList,as.matrix)
	if(!all(sapply(timeList,inherits,what="matrix"))){
		stop("timeList elements are not matrices")
		}
	if(!all(sapply(timeList,function(x) ncol(x)==2))){
		stop("timeList matrices do not have two columns")
		}
	firstInt<-t(sapply(timeList[[2]][,1],function(x) timeList[[1]][x,]))
	lastInt<-t(sapply(timeList[[2]][,2],function(x) timeList[[1]][x,]))
	fourDate<-cbind(firstInt,lastInt)
	return(fourDate)
	}

#' @name fourDate2timeList
#' @rdname timeList2fourDate
#' @export
fourDate2timeList<-function(fourDate){
	fourDate<-as.matrix(fourDate)
	if(!is(fourDate,"matrix")){if(ncol(fourDate)!=4){
		stop("fourDate must be a matrix with four columns")}}
	taxaFirst<-fourDate[,1:2,drop=FALSE]
	taxaLast<-fourDate[,3:4,drop=FALSE]
	#make interval list
	intTimes<-unique(rbind(taxaFirst,taxaLast))
	intTimes<-intTimes[order(-intTimes[,1],-intTimes[,2]),]
	#now assign taxa first and last intervals
	firstInt<-apply(taxaFirst,1,function(x)
		which(apply(intTimes,1,function(y) identical(y,x))))
	lastInt<-apply(taxaLast,1,function(x)
		which(apply(intTimes,1,function(y) identical(y,x))))
	taxonTimes<-cbind(firstInt,lastInt)
	#package it together
	dimnames(intTimes)<-list(NULL,c("startTime","endTime"))
	dimnames(taxonTimes)<-list(NULL,c("firstInt","lastInt"))
	res<-list(intTimes=intTimes,taxonTimes=taxonTimes)
	return(res)
	}