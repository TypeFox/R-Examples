#' @title Genetic Code Table
#' @description DNA Genetic Code Table
#' @param DNA if 'TRUE', DNA code will be used, else, 'RNA' code will be used. Default is 'TRUE'
#' @seealso \code{\link{codonToAAone}}
#' @seealso \code{\link{codonToAAthree}}
#' @seealso \code{\link{codonToAAoneRNA}}
#' @seealso \code{\link{codonToAAthreeRNA}}
#' @export
#' @examples
#' geneticCodeTable()
#' geneticCodeTable(DNA=FALSE)
geneticCodeTable <- function(DNA=TRUE){
	if(DNA){
		base=c('A','T','C','G')
		geneticCodes <-data.frame(FirstPosition=rep(base,each=16),SecondPosition=rep(base,each=4,times=4),ThirdPosition=rep(base,16))
		geneticCodes$GeneticCode=paste(geneticCodes$FirstPosition,geneticCodes$SecondPosition,geneticCodes$ThirdPosition,sep='')
		geneticCodes$AminoAcids=sapply(geneticCodes$GeneticCode,codonToAAthree)
		geneticCodes$AA=sapply(geneticCodes$GeneticCode,codonToAAone)
		return(geneticCodes)
	}else{
		base=c('A','U','C','G')
		geneticCodes <-data.frame(FirstPosition=rep(base,each=16),SecondPosition=rep(base,each=4,times=4),ThirdPosition=rep(base,16))
		geneticCodes$GeneticCode=paste(geneticCodes$FirstPosition,geneticCodes$SecondPosition,geneticCodes$ThirdPosition,sep='')
		geneticCodes$AminoAcids=sapply(geneticCodes$GeneticCode,codonToAAthreeRNA)
		geneticCodes$AA=sapply(geneticCodes$GeneticCode,codonToAAoneRNA)
		return(geneticCodes)
	}
}
