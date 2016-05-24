#' Function to detect whether a character is continuous or discrete
#'
#' @param dat A vector of data
#' @param repeatsAsDiscrete If TRUE, consider numeric variables that repeat values exactly as discrete; see cutoff
#' @param cutoff Cutoff value for deciding if numeric data might actually be discrete: if nlev is the number of levels and n the length of dat, then nlev / n should exceed cutoff, or the data will be classified as discrete
#' @return Either "discrete" or "continuous"
#' @examples
#' data(anolis)
#' detectCharacterType(anolis$dat[,1])
#' @export
detectCharacterType<-function(dat, repeatsAsDiscrete=TRUE, cutoff=0.1) {
	if(is.factor(dat)) {
			charType<-"discrete"
	} else if(nlevels(as.factor(dat))/length(dat) < cutoff) {
			warning("Guessing that this is a discrete character based on repeated values")
			charType<-"discrete"
	} else {
			charType<-"continuous"
	}
	return(charType)
}		# needless to say, this is not yet robust

#' Apply detectCharacterType over an entire matrix
#'
#' @param mat A matrix of data
#' @param repeatsAsDiscrete If TRUE, consider numeric variables that repeat values exactly as discrete; see cutoff
#' @param cutoff Cutoff value for deciding if numeric data might actually be descrete: if nlev is the number of levels and n the length of dat, then nlev / n should exceed cutoff, or the data will be classified as discrete
#' @return Vector of either "discrete" or "continuous" for each variable in matrix
#' @examples
#' data(anolis)
#' detectAllCharacters(anolis$dat)
#' @export
detectAllCharacters<-function(mat, repeatsAsDiscrete=TRUE, cutoff=0.1) {
  mat <- as.matrix(mat)
	nchar<-dim(mat)[2]
	result<-numeric(nchar)
	for(i in 1:nchar) {
		result[i]<-detectCharacterType(mat[,i], repeatsAsDiscrete, cutoff)
	}
	return(result)
}

#' Filter a matrix, returning either all continuous or all discrete characters
#'
#' @param mat A matrix of data
#' @param charType A vector of character types (perhaps from detectAllCharacters)
#' @param returnType Either discrete or continuous
#' @return Matrix with only discrete or continuous characters
#' @examples
#' data(anolis)
#' aType<-detectAllCharacters(anolis$dat)
#' filterMatrix(anolis$dat, aType, "discrete")
#' @export

filterMatrix<-function(mat, charType, returnType="discrete") {
  rType<-match.arg(returnType, c("discrete", "continuous"))
  columnFilter<-charType==rType
  result<-mat[,columnFilter]
  return(result)
}

#' Row and column name check
#'
#' @param dat A vector of data
#' @param nameType, either:
#' \describe{
#'     \item{"row"}{Rows}
#'   	 \item{"col"}{Columns}
#' 		 \item{"rowcol"}{Both rows and columns}
#'	}
#' @examples
#' data(anolis)
#' hasNames(anolis$dat, "row")
#' @export
hasNames <- function(dat, nameType="row") {
	nType = match.arg(nameType, c("row", "col", "rowcol"))
	if(nType=="row") {
		res<-!is.null(rownames(dat))
	}
	if(nType=="col") {
		res<-!is.null(colnames(dat))
	}
	if(nType=="rowcol") {
		res<-!is.null(rownames(dat)) & !is.null(colnames(dat))
	}
	names(res)<-nameType
	res
}

#' Force names for rows, columns, or both
#'
#' @param dat A vector of data
#' @param nameType, either:
#' \describe{
#'     \item{"row"}{Rows}
#' 		 \item{"col"}{Columns}
#' 		 \item{"rowcol"}{Both rows and columns}
#'	}
#' @examples
#' data(anolis)
#' forceNames(anolis$dat, "row")
#' @export
forceNames <- function(dat, nameType="row") {
	nType = match.arg(nameType, c("row", "col", "rowcol"))
	if(nType=="row" | nType=="rowcol") {
		if(!hasNames(dat, nameType="row")) {
			nrows<-dim(dat)[1]
			rownames(dat) <- paste("n", 1:nrows, sep="")
		}
	}
	if(nType=="col" | nType=="rowcol") {
		if(!hasNames(dat, nameType="col")) {
			ncols<-dim(dat)[2]
			colnames(dat) <- paste("n", 1:ncols, sep="")
		}
	}

	dat
}
