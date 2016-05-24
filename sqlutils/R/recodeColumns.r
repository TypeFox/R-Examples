#' Recodes factors with more than \code{maxLevels} to characters.
#' @param df the data frame to recode.
#' @param maxLevels the maximum number of levels a factor can have before being
#'        converted to a character.
recodeColumns <- function(df, maxLevels=20) {
	for(c in seq_len(ncol(df))) {
		if(class(df[,c])[1] == 'factor' & length(levels(df[,c])) > maxLevels) {
			df[,c] = as.character(df[,c])
		}
	}
	return(df)	
}
