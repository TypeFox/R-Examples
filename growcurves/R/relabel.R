#' Relabel user vector input to sequential numerical 
#'
#' Takes user input for vector in either character or numerical format and converts it to sequential numeric format.
#' For use in functions requiring sequential numerical format.  Returns new sequential numerical vector and 
#' vector of unique values inputted by user.  The latter is used to label plot variables returned to the user.
#'
#' @param label.input  A vector in numerical or character format that the user desires to convert to sequential numeric.
#' @param start An integer representing the starting value of the sequential sequence for the new label vector.
#' @return A list object containing the new sequential vector, \code{label.new}, and a vector of unique input values,
#'  		\code{labeli.u}.
#' @seealso \code{\link{dpgrowmm}}, \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases relabel
#' @export
relabel = function(label.input, start)
{

  ## check input
  ## if( is.factor(label.input) ) stop("Input vector must be either as.numeric or as.character")
  start		= as.integer(start)

  ## extract unique entries from user input label
  labeli.u			= unique(label.input) 
 
  ## if label.input is numeric, we align order of the new label. label.input bigger -> label.new bigger
  ## the lowest value of label.input is chosen as baseline
  ## otherwise (for character input), the first value of the new label is assumed as baseline - this is KEY
  if( is.numeric(labeli.u) ) 
  {
	labeli.u		= sort(labeli.u)
  }

  ## create new label, label.new, from start position input by user in a sequential fashion
  end				= start + (length(labeli.u) - 1)
  labeln.u			= start:end ## create sequential 0, 1, 2, ..
  dat.label			= data.frame(labeln.u,labeli.u, stringsAsFactors = FALSE)
  names(dat.label)		= c("label.new","label.input")

  ## merge label.input with data.frame map of unique values to reinflate to the label.new
  label.idat			= data.frame(1:length(label.input),label.input, stringsAsFactors = FALSE)
  names(label.idat)		= c("index","label.input")

  ## sort = FALSE produces the label.new in the order of the user input, label.input case vector, which is what we want
  dat.label			= merge(label.idat, dat.label, by = "label.input", all.x = T, sort = FALSE)
  dat.label			= dat.label[order(dat.label$index),] ## ensures order of dat.label is by label.input
  rownames(dat.label)		<- NULL
  dat.label$index		<- NULL
  label.new			= dat.label$label.new ## now in case format
  labeli.u			= as.matrix(labeli.u)

  ## label.new is the new label in the same length / case pattern as label.input
  ## labeli.u is the vector of unique labels extracted from the user input
  out	= list(label.new = label.new, labeli.u = labeli.u, dat.label = dat.label)
  return(out)

} ## end function relabel