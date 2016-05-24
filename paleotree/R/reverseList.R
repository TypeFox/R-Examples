#' Reverse List Structure
#' 
#' Takes a list and reverses the list structure, such that list composed of
#' five elements with eight sub-elements is restructured to have eight elements
#' with five sub-elements each, with the order of elements and sub-elements
#' being retained despite their reversal in hierarchical position.
#' 
#' @details 
#' The function will fail and return an error if all sub-elements are not
#' vectors or lists of equal length.
#' 
#' This function can be useful for instances when each element of a list is
#' by-sample, composed of multiple, different tests on that sample, but where
#' for further analysis/plotting, it would be beneficial to have a list where
#' each element represented values from the same test performed across multiple
#' samples (i.e. plotting a box-plot).
#' 
#' @param list A list composed of multiple elements, with each element a vector
#' or list of equal length

#' @param simplify Should the result be simplified, as the argument in sapply
#' @return Returns a list with a reversed structure relative to the input, see
#' above.

#' @author David W. Bapst
#' @examples
#' 
#' list1<-list(list(1:3),list(1:3),list(1:3))
#' reverseList(list1,simplify=FALSE)
#' reverseList(list1,simplify=TRUE)
#' 
#' @export reverseList
reverseList<-function(list,simplify=FALSE){
	#reverses primary and secondary list structure
		#i.e. if a list is 10 elements each 50 long, get 50 elements 10 long
	if(length(unique(sapply(list,length)))!=1){
		stop("Not all lists equally long")}
	list1<-list()
	for(i in 1:length(list[[1]])){
		list1[[i]]<-sapply(list,function(x) x[[i]],simplify=simplify)
		}
	names(list1)<-names(list[[1]])
	return(list1)
	}
