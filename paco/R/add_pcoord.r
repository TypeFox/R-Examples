#' Principle coordinate analysis of phylogenies
#' @param D A list with objects H, P, and HP, returned by prepare_paco_data
#' @param correction Choose the method with which to correct negative eigenvalues ('none', cailliez', 'lingoes'). Default is 'none'.
#' @return The input list with added objects for the principal coordinates of the objects
#' @note Internal function coordpcoa is a modified version of ape::pcoa, utilising vegan::eigenvals
#' @export
#' @examples
#' data(gopherlice)
#' library(ape)
#' gdist <- cophenetic(gophertree)
#' ldist <- cophenetic(licetree)
#' D <- prepare_paco_data(gdist, ldist, gl_links)
#' D <- add_pcoord(D)
add_pcoord <- function(D, correction='none')
{
HP_bin <- which(D$HP >0, arr.ind=TRUE)
H_PCo <- coordpcoa(D$H, correction)$vectors #Performs PCo of Host distances
P_PCo <- coordpcoa(D$P, correction)$vectors #Performs PCo of Parasite distances
D$H_PCo <- H_PCo[HP_bin[,1],] #Adjust Host PCo vectors
D$P_PCo <- P_PCo[HP_bin[,2],] #Adjust Parasite PCo vectors
return(D)
}
