#' Ordinates interaction matrix
#'
#' 'OrderMatrix' ordinates an interaction matrix scores derived from reciprocal
#' averaging (Gauch et al. 1977). These scores represent a latent environmental
#' gradient along which species distributions are structured.
#'
#'
#' @param comm community data in the form of a presence absence matrix
#' @param scores axis scores to ordinate matrix. 1: primary axis scores
#' (default) 2: secondary axis scores
#' @param outputScores logical. Default is to return the ordinated matrix. If
#' 'outputScores' is TRUE, the function returns the site and species scores.
#' @param binary logical argument indicating whether to ordinate the community
#' matrix based on abundance or binary (default) data.
#' @return 'OrderMatrix' returns either an ordinated matrix (outputScores =
#' FALSE) or the site and species scores (outputScores = TRUE) obtained from
#' reciprocal averaging. This function is already contained within functions
#' calculating coherence, species turnover & boundary clumping, but may be
#' useful for visualizations or for hypothesis testing concerning the important
#' variables associated with the site or species scores.
#' @note 'OrderMatrix', like many of these functions, relies heavily on the
#' 'vegan' package.
#' @author Tad Dallas
#' @export
#' @references Gauch, H.G., R.H. Whittaker, and T.R. Wentworth. 1977. A
#' comparative study of reciprocal averaging and other ordination techniques.
#' Journal of Ecology 65:157-174.
#'
#' Leibold, M.A. and G.M. Mikkelson. 2002. Coherence, species turnover, and
#' boundary clumping: elements of meta-community structure. Oikos 97: 237 -
#' 250.
#'
#' Oksanen,J., F.G. Blanchet, R. Kindt, P. Legendre, P.R. Minchin, R.B. O'Hara,
#' G.L. Simpson, P. Solymos, M.H.H. Stevens and H. Wagner (2012). vegan:
#' Community Ecology Package. R package version 2.0-4.
#' http://CRAN.R-project.org/package=vegan
#' @keywords ordination
#' @examples
#'
#' #define an interaction matrix
#' data(TestMatrices)
#' pres3c <- TestMatrices[[6]]
#'
#' #obtain an ordinated interaction matrix
#' OrderMatrix(pres3c, scores = 1, outputScores = FALSE)
#'
#' #obtain site and species scores from reciprocal averaging
#' OrderMatrix(pres3c, scores = 1, outputScores = TRUE)
#'
#'
OrderMatrix <-function(comm, scores=1, outputScores=FALSE, binary=TRUE){
	if(binary==TRUE){comm=(comm>0)+0}
	#reciprocal averaging
	temp <- decorana(comm,ira=1)
	#ordering of matrix
	if(outputScores == FALSE){
  	ret <- comm[order(temp$rproj[,scores], decreasing=FALSE),
                order(temp$cproj[,scores],decreasing=FALSE)]
	ret <- as.matrix(ret)
  ret <- apply(ret,2,rev)
	ret <- 0+(ret>0)
	return(ret)}

	if(outputScores == TRUE){
		return(list(speciesscores=temp$cproj[,scores],
                sitescores=temp$rproj[,scores]))}
}
