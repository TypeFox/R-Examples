#' Corrects correlations using Case III-R
#' 
#' Using Fife, Hunter, and Mendoza's (2013) Case III-R correction,
#' \code{caseIIIR} corrects the tp correlation for direct restriction on x (and
#' by implication, indirect restriction on x)
#' 
#' later
#' 
#' @param data.matrix a data matrix with rows representing observations and
#' columns representing variables. Note: no observations should be present when
#' using case III-R.
#' @param rtp the restricted correlation between t (the predictor variable) and
#' p (the outcome variable).
#' @param rsp the restricted correlation between s (the suitability variable)
#' and p (the outcome variable).
#' @param rts the restricted correlation between t (the predictor variable) and
#' s (the suitability variable).
#' @param SR the selection ratio (i.e., the proportion of applicants who were
#' selected)
#' @return a scalar that is the unbiased estimate of the correlation between T
#' and P.
#' @author Dustin Fife
#' @export
#' @seealso \code{\link{caseIV}}, \code{\link{caseIII}}, \code{\link{em}},
#' \code{\link{rel.correction}}
#' @references Thorndike, R. L. (1949). Personnel selection: Test and
#' measurement techniques. Oxford, England: Wiley.
#' 
#' Pearson, K. (1903). Mathematical contributions to the theory of evolution.
#' XI. On the influence of natural selection on the variability and correlation
#' of organs. Philosophical Transactions of the Royal Society of London. Series
#' A, Containing Papers of a Mathematical or Physical Character, 200, 1-66.
#' @examples
#' 
#' 	# load example data
#' data(selection.example.data)
#' 	# give me only those rows that have full information
#' new.dat = selection.example.data[!is.na(selection.example.data$Performance),]
#' 
#' 	# correct assuming direct selection on R, indirect on biodata, and a dv of performance
#' corrected = caseIIIR(data.matrix = new.dat, SR=.5)	
#' 
#' corrected = caseIIIR(rtp=.4, rsp=.3, rts=.5, SR=.5)	
#' 
caseIIIR = function(data.matrix=NULL, rtp, rsp, rts, SR){

		#### if a data matrix is supplied, compute the correlations
	if (!is.null(data.matrix)){
		cor.mat = cor(data.matrix, use="complete.obs")
		rtp = cor.mat[2,3]; rsp = cor.mat[1,3]; rts = cor.mat[1,2]
	}
		##### estimate the selection ratio by integrating
	zval = qnorm(SR)
	ud = normalVar(zval)
				
		##### compute case III				
	c.III.R = caseIII(rxy=rtp, rzy=rsp, rxz=rts, uz=ud)
	return(c.III.R)			
}
