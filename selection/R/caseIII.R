#' Corrects correlations using Case III
#' 
#' Using Thorndike's Case III correction, \code{caseIII} corrects the xy
#' correlation for direct restriction on z (and by implication, indirect
#' restriction on x)
#' 
#' The Case III correction is defined as follows insert later The result is an
#' unbiased estimate of the unattenuated correlation between X and Y
#' 
#' @param data a dataset containing the two incidentally restricted variables (X and Y) and complete
#' information on the selection variable (Z). Conversely, one could supply values for rxy, rzy, and rxz
#' instead. 
#' @param x The column index (or name) of the X variable 
#' @param y The column index (or name) of the Y variable 
#' @param z The column index (or name) of the Z variable (the one used for selection)
#' @param rxy the restricted correlation between x (the indirectly selected
#' variable) and y (the outcome variable).
#' @param rzy the restricted correlation between z (the selection variable) and
#' y (the outcome variable).
#' @param rxz the restricted correlation between x (the indirectly selected
#' variable) and z (the selection variable).
#' @param uz the ratio of restricted to unrestricted variance (i.e.,
#' sigmaz'/sigmaz).
#' @return a scalar that is the unbiased estimate of the correlation between X
#' and Y.
#' @author Dustin Fife
#' @export
#' @seealso \code{\link{caseIV}}, \code{\link{caseIIIR}}, \code{\link{em}},
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
#' cor.mat = cor(new.dat[,c("R", "Biodata", "Performance")])
#' 	# correct assuming direct selection on R, indirect on biodata, and a dv of performance
#' corrected = caseIII(rxy=cor.mat[1,3], rzy=cor.mat[2,3], 
##' 		rxz=cor.mat[1,2], uz = sd(new.dat$R)/sd(selection.example.data$R))	
#' corrected
#' ## do a simulation to show it works
#' cor.mat = matrix(c(1, .3, .4, 
#' 					.3, 1, .5,
#' 					.4, .5, 1), nrow=3)
#' data = mvrnorm(100000, mu=c(0,0,0), Sigma = cor.mat)					
#' ### restrict the data
#' data[data[,1]<.5, 2:3] = NA
#' caseIII(data=data, x=2, y=3, z=1)
caseIII = function(data=NULL, x=2, y=3, z=1, rxy, rzy, rxz, uz){
	if (!is.null(data)){
		uz = sd(data[,z])/sd(data[!is.na(data[,y]),z], na.rm=T)
		rxy = cor(data[,x], data[,y], use="pairwise.complete.obs"); rxz = cor(data[,x], data[,z], use="pairwise.complete.obs"); rzy = cor(data[,z], data[,y], use="pairwise.complete.obs")
	} else {
		uz = 1/uz
	}
		num = rxy-rxz*rzy + rxz*rzy*uz^2
		denom = sqrt((1-rxz^2 + rxz^2*uz^2)*(1-rzy^2 + rzy^2*uz^2))
		return(num/denom)
}