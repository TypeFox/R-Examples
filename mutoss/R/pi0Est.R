# Collection of elementary functions calculating an estimate
# of the number m0 resp. the proportion pi0 of true null hypotheses 
# in a finite family of hypotheses. 
# 
# Author: MarselScheer and WiebkeWerft
###############################################################################



storey_pi0_est = function(pValues, lambda) 
{
	pi0 = (sum(pValues > lambda) + 1) / (1 - lambda) / length(pValues)	
	return(list(pi0 = pi0, lambda = lambda))
}

mutoss.storey_pi0_est <- function() { return(new(Class="MutossMethod",
					label="Storey-Taylor-Siegmund (2004) Procedure",
					callFunction="storey_pi0_est",
					output=c("pi0", "lambda"),
					info="<h2>Storey-Taylor-Siegmund procedure</h2>\n\n\
							<p>The Storey-Taylor-Siegmund procedure for estimating pi0 is applied to pValues.\
							The formula is equivalent to that in Schweder and Spjotvoll (1982),\
							page 497, except the additional '+1' in the nominator that\
							introduces a conservative bias which is proven to be sufficiently large\
							for FDR control in finite families of hypotheses if the estimation\
							is used for adjusting the nominal level of a linear step-up test.</p>\n\
							<h3>Reference:</h3>\
							<ul>\
							<li>Storey, J. D., Taylor, J. E. and Siegmund, D. (2004). \"<i> Plots of P-values to evaluate many tests simultaneously. </i>\" Biometrika 69, 3, 493-502. </li>\n\
							<li>Huang, Y. and Hsu, J. (2007). \"<i> Strong control, conservative point estimation and simultaneous conservative consistency of false discovery rates: a unified approach.</i>\" JRSS B 66, 1, 187-205.</li>
							</ul>",
					parameters=list(pValues=list(type="numeric"), lambda=list(type="numeric"))
			)) }


ABH_pi0_est <- function(pValues)
{
	m <- length(pValues)
	index <- order(pValues)
	spval <- pValues[index]
	m0.m <- rep(0, m)
	for (k in 1:m) {
		m0.m[k] <- (m + 1 - k)/(1 - spval[k])
	}
	idx <- which(diff(m0.m, na.rm = TRUE) > 0)
	if (length(idx) == 0)
		grab <- 2
	else grab <- min(idx, na.rm = TRUE)+1
	pi0.ABH <- (ceiling(min(m0.m[grab], m))/m)
	return(list(pi0 = pi0.ABH))
}



mutoss.ABH_pi0_est <- function() { return(new(Class="MutossMethod",
					label="Hochberg-Benjamini (1990) lowest slope line method",
					callFunction="ABH_pi0_est",
					output=c("pi0"),
					info="<h2>Hochberg-Benjamini lowest slope line method</h2>\n\n\
							<p>The Lowest Slope Line (LSL) method of Hochberg and Benjamini for estimating pi0 is applied to pValues.\
							This method for estimating pi0 is motivated by the graphical approach proposed\
							by Schweder and Spjotvoll (1982), as developed and presented in Hochberg and Benjamini (1990).</p>\n\
							<h3>Reference:</h3>\
							<ul>\
							<li>Hochberg, Y. and Benjamini, Y. (1990). \"<i> More powerful procedures for multiple significance testing. </i>\" Statistics in Medicine 9, 811-818. </li>\n\
							<li>Schweder, T. and Spjotvoll, E. (1982). \"<i> Plots of P-values to evaluate many tests simultaneously.</i>\" Biometrika 69, 3, 493-502.</li>
							</ul>",
					parameters=list(pValues=list(type="numeric"))
			)) }


TSBKY_pi0_est <- function(pValues, alpha)
{	
	m <- length(pValues)
	adjp <- p.adjust(pValues,"BH")
	pi0.TSBKY <- ((m - sum(adjp < alpha/(1 + alpha), na.rm = TRUE)) / m)
	return(list(pi0=pi0.TSBKY))
}


mutoss.TSBKY_pi0_est <- function() { return(new(Class="MutossMethod",
					label="Benjamini, Krieger and Yekutieli (2006) two-step estimation method",
					callFunction="TSBKY_pi0_est",
					output=c("pi0"),
					info="<h2>Hochberg-Benjamini lowest slope line method</h2>\n\n\
							<p>The two-step estimation method of Benjamini, Krieger and Yekutieli for estimating pi0 is applied to pValues.
							It consists of the following two steps:</p>\
							<p>Step 1. Use the linear step-up procedure at level alpha' =alpha/(1+alpha). Let r1 be the number of\
							rejected hypotheses. If r1=0 do not reject any hypothesis and stop; if r1=m reject all m\
							hypotheses and stop; otherwise continue.</p>\n\
							<h3>Reference:</h3>\
							<ul>\
							<li>Benjamini, Y., Krieger, A. and Yekutieli, D. (2006). \"<i> Adaptive linear step-up procedures that control the false discovery rate. </i>\" Biometrika 93, 3, page 495. </li>\n\
							</ul>",
					parameters=list(pValues=list(type="numeric"),alpha=list(type="numeric"))
			)) }
# TODO: GB (MS) This will help will look terrible!
BR_pi0_est <- function(pValues, alpha, lambda=1, truncate = TRUE)
{	
	
	if ( lambda <= 0 || lambda >= 1/alpha) {
		stop('BR_pi0_est() : lambda should belong to (0, 1/alpha)')
	}
	
	m <- length(pValues)
	stage1 <- indepBR( pValues, alpha, lambda, silent = TRUE)
	pi0 <- ( m + 1 - sum(stage1$rejected) ) / ( m * ( 1 - lambda*alpha ) )
	
	if (truncate) {
		pi0 = min(1,pi0)
	}
	
	return(list(pi0=pi0))
}


mutoss.BR_pi0_est <- function() { return(new(Class="MutossMethod",
					label="Blanchard-Roquain (2009) estimation method",
					callFunction="BR_pi0_est",
					output=c("pi0"),
					info="<h2> Blanchard-Roquain estimation under independence </h2>\n\n\
							<h3>Reference:</h3>\
							<ul>\
							<li>Blanchard, G. and Roquain, E. \"<i> Adaptive False Discovery Rate Control under Independence and Dependence.</i>\" 
							Journal of Machine Learning Research 10:2837-2871, 2009. . </li>\n\
							</ul>\
							<p>The proportion of true\
							 nulls is estimated using the Blanchard-Roquain 1-stage procedure with parameter lambda,\
                             via the formula\n\
 
				 estimated pi<sub>0</sub> = ( m - R(alpha,lambda)  + 1) / ( m*( 1 - lambda * alpha ) )\n\
 
				 where R(alpha,lambda) is the number of hypotheses rejected by the BR 1-stage procedure,
				 alpha is FDR level control for this procedure and lambda a 
 				 parameter belonging to (0, 1/alpha) with default value 1. Independence of p-values is assumed.
 				 This estimate may in some cases be larger than 1; it is truncated to 1 if the parameter truncated=TRUE.
                 The estimate is used in the Blanchard-Roquain 2-stage step-up (with the non-truncated version)</p>",
					parameters=list(pValues=list(type="numeric"),alpha=list(type="numeric"),
									lambda=list(type="numeric", default=1), truncate=list(type="logical", default=TRUE))
			)) }

