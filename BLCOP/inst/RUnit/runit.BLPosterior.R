test.BLPosterior <- function()
{    
    
    #entries <- c(0.001005,0.001328,-0.000579,-0.000675,0.000121,0.000128,-0.000445,-0.000437 ,
    # 0.001328,0.007277,-0.001307,-0.000610,-0.002237,-0.000989,0.001442,-0.001535 ,
    # -0.000579,-0.001307,0.059852,0.027588,0.063497,0.023036,0.032967,0.048039 ,
    #-0.000675,-0.000610,0.027588,0.029609,0.026572,0.021465,0.020697,0.029854 ,
    # 0.000121,-0.002237,0.063497,0.026572,0.102488,0.042744,0.039943,0.065994 ,
    # 0.000128,-0.000989,0.023036,0.021465,0.042744,0.032056,0.019881,0.032235 ,
    #-0.000445,0.001442,0.032967,0.020697,0.039943,0.019881,0.028355,0.035064 ,
    #-0.000437,-0.001535,0.048039,0.029854,0.065994,0.032235,0.035064,0.079958 )
    #
    #myVarcov2 <- matrix(entries, ncol = 8, nrow = 8)
    #mu <- c(0.08, 0.67,6.41, 4.08, 7.43, 3.70, 4.80, 6.60) / 100
    #pick <- matrix(0, ncol = 8, nrow = 3, dimnames = list(NULL, letters[1:8]))
    #pick[1,7] <- 1
    #pick[2,1] <- -1; pick[2,2] <- 1
    #pick[3, 3:6] <- c(0.9, -0.9, .1, -.1)
    #confidences <- 1 / c(0.00709, 0.000141, 0.000866)
    #myViews <- BLViews(pick, c(0.0525, 0.0025, 0.02), confidences, letters[1:8])
    #myPosterior <- posteriorEst(myViews, tau = 0.025, mu = mu, myVarcov2 )
    
	BLEx <- get(load( file.path(BLCOPOptions("unitTestPath"), "BLExample.RData") ))
	myPosterior <- BLEx$posterior
	
    checkEquals(myPosterior@posteriorMean, c(0.0006845523,0.0049366147,0.0624770045,0.0412809077,0.0728581285,
        0.0375468518, 0.0469830868, 0.0655370385 ), checkNames = FALSE )
}    

test.posteriorFeasibility <- function()
{
	# at the moment this test requires corpcor, due to the expected availability of the shrinkage estimator.  
	# TODO: remove this dependency
	if(!require("corpcor", quiet = TRUE))
	{
		warning("Could not load the corpcor package, these tests will not be run")
		return()
	}
	pick <- matrix(0, 2, 6)
	pick[1,1:3] <- 1/3
	pick[2,4:6] <- 1/3
	views <- BLViews(pick, q = c(0.2, -0.1), confidences = c(100,100), colnames(monthlyReturns))
	post <- BLPosterior(monthlyReturns, views, 1, sp500Returns, US13wTB, 0.9, covEstimator = "cov.shrink")
	feasibility <- posteriorFeasibility(post)
	
	checkEquals(2.261890, feasibility$mahalDist, tolerance = 1e-03)
	checkEquals(0.8941064, feasibility$mahalDistProb, tolerance = 1e-06)
	
}   