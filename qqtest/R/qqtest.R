#' @title \code{qqtest} A self-calibrated quantile-quantile plot for assessing distributional shape.
#'
#' @description Draws a quantile-quantile plot for visually assessing whether the data come from a test distribution that has been defined in one of many ways.
#'    The vertical axis plots the data quantiles, the horizontal those of a test distribution.
#'    Interval estimates and exemplars provide different comparative information to assess the evidence provided by the qqplot against the hypothesis that the data come from the test distribution (default is normal or gaussian).  Interval estimates provide test information related to individual quantiles, exemplars provide test information related to the shape of the quantile quantile curve.
#'    Optionally, a visual test of significance (a lineup plot) can be displayed to provide a coarse level of significance for testing the null hypothesis that the data come from the test distribution.
#'    The default behaviour generates 1000 samples from the test distribution and overlays the plot with pointwise interval estimates for the ordered quantiles from the test distribution.  Various option choices are available to effect different visualizations of the uncertainty surrounding the quantile quantile plot.  These include overlaying independently generated exemplar test distribution sample quantile traces so as to assess the joint (as opposed to pointwise) distribution of quantiles. See argument descriptions and examples for more details.
#'
#' @export qqtest
#'
#' @importFrom graphics Axis abline box lines par plot points polygon
#' @importFrom grDevices adjustcolor
#' @importFrom stats dchisq fivenum pchisq ppoints qchisq qexp qlnorm qnorm qt quantile qunif rbeta rchisq rexp rlnorm rnorm rt runif
#'
#' @param data A univariate dataset to be tested. If data has more than one column, the first is used.
#' @param dist The name of the distribution against which the comparison is made, the test distribution for a few built-in distributions.  One of \code{"gaussian"} (or \code{"normal"}), \code{"log normal"},\code{"half normal"}, \code{"uniform"},\code{"exponential"},\code{"student"}, \code{"chi-squared"}, or \code{"kay"}.  Only the first three characters of any of these is needed to specify the dist.  If dist is \code{"student"},  \code{"chi-squared"}, or \code{"kay"}, then a value for the degrees of freedom argument (\code{df} below) is also required.
#' @param df Degrees of freedom of \code{dist} to be used when \code{dist} is either \code{"student"} or \code{"chi-squared"}.
#' @param qfunction If non-\code{NULL}, this must be a function of a single argument (a proportion, p say) which will be used to calculate the quantiles for the test distribution.  If non-\code{NULL}, the \code{rfunction} should also be non-\code{NULL}.  The value of the \code{dist} argument will be ignored when this is the case.
#' @param rfunction If non-\code{NULL}, this must be a function of a single argument (a count, n say) which will be used to randomly select a sample of size n from the test distribution.   If non-\code{NULL}, the \code{qfunction} must also be non-NULL.  If  \code{qfunction} is non-\code{NULL} and  \code{rfunction} is \code{NULL}, then \code{qfunction} will be applied to the output of a call to \code{runif} in place of the \code{NULL} \code{rfunction} (i.e. a probability integral transform is used to generate a random sample).  The value of the \code{dist} argument will be ignored whenever \code{qfunction} is a function.
#' @param dataTest If non-\code{NULL}, this must be a second data set.  The empirical distribution given by this data will be used as the test distribution against which the value of data will be tested.   If non-\code{NULL}, the values of the arguments \code{dist}, \code{qfunction}, and \code{rfunction} will all be ignored in favour of using this empirical distribution as the test distribution.
#' @param p If non-\code{NULL}, this must be a vector containing the probability points at which the quantiles are to be calculated.  If the length of this vector is the same as that of the data, then \code{p} is taken to be the probabilities that correspond to the data; otherwise the data are taken to provide and empirical quantile function values of which are  taken at \code{p} to produce the plot.  If \code{p} is \code{NULL} (the default), then \code{p} is determined by the function \code{ppoints(n,a)} where \code{n} is the number of data points in \code{data} and \code{a} is given by the argument \code{a} below.
#' @param a This is the second parameter given to the \code{ppoints} call to determine \code{p} as necessary.  If \code{a} is \code{NULL}, then default values are chosen for each specific distribution (e.g. 3/8 for gaussian, 0 for uniform, etc.) and as 2/5 otherwise, following the recommendations of C. Cunnane (1978). Users may supply their own value such as the original Hazen's 1/2 or Tukey's 1/3.
#' @param np This is required if the vector \code{p} is provided.  Because \code{p} could take any values, we need to know \code{np} the sample size that was used  to construct the  quantiles at the provided vlaues \code{p}.  When \code{p} and \code{np} are provided all simulation is based on simulating values from the distribution of the order statistics \code{p*np}.
#' @param xAxisAsProbs If \code{TRUE} (the default is FALSE) the horizontal axis will be labelled as probabilities.  These are the cumulative probabilities according to the test distribution.  They are located at the corresponding quantile values.  They are handy in comparing percentiles of the test and data distributions as well as giving some measure of the symmetry and tail weights of the test distribution by their location.  If \code{FALSE} the axis is labelled according to the quantile values.
#' @param yAxisAsProbs If \code{TRUE} (the default is FALSE) the vertical axis will be labelled as probabilities.  These are the cumulative probabilities according to the empirical distribution of the \code{data}.  They are located at the corresponding quantile values.  They are handy in comparing percentiles of the test and data distributions as well as giving some measure of the symmetry and tail weights of the \code{data} distribution by their location. If \code{FALSE} the axis is labelled according to the quantile values.
#' @param xAxisProbs A vector of probabilities to be used to label the x axis ticks when \code{xAxisAsProbs} is \code{TRUE}.  Default is \code{c(0.05, 0.25, 0.50, 0.75, 0.95)}.  Ignored if \code{xAxisAsProbs} is \code{FALSE}.
#' @param yAxisProbs A vector of probabilities to be used to label the y axis ticks when \code{yAxisAsProbs} is \code{TRUE}.  Default is \code{c(0.05, 0.25, 0.50, 0.75, 0.95)}.  Ignored if \code{yAxisAsProbs} is \code{FALSE}.
#' @param nreps The number of replicate samples to be taken from the test distribution to construct the pointwise intervals for each quantile.  Default is 1000.  From these samples, an empirical distribution is generated from the test distribution for the ordered quantiles corresponding to the values of\code{ppoints(length(data))}.  These are used to construct central intervals of whatever proportions are given by \code{centralPercents}.
#' @param centralPercents The vector of proportions determining the central intervals of the empirical distribution of each ordered quantile from the test distribution.  Default is \code{c(0.90, 0.95, 0.99)} corresponding to central 90, 95, and 99\% simulated pointwise confidence intervals for each quantile coming from the test distribution for a sample the same size as \code{data}.The quality of these interval locations typically increases with \code{nreps} and decreases with the probability used for each interval.
#' @param envelope If \code{TRUE} (the default), a grey envelope is plotted showing the central intervals for each quantile as a shade of grey.  The higher is the corresponding probability associated with the interval, the lighter is the shade.  The outermost edges of the envelope are the range of the simulated data from the test distribution.  The envelope thus provides a (pointwise) density estimate of the quantiles drawn from the test distribution for this sample size.   If \code{FALSE} no envelope is drawn.
#' @param drawPercentiles If \code{TRUE}, a pair of curves is plotted to show each of the central intervals as a different line type.  These are plotted over the envelope if \code{envelope} is \code{TRUE}. If \code{FALSE} (the default) no simulated percentile curves are drawn.
#' @param drawQuartiles If \code{TRUE}, a pair of curves is plotted to show the quartiles (central 50\% region) of the ordered quantiles simulated from the test distribution.  The median of these is also plotted as a solid line type. These are plotted over the envelope if \code{envelope} is \code{TRUE}.  If \code{FALSE} (the default) none of these curves are drawn.
#' @param legend If \code{TRUE} (the default is \code{NULL}) with \code{nreps>0} a legend for the appearance of the simulated ranges of the central intervals is added to the plot.  If \code{FALSE}, no legend appears. If \code{NULL}, legend always appears except when \code{lineup = TRUE}.
#' @param nexemplars (default is 0) The number of replicate samples to be taken from the test distribution and plotted as a coloured trail on the qqplot.  Each such trail is a sample of the same size as \code{data} but truly coming from the test distribution.  Each trail gives some idea of what the shape of a qqplot would be for a sample of that size from the test distribution. Together, they give some sense of the variability in the plot's shape.
#' @param typex (default is "o", or match \code{type} if supplied) The \code{plot} type to be used in the plotting of the exemplars, legal values are the same as for the  \code{type} argument of \code{plot}.
#' @param plainTrails If \code{TRUE}, then a single grey colour is used for all exemplar trails.  If \code{FALSE} (the default), each exemplar trail is shown in a different colour.
#' @param colTrails Colours to be used for the trails (default will be multi-colours).
#' @param alphaTrails The alpha transparency to be used in plotting all exemplar trails. The default is 0.25.  Because the trails will over plot, a smaller \code{alphaTrails} value is recommended as \code{nexemplars} increases.
#' @param lwdTrails The graphical line width (\code{lwd}) to be used in plotting all exemplar trails.   The default is 1.  Because the trails will over plot, combining a larger \code{lwdTrails} with \code{envelope = FALSE}, a lower \code{alphaTrails} value larger \code{nexemplars} can give a truer sense of the density of qqplot configurations than with \code{envelope = TRUE}.
#' @param lineup If \code{TRUE} (default is \code{FALSE}) the qqplot of \code{data} is randomly located in a grid of \code{nsuspects} plots.  Identical arguments are given to construct all \code{qqtest} plots in the grid.  Assuming the viewer has not seen the qqplot of this \code{data} before, a successful selection of the true \code{data} plot out of the grid of plots corresponds to evidence against the hypothesis that the \code{data} come from the test distribution.  Significance level is 1/\code{nsuspects}.  Each plot is given a suspect number from 1 to \code{nsuspects} (left to right, top to bottom).  The suspect number of the plot corresponding to the actual \code{data} is returned, slightly obfuscated to help keep the test honest.
#' @param nsuspects The total number of plots (default is 20) to be viewed in the lineup display when \code{lineup} is \code{lineup}.
#' @param col If non-\code{NULL}, \code{col} must be colour to be used for the points in the plot.  If \code{NULL} (the default), an \code{hcl} colour will be used from the values of the arguments \code{h}, \code{c}, \code{l}, and \code{alpha}.
#' @param h The hue of the colour of the points.  Specified as an angle in degrees from 0 to 360 around a colour wheel.  E.g. 0 is red, 120 green, 240 blue,  Default is 260 (a bluish).
#' @param c The chroma of the colour of the points.  Takes values from 0 to an upper bound that is a function of hue, \code{h},  and luminance, \code{l}.  Roughly, for fixed \code{h} and \code{l} the higher the value of \code{c} the greater the intensity of colour.
#' @param l The luminance of the colour of the points.  Takes values from 0 to 100. For any given combination of hue, \code{h},  and chroma, \code{c}, only a subset of this range will be possible.  Roughly, for fixed \code{h} and \code{c} the higher the value of \code{l} the lighter is the colour.
#' @param alpha The alpha transparency of the colour of the points.  Takes values from 0 to 1. Values near 0 are more transparent, values near 1 (the default) are more opaque. Alpha values sum when points over plot, giving some indication of density.
#' @param cex The graphical parameter \code{cex} for the size of the points.
#' @param pch The graphical parameter \code{pch} for the point character to be used for the points.  Default is 19, a filled circle.
#' @param type The graphical parameter \code{type} for the points of the qqplot. Default is "p".
#' @param main The graphical parameter \code{main} providing a title for the plot.  If \code{NULL} (the default), the title will be "qqtest" when \code{lineup = FALSE} and "Suspect: " followed by the plot number when \code{lineup = TRUE}.  An empty string will suppress the title.
#' @param xlab The graphical parameter \code{xlab} labelling the x axis of the plot.  If \code{NULL} (the default), an \code{xlab} is created based on the information available from the other arguments to \code{qqtest} about the test distribution.  An empty string will suppress the labelling.
#' @param ylab The graphical parameter \code{ylab} labelling the y axis of the plot.  If \code{NULL} (the default), a \code{ylab} is created based on the information available from the other arguments to \code{qqtest}. An empty string will suppress the labelling.
#' @param xlim The graphical parameter \code{xlim} determining the display limits of the x axis.
#' @param ylim The graphical parameter \code{ylim} determining the display limits of the y axis.
#' @param axes The graphical parameter \code{axes} determining whether axes are to be displayed (default is \code{NULL} which the same as \code{TRUE} except when \code{lineup=TRUE}, then \code{axes} is \code{FALSE}).
#' @param bty The graphical parameter \code{bty} determining the type of box to be used to enclose the plot (default is \code{"o"}, set \code{bty="n"} for no box).
#' @param ... Any further graphical parameters to be passed to the \code{plot} function.
#'
#' @return Displays the qqplot.  If \code{lineup} is \code{TRUE}, it returns a list with the location (\code{TrueLoc}) of the plot that corresponds to \code{data} encoded as a string whose contents need to be evaluated.  This provides some simple obfuscation of the true location so that the visual assessment can be honest.
#' @source
#' "Self calibrating quantile-quantile plots",
#' R. Wayne Oldford, The American Statistician, 70, (2016)
#'
#' "Unbiased Plotting Positions -- A Review",
#' C. Cunnane, Journal of Hydrology, Vol. 37 (1978), pp. 205-222.
#'
#' @examples
#' #
#' # default qqtest plot
#' qqtest(precip, main = "Precipitation (inches/year) in 70 US cities")
#' #
#' # qqtest to compare to qqnorm
#' op <- par(mfrow=c(1,2))
#' qqnorm(precip, main="qqnorm")
#' qqtest(precip, main="qqtest",
#'        xAxisAsProbs=FALSE, yAxisAsProbs=FALSE)
#' par(op)
#' #
#' #  Use lines instead of envelope
#' qqtest(precip, envelope=FALSE, drawPercentiles=TRUE,
#'        main = "Precipitation (inches/year) in 70 US cities")
#' #
#' #  Use quartiles instead of envelope
#' qqtest(precip, envelope=FALSE, drawQuartiles=TRUE,
#'        main = "Precipitation (inches/year) in 70 US cities")
#' #
#' #  Use coloured exemplars (qqplot of data simulated from the test distribution)
#' #  and suppress the envelope.  Where the envelope, percentiles, and quartiles are
#' #  simulated pointwise bands, exemplars give some sense of what the (joint) shape of the
#' #  quantile-quantile plot should look like (for data from the test distribution).
#' #  Each simulated sample is a different colour.
#' qqtest(precip, nexemplars=10, typex="o", envelope=FALSE, type="p",
#'        main = "Precipitation (inches/year) in 70 US cities")
#' #
#' #  Alternatively, the trail of each exemplar could be plain (the identical grey).
#' #  Making each trail wide and assigning it some transparency (alpha near 0)
#' #  allows the trails to give a sense of the density through the darkness of the grey.
#' #
#' qqtest(precip, nexemplars=20, envelope=FALSE,
#'        lwdTrails=3, plainTrails=TRUE, alphaTrail=0.4, typex="o", type="o",
#'        main = "Precipitation (inches/year) in 70 US cities")
#' #
#' #  Wide coloured exemplars with some transparency provide an indication of
#' #  density and allow some trails to be followed by colour.
#' #
#' qqtest(precip, nexemplars=20, envelope=FALSE,
#'        lwdTrails=3,  alphaTrail=0.4, typex="o", type="o", col="black",
#'        main = "Precipitation (inches/year) in 70 US cities")
#'
#'
#' #  Envelope and exemplars with coloured trails to be followed.
#' #
#' qqtest(precip, nexemplars=5,
#'        lwdTrails=2, alphaTrail=0.6, alpha=0.8,
#'        main = "Precipitation (inches/year) in 70 US cities")
#' #
#' #
#' #  gaussian - qqplot, but now showing in the line up
#' result <- qqtest(precip, lineup=TRUE, main="Suspect", legend=FALSE,
#'                  cex=0.75, col="grey20", ylab="", pch=21)
#' # the location of the real data in the line up can be found by evaluating
#' # the contents of the string
#'  result$TrueLoc
#' #
#' # Cut and paste the string contents into the R console, or evaluate
#'  eval(parse(text=result$TrueLoc))
#' #
#' #
#' # log-normal ... using the bacteria data from Whipple(1916)
#' data(bacteria, package="qqtest")
#' # Note that these are selected percentiles from a sample of 365 days in a year
#' with(bacteria,
#'     qqtest(count, dist = "log-normal", p=percentTime/100, np=365, type="o",
#'  		  yAxisAsProbs=FALSE, ylab="bacteria per cc",
#'            xAxisProbs = c(0.01, 0.50,0.75, 0.90, 0.95, 0.99, 0.995),
#'            xlab="Percentage of days in 1913",
#'            main = "Number of bacteria from the Delaware river in 1913")
#'     )
#' ptics <- c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99 )
#' axis(1,at=qnorm(ptics), labels=floor(ptics*100))
#' yvals <- c(100, 1000, 10000, 100000)
#' axis(2, at=log(yvals,10),
#'      labels=c("100", "1,000", "10,000", "100,000"))
#' #
#' # compare this to the log-scaled normal qqplot
#' #
#' #
#' with(bacteria,
#'     qqtest(log(count, 10), dist = "normal",
#'            p=percentTime/100, np=365,
#'  		  type="o", axes=FALSE,
#'            ylab="bacteria per cc",
#'            xlab="Proportion of days in 1913",
#'            main = "Number of bacteria from the Delaware river in 1913")
#'     )
#' #

#' #
#' # Half normal ... using the penicillin data from Daniel(1959)
#' data(penicillin)
#'
#' qqtest(penicillin, dist = "half-normal")
#'
#' # Or the same again but with significant contrast labelled
#'
#'
#' with (penicillin,
#' 	{qqtest(value, yAxisProbs=c(0.1, 0.75, 0.90, 0.95),
#'          dist="half-normal",
#' 			ylab="Sample cumulative probability",
#'          xlab="Half-normal cumulative probability")
#' 	 ppAdj <- (1+ppoints(31))/2  # to get half-normals from normal
#' 	 x <- qnorm(ppAdj)
#' 	 valOrder <- order(value)    # need data and rownames in increasing order
#' 	 y <- value[valOrder]
#' 	 tags <- rownames(penicillin)[valOrder]
#' 	 selPoints <- 28:31          # going to label only the largest effects
#'   xoffset <- c(0.01, 0.02, 0.03, 0.075)  # text function is a bit off
#' 	 text(x[selPoints]-xoffset, y[selPoints],
#'        tags[selPoints],
#'        pos=2, cex=0.75)
#' 	}
#' )
#' #
#' # K on 9 df  ... see help(dkay)
#' # Use data on primer paint thickness (standard deviations on n=10)
#' data(primer, package="qqtest")
#' with (primer,
#' 	     qqtest(s,  dist="kay", df=9,
#' 		        yAxisAsProbs=FALSE,
#' 			    ylab="Standard deviation of primer thickness (in mils)")
#' 		)
#' #
#' # chi-squared on 3 df
#' # Use robust covariance matrix in calculation Mahalanobis distances
#' # for the classical Brownlee stackloss data.
#' data(stacklossDistances, package="qqtest")
#' with(stacklossDistances,
#'      qqtest(robust, dist="chi", df=3, ylab="Robust Mahalanobis distances"))
#' #
#' #
#' # user supplied qfunction and rfunction -- compare to beta distribution
#' qqtest(precip,
#'        qfunction=function(p){qbeta(p, 2, 2)},
#'        rfunction=function(n){rbeta(n, 2, 2)},
#'        main = "Precipitation (inches/year) in 70 US cities")
#' #
#' #
#' # user supplied qfunction only -- compare to beta distribution
#' qqtest(precip,
#'        qfunction=function(p){qbeta(p, 2, 2)},
#'        main = "Precipitation (inches/year) in 70 US cities")
#' #
#' # comparing data samples
#' #
#' # Does the sample of beaver2's temperatures look like they
#' # could have come from a distribution shaped like beaver1's?
#' #
#'  	qqtest(beaver2[,"temp"],
#' 		       dataTest=beaver1[,"temp"],
#' 		       ylab="Beaver 2", xlab="Beaver 1",
#' 		       main="Beaver body temperatures")
#' #
#' #
#' # For the famous iris data, does the sample of iris versicolor
#' # appear to have the same (marginal) distributional shape
#' # as does that of iris virginica (to which it is more closely related)?
#' #
#' op <- par(mfrow=c(2,2))
#' with(iris, {
#' 	qqtest(Sepal.Length[Species=="versicolor"],
#' 		   dataTest= Sepal.Length[Species=="virginica"],
#' 		   ylab="versicolor", xlab="virginica",
#' 		   main="Sepal length")
#' 	qqtest(Sepal.Width[Species=="versicolor"],
#' 		   dataTest= Sepal.Width[Species=="virginica"],
#' 		   ylab="versicolor", xlab="virginica",
#' 		   main="Sepal width", legend=FALSE)
#' 	qqtest(Petal.Length[Species=="versicolor"],
#' 		   dataTest=Petal.Length[Species=="virginica"],
#' 		   ylab="versicolor", xlab="virginica",
#' 		   main="Petal length", legend=FALSE)
#' 	qqtest(Petal.Width[Species=="versicolor"],
#' 		   dataTest= Petal.Width[Species=="virginica"],
#' 		   ylab="versicolor", xlab="virginica",
#' 		   main="Petal width", legend=FALSE)
#' 	}
#' 	)
#' par(op)



qqtest <- function (data,
			        dist=c("gaussian", "normal",
			        	   "log-normal", "half-normal",
			        	   "uniform", "exponential",
			        	   "chi-squared", "kay",
			        	   "student","t"),
			        df=1,
			        qfunction=NULL,
			        rfunction=NULL,
			        dataTest=NULL,
				      p=NULL,
				      a=NULL,
				      np = NULL,
			        xAxisAsProbs = FALSE,
			        yAxisAsProbs = FALSE,
			        xAxisProbs =  c(0.05, 0.25, 0.50, 0.75, 0.95),
			        yAxisProbs = c(0.05, 0.25, 0.50, 0.75, 0.95),
			        nreps=1000,
			        centralPercents=c(0.90, 0.95, 0.99),
			        envelope=TRUE,
			        drawPercentiles = FALSE,
			        drawQuartiles = FALSE,
			        legend=NULL,
			        nexemplars=0,
			        typex=NULL,
			        plainTrails=FALSE,
				      colTrails=NULL,
				      alphaTrails=0.25,
			        lwdTrails=1,
			        lineup=FALSE,
			        nsuspects=20,
			        col=NULL,
			        h=260,
			        c=90,
			        l=60,
			        alpha=1.0,
			        cex=NULL,
			        pch=19,
			        type=NULL,
				      main=NULL,
				      xlab=NULL,
			        ylab=NULL,
			        xlim=NULL,
			        ylim=NULL,
			        axes=NULL,
			        bty="o",
			        ...
                      ) {

    if (is.matrix(data) || is.data.frame(data)) {
            data <- data.matrix(data)[,1]
            }

    if (is.null(typex)){
    	typex <- if (is.null(type))  "o" else type
    }
    if (is.null(type)) type <- "p"

    if (!is.numeric(alphaTrails)) {stop("alphaTrails must be a number from 0 to 1.")} else {
      if ((alphaTrails < 0) | (alphaTrails > 1)) stop("alphaTrails must be a number from 0 to 1.")
    }

    if (nexemplars > 0) {
      if (is.null(colTrails)) {
        colTrails <- sapply(1:nexemplars, function(i) {grDevices::hcl(h=i*360/nexemplars, c=90, l=60)})
      }
      colTrails <- adjustcolor(colTrails, alpha.f = alphaTrails)
      alphaTrails <- 1
      colTrails <- rep(colTrails, length.out=nexemplars)
    }


	generateUniformOrderStat <- NULL
     # If probabilities p are provided, only these values will be used.
	if (!is.null(p) ){
		# first check whether these are the same as the number of data points
     	# if so, then assume that the data points *are* the values at these places.
     	# if not, then extract from the data the empirical quantiles at the given p values
		if (length(p)!=length(data)) {data <- quantile(data,p)}
		#
		if (envelope||drawPercentiles||drawQuartiles||lineup||(nexemplars>0)) {
			# This means we are going to have generate some values but cannot count
			# on the length(data) or length(p) to determine the sample size on which
			# to base the generation.
			# In this case, the sample size must be provided as the value of np:

			if (is.null(np) || !is.numeric(np)) {
				stop("The sample size np, on which values of p are to be based, must be given whenever the values of p are specified.")}


			# Because we have specified probabilities, we are going to have to generate
			# the corresponding order statistics directly.
			# Plan is to generate these by applying the appropriate quantile function to
			# the corresponding order statistic from a uniform ... Q(U_(i))
			# Need to determine an appropriate index.
			# We will take the sample size N to be large enough to have each prob in p
			# be i/N for large enough N (Note fractional indices are possible)
			orderIndices <- np * p
			# We generate the uniform order statistics
			generateUniformOrderStat <- function(index, thismany) {
				rbeta(thismany, index, np + 1 - index)}
			}
		}


     n <- length(data)

	# Get the distribution argument and fix synonyms
	dist <- match.arg(dist)
	# normal and gaussian are the same
	if (dist=="normal") {dist <- "gaussian"}
	# t and student are the same
	if (dist=="t") {dist <- "student"}

	# The following stay null if there are none to create
	reps <- NULL
	exemplars <- NULL

	# Has the user supplied some axis labels?
	makeXlabel <- is.null(xlab)
	makeYlabel <- is.null(ylab)
	#
	# Has the user supplied a main title?
	makeTitle <- is.null(main)
	#

	# If lineup, don't do all reps yet
	if (lineup) {
		if (!(nsuspects > 1)) {
			warning("nsuspects must be > 1, set to 20 instead")
			nsuspects <- 20
			}
		nrepsInput <- nreps
		nreps <- nsuspects
		nexemplarsInput <- nexemplars
		nexemplars <- 0   ##### Should be zero? or 1?
		if (is.null(axes)) axes <- FALSE
		if (makeYlabel) ylab <- ""
		if (is.null(legend)) legend <- FALSE
		if (is.null(cex)) cex <- 2 / floor(sqrt(nsuspects))
	} else {
		# No lineup? Default axes = TRUE
		if(is.null(axes)) axes <- TRUE
	}


 	#
	# Getting some reps from the appropriate distribution

	if (!is.null(dataTest)) {
		## plan is to test the data against the
		## empirical distribution from dataTest
		if (is.null(p)){
			if (is.null(a)){
			## use general purpose a=2/5 as per Cunnane (1978)
				p <- ppoints(n,2/5)
			}
			else {p <- ppoints(n,a)}
		}

		qfunction <- function(p){quantile(dataTest,p)}
		rfunction <- function(n){sample(dataTest,n,replace=TRUE)}
	    if (makeXlabel) {xlab <- "Empirical distribution"}
	}
	else {
	if (is.null(qfunction)){
		switch(dist,
				# Gaussian or Normal
				"gaussian"={if (is.null(p)){
								if (is.null(a)){
								##  a=3/8 as per Cunnane (1978) for Normal
									p <- ppoints(n,3/8)
									}
							else { p <- ppoints(n,a) }
							}
							qfunction <- qnorm
							rfunction <- rnorm
							if (makeXlabel) {xlab <- "Gaussian"}
							},
				"log-normal"={if (is.null(p)){
								if (is.null(a)){
									##  a=3/8 as per Cunnane (1978) for Normal
									p <- ppoints(n,3/8)
									}
								else { p <- ppoints(n,a) }
								}
							qfunction <- qlnorm
							rfunction <- rlnorm
							if (makeXlabel) {xlab <- "Lognormal"}
							},
				"half-normal"={if (is.null(p)){
								if (is.null(a)){
									##  a=3/8 as per Cunnane (1978) for Normal
									p <- ppoints(n,3/8)
									}
								else { p <- ppoints(n,a) }
								}
							qfunction <- function(p) {qnorm((1 + p)/2)}
							rfunction <- function(n) {abs(rnorm(n))}
							if (makeXlabel) {xlab <- "Half-normal"}
							},
				"uniform"={if (is.null(p)){
								if (is.null(a)){
									##  a=0 as per Cunnane (1978) for uniform
									p <- ppoints(n,0)
									}
								else { p <- ppoints(n,a) }
								}
							qfunction <- qunif
							rfunction <- runif
							if (makeXlabel) {xlab <- "Uniform"}
							},
				"exponential"={if (is.null(p)){
								if (is.null(a)){
									##  a=0.44 as per Cunnane (1978) for exponential
									p <- ppoints(n,0.44)
									}
								else { p <- ppoints(n,a) }
								}
							qfunction <- qexp
							rfunction <- rexp
							if (makeXlabel) {xlab <- "Exponential(1)"}
							},
				# K distribution (= sqrt(chisquared/df))
				"kay"=		{if (is.null(p)){
								if (is.null(a)){
									if (df <= 1){
										##  a=1/2 as per Cunnane (1978) for (Pearson III or chi)
										p <- ppoints(n,1/2)
										}
									else {
										##  a=2/5 as per Cunnane (1978) for (Pearson III or chi)
										p <- ppoints(n,2/5)
										}
									}
								else { p <- ppoints(n,a) }
								}
							qfunction <- function(p) {qkay(p,df)}
							rfunction <- function(n) {rkay(n,df)}
							if (makeXlabel) {xlab <- paste("K(",df,")", sep="")}
							},
				"chi-squared"={if (is.null(p)){
								if (is.null(a)){
									if (df <= 1){
										##  a=1/2 as per Cunnane (1978) for (Pearson III or chi)
										p <- ppoints(n,1/2)
										}
									else {
										##  a=2/5 as per Cunnane (1978) for (Pearson III or chi)
										p <- ppoints(n,2/5)
										}
									}
								else { p <- ppoints(n,a) }
								}
							qfunction <- function(p) {qchisq(p,df)}
							rfunction <- function(n) {rchisq(n,df)}
							if (makeXlabel) {xlab <- paste("Chi-squared(",df,")", sep="")}
							},
				# Student's t
				"student"=	{if (is.null(p)){
								if (is.null(a)){
									## use general purpose a=2/5 as per Cunnane (1978)
									p <- ppoints(n,2/5)
									}
								else { p <- ppoints(n,a) }
								}
							qfunction <- function(p) {qt(p,df)}
							rfunction <- function(n) {rt(n,df)}
							if (makeXlabel) {xlab <- paste("Student t(",df,")", sep="")}
							},
				stop(paste("Unimplemented distribution:", dist,
						  "Could instead call qqtest by",
				          "supplying a function for qfunction",
				          "(and optionally another for rfunction)."))
				)  # end of switch

	} else if(is.function(qfunction)) {
		if (is.null(p)){
			if (is.null(a)){
			## use general purpose a=2/5 as per Cunnane (1978)
				p <- ppoints(n,2/5)
			} else
			{
				p <- ppoints(n,a)
			}
		}
	    if (!is.function(rfunction)) {
	    	rfunction <- function(n) {
	    		sapply(runif(n),qfunction)
	    	}
	    }
	    if (makeXlabel) {xlab <- "Hypothetical distribution"}

	    } else {warning("qfunction must be a function that returns quantiles")}

	}


 q <- qfunction(p)
 xAxisTicksAt <- qfunction(xAxisProbs)

 # Get the location and scale changes to use for the simulations (and to return).
 # Use a robust line fit to get the location scale correction for the
 # data sampled from the test distribution	(tukeyline not good enough)

 line <- robust::lmRob(sort(data) ~ 1 + q)
 loc <- line$coefficients[1]
 scale <- line$coefficients[2]

 if (nreps > 0) {
 	if (is.null(generateUniformOrderStat)) {
 		reps <-  rfunction(n*nreps)
 	}
 	else
 	{
 		reps <- qfunction(sapply(orderIndices,
 								function(index){generateUniformOrderStat(index,nreps)}
 								)
 						 )
 	}
 	reps <- loc + scale * array(reps,dim=c(nreps,n))
 }

 if (nexemplars > 0) {
 	if (is.null(generateUniformOrderStat)) {
 		exemplars <-  rfunction(n*nexemplars)
 	}
 	else
 	{
 		exemplars <- qfunction(sapply(orderIndices,
 								function(index){
 									generateUniformOrderStat(index, nexemplars)}))
 	}
 	exemplars <- loc + scale *  array(exemplars,dim=c(nexemplars,n))
 	}


 # Do lineup plot if asked. ... recursive
 if (lineup){
 	suspects <- reps[1:nsuspects,]
 	trueLoc <- sample(1:nsuspects,1)
 	suspects[trueLoc,] <- data
 	nrow <- floor(sqrt(nsuspects))
 	ncol <- ceiling(sqrt(nsuspects))
 	parOptions <- par(mfrow=c(nrow,ncol), mar=c(0,0,1.5,0))

 	if (makeTitle){main <- "Suspect"}

 	for (i in 1:nsuspects) {
 		qqtest(data = suspects[i,],
                      dist=dist,
                      df=df,
                      qfunction = qfunction,
                      rfunction = rfunction,
                      dataTest = dataTest,
                      xAxisProbs = xAxisProbs,
                      yAxisProbs = yAxisProbs,
                      xAxisAsProbs = xAxisAsProbs,
                      yAxisAsProbs = yAxisAsProbs,
                      nreps= nrepsInput,
                      centralPercents = centralPercents,
                      envelope = envelope,
                      drawPercentiles = drawPercentiles,
                      drawQuartiles = drawQuartiles,
                      legend = legend,
                      nexemplars = nexemplarsInput,
                      typex = typex,
                      plainTrails = plainTrails,
 		                  colTrails = colTrails,
 		                  alphaTrails = alphaTrails,    # Note this will be 1 here so that repeated application doesn't make the colours disappear in lineup
                      lwdTrails = lwdTrails,
                      lineup=FALSE,
                      nsuspects=1,
                      col=col,
                      h=h,
                      c=c,
                      l=l,
                      alpha=alpha,
                      cex=cex,
                      pch=pch,
                      type=type,
                      main =  if(main=="") {paste0(i)} else {paste(main, i, sep=": ")},
 		                  xlab=xlab,
                      ylab = ylab,
                      xlim = xlim,
                      ylim = ylim,
                      axes=axes,
                      bty=bty,
                      ...)
 	}
 	par(parOptions)
 	base <- sample(2:2*nsuspects,1)
 	offset <- sample(5:5*nsuspects,1)
 	list(TrueLoc = paste("log(",base^(trueLoc + offset),", base=",base,") - ",offset,sep=""))
 } else {
 	# Otherwise construct the plot.

 if (!is.null(reps)) {
 	# Sort the reps into the right form and get summary statistics
 	reps <- apply(reps, 1, sort)      # Note this flips the dimensions of reps
 	Nums <- apply(reps, 1, fivenum)

	 if (envelope||drawPercentiles){
 		nLevels <- length(centralPercents)
 		centralPercents <- sort(centralPercents)
 		SymmetricAdjust <-(1-centralPercents)/2
		bottomPcts <- (1-centralPercents)-SymmetricAdjust
 		topPcts <-  centralPercents + SymmetricAdjust
 		Pcts <- c(rev(bottomPcts),topPcts)
 		NumsPct <- apply(reps, 1, function(x){quantile(x,Pcts)})

		# defining grey
 		greyhue <- 260
 		greyluminance <- 65
 		greychroma <- 0
 		greyalpha <- 0.6/(nLevels +1)  # plus 1 for the range
 		grey <- grDevices::hcl(h= greyhue, c= greychroma, l= greyluminance,
 								alpha= greyalpha)

	 	}
 }

 ## Get the base plot
 xlim <- if(!is.null(xlim)) {xlim} else {range(q)}
 ylim <- if(!is.null(ylim)) {
 				ylim
 				} else {
 				    yrange <- range(data)
 					if (!is.null(reps)) {
 						yrange <- range(c(yrange, range(reps)))
 						}
 					if (!is.null(exemplars)) {
 						yrange <- range(c(yrange, range(exemplars)))
 						}
 					yrange
 				  }


 plot(0,0,
       main = if(makeTitle) {"qqtest"} else main,
       xlab = if(makeXlabel) {
       			paste(xlab, if(xAxisAsProbs) {
       						"cumulative probability (on quantile scale)"
       						} else {"quantiles"})
       			} else xlab,
       ylab = if(makeYlabel) {
       			paste(if(yAxisAsProbs) {
       						"Sample cumulative probability (on quantile scale)"
       						} else {"Sample quantiles"})
       			} else ylab,
       col = "white",
       xlim = xlim,
       ylim = ylim,
       axes = FALSE,
       type="n",
       bty=bty,
       ...             # other plot parameters
       )

#  Draw the axes
if(axes) {
	if (xAxisAsProbs) {
	 Axis(side=1,
	 	  labels = paste(xAxisProbs),
	      at = xAxisTicksAt)
	      }   else  {
     Axis(side=1)
		}
	if (yAxisAsProbs) {
	 Axis(side=2,
	 	  labels = paste(yAxisProbs),
	      at = quantile(data,yAxisProbs))
	      } else {
     Axis(side=2)
		}
}




   # put the box around to look standard
if (!bty=="n") {box(bty=bty)}


 ## Add the envelope for the replicates

if(envelope && !is.null(reps)) {
 # first the range
 polygon(c(q,rev(q)),
         c(Nums[1,],rev(Nums[5,])),
         border=grey,
         col=grey)
 # now the central pointwise "confidence" intervals
 nLevels <- length(centralPercents)
 for (i in 1:nLevels){
 	iLower <- i
 	iUpper <- nLevels*2 - iLower +1
 	polygon(c(q,rev(q)),
         	c(NumsPct[iLower,],rev(NumsPct[iUpper,])),
         	border=grey,
         	col=grey)
	}
  }
 ## Add exemplar trails
 if (!is.null(exemplars)&&nexemplars>=1) {

 	if(!plainTrails){
 			for (i in 1:nexemplars){
    	 		points(q,sort(exemplars[i,]),
              		  col=colTrails[i],
              		  type=typex,
               		  lwd=lwdTrails)
     			}
 		} else { # grey trails
 			plainCol <- grDevices::hcl(h=260, c=0, l= 90, alpha=alphaTrails)
    			for (i in 1: nexemplars){
    	 			points(q,sort(exemplars[i,]),
               		 col=plainCol,
              		 type=typex,
               		 lwd=lwdTrails)
                }
   		}
   }

if(drawPercentiles && !is.null(reps)){
 # and the various percentiles
 if (drawQuartiles) {
 	 lineTypes <- c(3:(nLevels +2))
 	} else{
 		lineTypes <- c(1:nLevels)
 		}

  lineCols <- rep("darkgrey",nLevels)
  for(i in 1:nLevels) {
 	# central intervals
 	lineCols[i] <- "darkgrey"
 	# draw the lower and upper lines
 	lines(q,NumsPct[i,],col= lineCols[i],lty= lineTypes[i], lwd=2)
    lines(q,NumsPct[2*nLevels-i+1,],col= lineCols[i],lty= lineTypes[i], lwd=2)
    }
 }
 # quartiles
 if(drawQuartiles && !is.null(reps)){
   	# draw quartiles
 	lines(q,Nums[2,],col="black",lty=2)
 	lines(q,Nums[4,],col="black",lty=2)
 	# draw median
 	lines(q,Nums[3,],col="black",lty=1)
 	}


 # and optional legend
 if (is.null(legend)) legend <- TRUE
 if (legend && (drawQuartiles||drawPercentiles) && !is.null(reps)){
 	if (drawPercentiles) {
 		if (drawQuartiles) {
 			legendString <- c(paste(signif(100*rev(centralPercents),3), "% central range", sep=""),
 							  "quartiles", "median")
 			legendLineTypes <- c(lineTypes,2,1)
 			legendCols <- c(lineCols,"black","black")
 		} else {
 			# no quartiles
 			legendString <- paste(signif(100*rev(centralPercents),3), "% central range", sep="")
 			legendLineTypes <- lineTypes
 			legendCols <- lineCols
 		   }
 	} else {
 			# only quartiles
 			legendString=c("quartiles", "median")
 			legendLineTypes <- c(2,1)
 			legendCols <- c("black","black")
 		}
    # Draw it
    if (envelope && !drawQuartiles) {
    fillCols <- rep(grey, nLevels+1)
    for (i in 2:(nLevels+1)) {
    	fillCols[i] <- grDevices::hcl(h= greyhue, c= greychroma, l= greyluminance, alpha=i * greyalpha)
    }

    	legend("topleft",
      	  	legend = legendString,
             lwd=1,
             cex=0.8,
             bty="n",
             lty= legendLineTypes,
             col= legendCols,
             fill=fillCols,
             border=fillCols,
             text.col="darkgrey",
             title=paste("Simulated ranges", "n =", nreps))

    } else {
    	legend("topleft",
      	  	legend = legendString,
             lwd=1,
             cex=0.8,
             bty="n",
             lty= legendLineTypes,
             col= legendCols,
             text.col="darkgrey",
             title=paste("Simulated ranges", "n =", nreps))
             }
   }

 if (legend && !drawQuartiles && !drawPercentiles && envelope && !is.null(reps)){
    # Draw it
    fillCols <- rep(grey, nLevels+1)
    for (i in 2:(nLevels+1)) {
    	fillCols[i] <- grDevices::hcl(h= greyhue, c= greychroma, l= greyluminance, alpha=i * greyalpha)
    }
 	legend("topleft",
      	  	legend = c("Range",paste(signif(100*rev(centralPercents),3),
      	  							 "% central range",
      	  							 sep="")),
             col=fillCols,
             fill=fillCols,
             border=fillCols,
             cex=0.8,
             bty="n",
             text.col="darkgrey",
             title=paste("Simulated ranges", "n =", nreps))
             }
   ## And finally the points
if (is.null(col)){
	plot_colour <- grDevices::hcl(h=h, c=c, l=l, alpha=alpha)
	} else {plot_colour <- col}

 points(q,sort(data),
        col=plot_colour,
        type=type,
        pch=pch,
        cex=if(is.null(cex)) 1 else cex
        )

    } #end of else from #lineup


 }