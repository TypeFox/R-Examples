#### Distances ####


#' Distance measures (between constructs or elements).
#'
#' Various distance measures between elements or constructs are calculated.
#'
#' @param x           \code{repgrid} object.
#' @param along       Whether to calculate distance for 1 = constructs (default) 
#'                    or for 2= elements.
#' @param dmethod     The distance measure to be used. This must be one of 
#'                    "euclidean", "maximum", "manhattan", "canberra", "binary" 
#'                    or "minkowski". Any unambiguous substring can be given. 
#'                    For additional information on the different types type
#'                    \code{?dist}. 
#' @param  p          The power of the Minkowski distance, in case \code{"minkowski"}
#'                    is used as argument for \code{dmethod}.
#' @param trim        The number of characters a construct or element is trimmed to (default is
#'                    \code{20}). If \code{NA} no trimming occurs. Trimming
#'                    simply saves space when displaying correlation of constructs
#'                    with long names.
#' @param index       Whether to print the number of the construct or element 
#'                    in front of the name (default is \code{TRUE}). This is useful to avoid
#'                    identical row names, which may cause an error.
#' @param ...         Additional parameters to be passed to function \code{dist}.
#'                    Type \code{dist} for further information. 
#' @return            \code{matrix} object.
#'
#' @author            Mark Heckmann
#' @export
#' @examples \dontrun{
#'
#'    # between constructs
#'    distance(bell2010, along=1)
#'    # between elements
#'    distance(bell2010, along=2)
#'  
#'    # several distance methods
#'    distance(bell2010, dm="man")         # manhattan distance
#'    distance(bell2010, dm="mink", p=3)   # minkowski metric to the power of 3
#'
#'    # to save the results without printing to the console
#'    d <- distance(bell2010, trim=7)
#'    d
#'    
#'    # some more options when printing the distance matrix
#'    print(d, digits=5)
#'    print(d, col.index=FALSE)
#'    print(d, upper=FALSE)
#'    
#'    # accessing entries from the matrix
#'    d[1,3]
#'
#' }
#'
distance <- function(x, along=1, dmethod="euclidean", 
                     p=2, trim=20, index=TRUE, ...)
{  
  dmethods <- c("euclidean", "maximum", "manhattan",    # possible distance methods
                "canberra", "binary", "minkowski")
  dmethod <- match.arg(dmethod, dmethods)               # match method
  if (!inherits(x, "repgrid"))   						            # check if x is repgrid object
    stop("Object x must be of class 'repgrid'")
  
  r <- getRatingLayer(x, trim = trim)
  if (along == 2)
    r <- t(r)
  d <- dist(r, method = dmethod, p = p, ...)
  d <- as.matrix(d) 
  d <- addNamesToMatrix2(x, d, index=index, trim=trim, along=along)  
  class(d) <-  c("distance", "matrix")
  attr(d, "arguments") <- list(along=along, dmethod=dmethod, p=p, 
                               notes=NULL, cutoff=NULL)
  return(d)
}


#' Print method for class distance.
#' 
#' @param x           Object of class distance.
#' @param digits      Numeric. Number of digits to round to (default is 
#'                    \code{2}).
#' @param col.index   Logical. Whether to add an extra index column so the 
#'                    column names are indexes instead of construct names. This option 
#'                    renders a neater output as long construct names will stretch 
#'                    the output (default is \code{TRUE}).
#' @param upper       Whether to display upper triangle of correlation matrix only 
#'                    (default is \code{TRUE}).
#' @param cutoffs     Cutoff values. Values below or above this interval are not 
#'                    printed. For Slater distances \code{c(.8, 1.2)} are common 
#'                    values.
#' @param diag        Whether to show the matrix diagonal.                    
#' @param ...         Not evaluated.
#' @export
#' @method            print distance
#' @keywords          internal
#'
print.distance <- function(x, digits=2, col.index=TRUE,
                           upper=TRUE, diag=FALSE, cutoffs=NA, ...)
{
  diag <- !diag                   # convert as used in upper.tri
  args <- attr(x, "arguments")
  d <- x
  class(d) <- "matrix" 
  d <- round(d, digits)     
  e <- format(d, nsmall=digits)   # convert to characters for printing
  
  ## console output ##  
  blank <- paste(rep(" ",  max(nchar(as.vector(e)))), collapse="", sep="")
  
  # remove values above or below explicit cutoff
  if (!is.na(cutoffs[1])) { 
    n.s. <- min(cutoffs) < d & d < max(cutoffs)
    e[n.s.] <- blank    
  }
  
  if (upper)
    e[lower.tri(e, diag=diag)] <- blank
  # make index column for neater colnames
  if (col.index)                                   
    e <- addIndexColumnToMatrix(e) else
      colnames(e) <- seq_len(ncol(e))
  e <- as.data.frame(e)
  if (args$along == 1) { 
    cat("\n############################")
    cat("\nDistances between constructs") 
    cat("\n############################")
  } else if (args$along == 2) {
    cat("\n##########################")
    cat("\nDistances between elements") 
    cat("\n##########################")
  }
  cat("\n\nDistance method: ", args$dmethod, "\n")
  if (args$dmethod == "minkowski")
    cat("power p:", args$p, "\n")
  cat("\n")
  print(e) 
  if (!is.null(args$notes))
    cat(args$notes)
  cat("\n")
}


#### Slater distance ####

#' Internal workhorse for Slater standardization.
#' 
#' Function uses a matrix as input. All overhead
#' of \code{repgrid} class is avoided. Needed for speedy simulations.
#'
#' @param x   A matrix.
#' @keywords internal
#' @export
#' 
slaterStandardization <- function(x)
{
  E <- dist(t(x), diag=TRUE, upper=TRUE)  # euclidean distance
  E <- as.matrix(E)
  m <- ncol(x)                            # number of elements  
  D <- sweep(x, 1, apply(x, 1, mean))     # row-center data
  S <- sum(diag(t(D) %*% D))
  U <- (2 * S/(m - 1))^0.5
  E/U                                     # devide by expected distance unit
}


#' Calculate Slater distance.
#'
#' The euclidean distance is often used as a measure of similarity between
#' elements (see \code{\link{distance}}. A drawback of this measure is that it 
#' depends on the range of the rating scale and the number of constructs used,
#' i. e. on the size of a grid. \cr 
#' An approach to standardize the euclidean distance to make it independent from
#' size and range of ratings and was proposed by Slater (1977, pp. 94). The
#' 'Slater distance' is the Euclidean distance divided by the expected distance.
#' Slater distances bigger than 1 are greater than expected, lesser than 1 are
#' smaller than expected. The minimum value is 0 and values bigger than 2 are
#' rarely found. Slater distances have been be used to compare inter-element
#' distances between different grids, where the grids do not need to have the
#' same constructs or elements. Hartmann (1992) showed that Slater distance is
#' not independent of grid size. Also the distribution of the Slater distances
#' is asymmetric. Hence, the upper and lower limit to infer 'significance' of
#' distance is not symmetric. The practical relevance of Hartmann's findings
#' have been demonstrated by Schoeneich and Klapp (1998). To calculate
#' Hartmann's version of the standardized distances see
#' \code{\link{distanceHartmann}}
#' 
#' 
#' @section Calculation: The Slater distance is calculated as follows. 
#' For a derivation see Slater (1977, p.94).       \cr
#' Let matrix \eqn{D}{D} contain the row centered ratings. Then
#'    \deqn{P = D^TD}{P = D^TD} and
#'    \deqn{S = tr(P)}{S = tr(P)}
#' The expected 'unit of expected distance' results as     \cr
#' \deqn{U = (2S/(m-1))^{1/2}}{U = (2S/(m-1))^.5} 
#' where \eqn{m}{m} denotes the number of elements of the grid.
#' The standardized Slater distances is the matrix of Euclidean distances
#' \eqn{E}{E} devided by the expected distance \eqn{U}{U}. 
#' \deqn{E/U}{E/U}
#'
#' 
#' @title 'Slater distances' (standardized Euclidean distances).
#'
#' @inheritParams     distance
#' @return            A matrix with Slater distances.
#'
#' @references        Hartmann, A. (1992). Element comparisons in repertory 
#'                    grid technique: Results and consequences of a Monte 
#'                    Carlo study. \emph{International Journal of Personal 
#'                    Construct Psychology, 5}(1), 41-56.
#'
#'                    Schoeneich, F., & Klapp, B. F. (1998). Standardization
#'                    of interelement distances in repertory grid technique
#'                    and its consequences for psychological interpretation 
#'                    of self-identity plots: An empirical study. 
#'                    \emph{Journal of Constructivist Psychology, 11}(1), 49-58.
#'
#'                    Slater, P. (1977). \emph{The measurement of intrapersonal 
#'                    space by Grid technique.} Vol. II. London: Wiley.
#'
#' @author            Mark Heckmann
#' @export
#' @seealso \code{\link{distanceHartmann}}
#' @examples 
#'
#'    distanceSlater(bell2010)
#'    distanceSlater(bell2010, trim=40)
#'
#'    d <- distanceSlater(bell2010)
#'    print(d)
#'    print(d, digits=4)
#'    
#'    # using Norris and Makhlouf-Norris (problematic) cutoffs
#'    print(d, cutoffs=c(.8, 1.2))
#'
distanceSlater <- function(x, trim=20, index=TRUE) {
  if (!inherits(x, "repgrid")) 
    stop("Object must be of class 'repgrid'")  
  E <- distance(x, along=2, index=index)
  E.sl <- slaterStandardization(E)
  notes <- c("\nNote that Slater distances cannot be compared across grids",
             "with a different number of constructs (see Hartmann, 1992).\n")
  attr(E.sl, "arguments") <- list(along=2, dmethod="Slater (standardized Euclidean)", 
                                  p=2, notes=notes)                                
  class(E.sl) <- c("distance", "matrix")
  E.sl
}


#### Hartmann distance ####

# helper functions for Slater distribution simulation
#
generate_quasi <- function(nc=5, ne=10, r=1:5, prob= rep(1, length(r))) 
{
  matrix(sample(r, size=nc*ne, replace=T, prob=prob), ncol=ne)
}


generate_quasis <- function(n, nc=5, ne=10, r=1:5, prob= rep(1, length(r)))
{
  replicate(n, generate_quasi(nc=nc, ne=ne, r=r, prob=prob), simplify=FALSE)
}


get_upper_triangle <- function(x)
{
  x[upper.tri(x, diag=FALSE)]
}


quasiDistributionDistanceSlater <- function(reps, nc, ne, range, 
                                            prob=NULL, progress=TRUE)
{
  q <- generate_quasis(reps, nc=nc, ne=ne, r=range, prob= NULL)
  fun <- lapply
  if (progress)
    fun <- lapply_pb 
  dist.sl <- fun(q, slaterStandardization)
  dist.sl <- lapply(dist.sl, get_upper_triangle)
  unlist(dist.sl)
}


# quasiDistributionDistanceSlater <- function(rep, nc, ne, range, prob=NULL, progress=TRUE)
# {
#   quasis <- randomGrids(rep, nc=nc, ne=ne, range=range, prob=prob, options=0)
#   if (progress)                 # whether to show progress bar
#     lapply_fun <- lapply_pb else
#       lapply_fun <- lapply
#   quasis.sd <- lapply_fun(quasis, function(x){
#     ds <- distanceSlater(x)
#     ds[lower.tri(ds, diag=FALSE)]
#   })
#   unlist(quasis.sd)         # return as vector
# }



# Return a list with the mean and sd as indicated in Hartmann's (1992) paper.
# 
getSlaterPaperPars <- function(nc)
{
  constructs <- NULL    # dummy to avoid R CMD CHECK non-visible variable binding NOTE
  
  # hartmann only provides values for a number of constructs between 7 and 21.
  # Smaller and bigger grids use the parameters with the next best number of
  # constructs from Hartmann
  if (nc < 7)
    nc <- 7
  if (nc > 21)
    nc <- 21
  
  ## constants ##
  # parameters for Slater distance distributions used to calculated Hartmann
  # distances as supplied by Hartmann, 1992,  p. 51. (only defined for 7 to 21
  # constructs)
  #
  hartmann.pars <- data.frame(constructs=7:21,
    mean=c(.97596, .97902, 0.98236, .98322, .98470, 
           .98643, .98749, .98811, .98908, .98972, 
           .99034, .99092, .99135, .99193, .99228),
    sd=c(.21910, .20376, .19211, .18240, .17396, .16416, .15860, .15374, .14700, 
         .14303, .13832, .13444, .13082, .12676, .12365))
  
  # "Therefore the means of all Z-transformed percentiles [avering three scale
  # range, MH] are suggested to be used as cutoffs for distance interpretation.
  # The use of the 5th and the 95th percentiles is recommended (see Table 5)"
  # (Hartmann, 1992, p.52).
  # 
  hartmann.cutoffs <- c(p01=2.492, p05=1.777, p10=1.387, 
                        p90=-1.186, p95=- 1.519, p99=- 2.129)
  slater.mean <- unlist(subset(hartmann.pars, constructs == nc, mean))
  slater.sd <- unlist(subset(hartmann.pars, constructs == nc, sd))
  list(mean=slater.mean, sd=slater.sd)  
}


simulateSlaterAndHartmannDistribution <- function(reps=1000, nc, ne, range, 
                                                  prob=NULL,
                                                  progress=TRUE)
{
  # this operation takes some time
  slater.vals <- quasiDistributionDistanceSlater(reps=reps, nc=nc, 
                                                 ne=ne, range=range,
                                                 prob=prob,                                      
                                                 progress=progress)  
  # mean and sd of Slater distribution
  mean.slater <- mean(slater.vals, na.rm=TRUE)   
  sd.slater <- sd(slater.vals, na.rm=TRUE)
  
  # conversion to Hartmann values
  hartmann.vals <- -1 * (slater.vals - mean.slater) / sd.slater  
  list(slater = slater.vals, 
       hartmann = hartmann.vals)
}


# NOT USED
# # caclulate coverage probability_ 
# coverageProbability <- function(x, cutoffs)
# {   
#   l  <- length(x)
#   oneCoverProb <- function(cutoff)
#     sum(x < cutoff) / l
#   sapply(cutoffs, oneCoverProb)
# }


getDistributionParameters <- function(x, probs=c(.01, .025, .05, .1, .9, .95, .975, .99), 
                                      na.rm=TRUE)
{
  pars <- describe(x)
  qs <- quantile(x, probs = probs, na.rm = na.rm)
  #cover.probs <- coverageProbability(x, cutoffs)       # get coverage probabalities for cutoffs
  list(pars=pars, quantiles=qs)
}


#' Calculate Hartmann distance
#' 
#' Hartmann (1992) showed in a simulation study that Slater distances (see 
#' \code{\link{distanceSlater}}) based on random grids, for which Slater coined 
#' the expression quasis, have a skewed distribution, a mean and a standard 
#' deviation depending on the number of constructs elicited. He suggested a 
#' linear transformation (z-transformation) which takes into account the 
#' estimated (or expected) mean and the standard deviation of the derived 
#' distribution to standardize Slater distance scores across different grid 
#' sizes. 'Hartmann distances' represent a more accurate version of 'Slater 
#' distances'. Note that Hartmann distances are multiplied by -1. Hence, 
#' negative Hartmann values represent dissimilarity, i.e. a big Slater distance.
#' \cr
#' 
#' There are two ways to use this function. Hartmann distances can either be 
#' calculated based on the reference values (i.e. means and standard deviations 
#' of Slater distance distributions) as given by Hartmann in his paper. The
#' second option is to conduct an instant simulation for the supplied grid 
#' size for each calculation. The second option will be more accurate when
#' a big number of quasis is used in the simulation. \cr
#' 
#' It is also possible to return the quantiles of the sample distribution and
#' only the element distances consideres 'significant' according to the
#' quantiles defined.
#'
#'
#' @section Calculation:
#' 
#' The 'Hartmann distance' is calculated as follows (Hartmann 1992, p. 49).  \cr
#' \deqn{D = -1 (\frac{D_{slater} - M_c}{sd_c})}{D = -1 (D_slater - M_c / sd_c)}
#' Where \eqn{D_{slater}}{D_slater} denotes the Slater distances of the grid,
#' \eqn{M_c}{M_c} the sample distribution's mean value and 
#' \eqn{sd_c}{sd_c} the sample distributions's standard deviation.
#'
#' @title 'Hartmann distance' (standardized Slater distances).
#' @param x           \code{repgrid} object.
#' @param method      The method used for distance calculation, on of 
#'                    \code{"paper", "simulate", "new"}. \code{"paper"} uses the 
#'                    parameters as given in Hartmann (1992) for caclulation.
#'                    \code{"simulate"} (default) simulates a Slater distribution
#'                    for the caclulation. In a future version the time consuming
#'                    simulation will be replaced by more accurate parameters for
#'                    Hartmann distances than used in Hartmann (1992).    
#' @param reps        Number of random grids to generate sample distribution for 
#'                    Slater distances (default is \code{10000}). Note that
#'                    a lot of samples may take a while to calculate.
#' @param prob        The probability of each rating value to occur. 
#'                    If \code{NULL} (default) the distribution is uniform.
#'                    The number of values must match the length of the rating scale.
#' @param progress    Whether to show a progress bar during simulation
#'                    (default is \code{TRUE}) (for \code{method="simulate"}).
#'                    May be useful when the distribution is estimated on the basis
#'                    of many quasis.
#' @param distributions Wether to additionally return the values of the simulated
#'                    distributions (Slater etc.) The default is \code{FALSE}
#'                    as it will quickly boost the object size.
#' @return            A matrix containing Hartmann distances. \cr
#'                    In the attributes several additional parameters can be found: \cr
#'                    \item{\code{"arguments"}}{A list of several parameters 
#'                           including \code{mean} and \code{sd} of Slater distribution.}
#'                    \item{\code{"quantiles"}}{Quantiles for Slater and Hartmann 
#'                          distance distribition.}
#'                    \item{\code{"distributions"}}{List with values of the 
#'                          simulated distributions.}
#'                                    
#' @references        Hartmann, A. (1992). Element comparisons in repertory 
#'                    grid technique: Results and consequences of a Monte 
#'                    Carlo study. \emph{International Journal of Personal 
#'                    Construct Psychology, 5}(1), 41-56.
#' @export
#' @author            Mark Heckmann
#' @seealso \code{\link{distanceSlater}}
#' @examples \dontrun{
#'
#'    ### basics  ###
#'    
#'    distanceHartmann(bell2010)
#'    distanceHartmann(bell2010, method="simulate")
#'    h <- distanceHartmann(bell2010, method="simulate")
#'    h
#'    
#'    # printing options
#'    print(h)
#'    print(h, digits=6)
#'    # 'significant' distances only
#'    print(h, p=c(.05, .95))
#'
#'    # access cells of distance matrix
#'    h[1,2]
#'    
#'    ### advanced  ###
#'    
#'    # histogram of Slater distances and indifference region
#'    h <- distanceHartmann(bell2010, distributions=TRUE)
#'    l <- attr(h, "distributions") 
#'    hist(l$slater, breaks=100)
#'    hist(l$hartmann, breaks=100)
#' }
#' 
distanceHartmann <- function(x, method="paper", reps=10000,  
                             prob=NULL, progress=TRUE, distributions=FALSE)
{
  if (distributions == TRUE & method != "simulate") {
    method <- "simulate"
    warning("'method' set to 'simulate' to return distributions", call.=FALSE)
  }
  if (!inherits(x, "repgrid")) 
    stop("Object must be of class 'repgrid'")
  ps <- seq(0, 1, .001)   # probabilty for quantiles to return
  
  # select parameter derivation for transformation
  method <- match.arg(method, c("paper", "simulate", "new"))

  ## get grid parameters ## 
  mm <- getScale(x)            # get min and max scale value
  range <- mm[1]:mm[2]
  nc <- getNoOfConstructs(x)
  ne <- getNoOfElements(x)
   
  # derive parameters mean and sd by simulation of Slater distance distributions
  # for given grid size
  if (method == "simulate") {
    v <- simulateSlaterAndHartmannDistribution(reps=reps, nc=nc, ne=ne, range=range, 
                                               prob=prob, progress=progress)
    sl <- getDistributionParameters(v$slater)
    hm <- getDistributionParameters(v$hartmann)    
  }
  
  # use parameters mean and sd from Hartmann paper or simulated 
  # TODO: use my own simulated parameters with bigger sample size
  #       for more accuracy than Hartmann. This will give accurate results
  #       while still avoiding time-consuming simulations.
  notes <- NULL
  if (method == "paper") {
    p <- getSlaterPaperPars(nc) 
    notes <- c("\nFor calculation the parameters from Hartmann (1992) were used.",
               "Use 'method=new' or method='simulate' for a more accurate version.\n")
  } else 
  if (method == "simulate")
    p <- list(mean=sl$pars$mean, sd=sl$pars$sd) else
  if (method == "new")
    stop("method 'new' has not yet been implemented", call.=FALSE) else
  stop("'method' must be 'paper', 'simulate' or 'new'", call.=FALSE)    
    
  # linear transformation to derive Hartmann distance from Slater distances 
  # (Hartmman, 1992, p. 49)
  D <- distanceSlater(x)
  H <- -1 * (D - p$mean) / p$sd
  #diag(H) <- 0      # replace diagonal as zeros (Hmmm?)
  
  # prepare output
  attr(H, "arguments") <- list(along=2, dmethod="Hartmann (standardized Slater distances)", p=2,
                               notes=notes,
                               parameters=p,
                               method=method)
  if (method == "simulate") {
    quantiles.slater <- quantile(v$slater, probs=ps)
    quantiles.hartmann <- quantile(v$hartmann, probs=ps)  
  } else {
    quantiles.slater <- NULL
    quantiles.hartmann <- NULL
  }
 
  attr(H, "quantiles") <- list(slater=quantiles.slater,
                               hartmann=quantiles.hartmann)
  # return Slater and Hartmann simulated distribution values if requested
  # caution: objects get quite big ~ 6mb with 10000 reps
  if (distributions){
    attr(H, "distributions") <- v       
  }  
  class(H) <- c("hdistance", "distance", "matrix")
  H
}


#' Print method for class hdistance (Hartmann distance objects).
#' 
#' Additionally allows the specification of p-values. The corresponding 
#' quantiles are calculated.  These are contained as attributes in the
#' \code{hdistance} object as returned by distanceHartmann. The lowest and
#' highest values are used as cutoffs in \code{print.distance}.
#' 
#' @param x           Object of class hdistance. 
#' @inheritParams      print.distance
#' @param p           Quantiles corresponding to probablities are used as cutoffs. 
#'                    Currently only works for Hartmann distances. If used 
#'                    \code{cutoffs} is overwritten.
#' @param ...         Not evaluated.                  
#' 
#' @export
#' @method            print hdistance
#' @keywords          internal
#' 
print.hdistance <- function(x, digits=2, col.index=TRUE,
                            upper=TRUE, diag=FALSE, cutoffs=NA, 
                            p=NA, ...)
{ 
  # only calculate quantiles when they are supplied in the object.
  # this is currently only the case for method="simulate". Will change
  # when bigger simulations are finished.
  do.quantiles <- FALSE
  if (attr(x, "arguments")$method == "simulate")
    do.quantiles <- TRUE
  
  # select quantiles to use (Hartmann or Power-transformed Hartmann distances)
  dmethod <- attr(x, "arguments")$dmethod 
  if (dmethod == "Hartmann (standardized Slater distances)")
    dm <- "hartmann" else dm <- "bc"
  
  # calculate quantiles (cutoff values) from lowest and biggest p value. get
  # quantiles from distance object (only for Hartmann distance). not the best 
  # approach, maybe change this in the future to use a lookup table when solid
  # cutoffs are available.
  if (!is.na(p[1]) & do.quantiles) { 
    which.p <- round(seq(0, 1, .001), 6) %in% round(p, 6)  # to avoid floating number representation inequalities
    qs <- attr(x, "quantiles")[[dm]][which.p]   
    cutoffs <- qs
  }
  
  # call standard printing function for distance objects
  print.distance(x=x, digits=digits, col.index=col.index, upper=upper,
                 diag=diag, cutoffs=cutoffs)
  
  # add quantiles used as cutoffs
  if (!is.na(p[1]) & do.quantiles) {
    cat("Quantiles:\n")
    print(round(qs, 4))  
    cat("\nThe smallest and biggest quantiles are used as cutoffs.\n")
  }
}


#### Power transformed Hartmann distance ####


#' Calculate power-transformed Hartmann distances.
#'
#' Hartmann (1992) suggested a transformation of Slater (1977) distances to make
#' them independent from the size of a grid. Hartmann distances are supposed to
#' yield stable cutoff values used to determine 'significance' of inter-element 
#' distances. It can be shown that Hartmann distances are still affected by grid
#' parameters like size and the range of the rating scale used (Heckmann, 2012).
#' The function \code{distanceNormalize} applies a Box-Cox (1964) transformation to the
#' Hartmann distances in order to remove the skew of the Hartmann distance
#' distribution. The normalized values show to have more stable cutoffs
#' (quantiles) and better properties for comparison across grids of different
#' size and scale range.  \cr
#' 
#' The function \code{distanceNormalize} can also return
#' the quantiles of the sample distribution and only the element distances
#' consideres 'significant' according to the quantiles defined.
#' 
#' @section Calculations:
#' The 'power tranformed Hartmann distance' are calulated as
#' follows: The simulated Hartmann distribution is added a constant as the
#' Box-Cox transformation can only be applied to positive values. Then a range
#' of values for lambda in the Box-Cox transformation (Box & Cox, 1964) are
#' tried out. The best lambda is the one maximizing the correlation of the
#' quantiles with the standard normal distribution. The lambda value maximizing
#' normality is used to transform Hartmann distances. As the resulting scale of
#' the power transformation depends on lambda, the resulting values are
#' z-transformed to derive a common scaling.
#' 
#' The code for the calculation of the optimal lambda was written by Ioannis 
#' Kosmidis.
#' 
#' @title             Standardized inter-element distances' (power transformed 
#'                    Hartmann distances).
#' @param x           \code{repgrid} object.
#' @param reps        Number of random grids to generate to produce
#'                    sample distribution for Hartmann distances
#'                    (default is \code{1000}). Note that
#'                    a lot of samples may take a while to calculate. 
#' @inheritParams     distanceHartmann
#'
#' @return            A matrix containing the standardized distances. \cr
#'                    Further data is contained in the object's attributes: \cr             
#'                    \item{\code{"arguments"}}{A list of several parameters 
#'                           including \code{mean} and \code{sd} of Slater distribution.}
#'                    \item{\code{"quantiles"}}{Quantiles for Slater, Hartmann 
#'                          and power transformed distance distribitions.}
#'                    \item{\code{"distributions"}}{List with values of the 
#'                          simulated distributions, if \code{distributions=TRUE}.}
#'                                    
#' @references        Box, G. E. P., & Cox, D. R. (1964). An Analysis of Transformations. 
#'                    \emph{Journal of the Royal Statistical Society. 
#'                    Series B (Methodological), 26}(2), 211-252.
#'
#'                    Hartmann, A. (1992). Element comparisons in repertory 
#'                    grid technique: Results and consequences of a Monte 
#'                    Carlo study. \emph{International Journal of Personal 
#'                    Construct Psychology, 5}(1), 41-56.
#'                    
#'                    Heckmann, M. (2012). Standardizing inter-element distances 
#'                    in grids - A revision of Hartmann's distances, 11th Biennal
#'                    Conference of the European Personal Construct Association 
#'                    (EPCA), Dublin, Ireland, Paper presentation, July 2012.
#' 
#'                    Slater, P. (1977). \emph{The measurement of intrapersonal space 
#'                    by Grid technique}. London: Wiley.
#'
#'
#' @export
#' @author            Mark Heckmann
#' @seealso           \code{\link{distanceHartmann}} and \code{\link{distanceSlater}}.
#' @examples \dontrun{
#'
#'    ### basics  ###
#'    
#'    distanceNormalized(bell2010)
#'    n <- distanceNormalized(bell2010)
#'    n
#'    
#'    # printing options
#'    print(n)
#'    print(n, digits=4)
#'    # 'significant' distances only
#'    print(n, p=c(.05, .95))
#'
#'    # access cells of distance matrix
#'    n[1,2]
#'    
#'    ### advanced  ###
#'    
#'    # histogram of Slater distances and indifference region
#'    n <- distanceNormalized(bell2010, distributions=TRUE)
#'    l <- attr(n, "distributions") 
#'    hist(l$bc, breaks=100)
#'  
#' }
#'
distanceNormalized <- function(x, reps=1000, prob=NULL, progress=TRUE,
                               distributions=TRUE)
{
  if (!inherits(x, "repgrid")) 
    stop("Object must be of class 'repgrid'")
  
  ps <- seq(0, 1, .001)   # probabilty for quantiles to return
  
  # calculate Hartmann and Slater distances
  h <- distanceHartmann(x, reps=reps, prob=prob, progress=progress, 
                        method="simulate", distributions=distributions)
  
  # optimal lambda for Box-Cox transformation. Add constant as only defined 
  # for positive values. Use offest 1 for same treatment of tails.
  d <- attr(h, "distributions")
  constant <- abs(min(c(d$hartmann, h))) + 1.00001  
  bc <- optimal.boxcox(d$hartmann + constant)  
  
  # parameters to standardize power transformed Hartmann values
  lambda.max <- bc$lambda
  sd.bc <- sd(bc$x)
  mean.bc <- mean(bc$x)
  
  # function to perform Box-Cox tranformation plus standardization
  bc.transform <- function(x){ # , constant, lambda.max, sd.bc, mean.bc){
    res <- ((x + constant)^lambda.max - 1) / lambda.max   # power transformation
    (res-mean.bc) / (sd.bc)           # z-transformation
  }
  
  # make bc transformations for all Hartmann data
  bc.dist <-  bc.transform(d$hartmann)    # transform simulated Hartmann distribition
  bc.vals <- bc.transform(h)         # transform grid data
  bc.qs <- quantile(bc.dist, ps, na.rm=TRUE)
   
  # prepare output
  notes <- NULL
  l <- list(dmethod="Power transformed Hartmann distances", notes=notes)
  attr(h, "arguments") <- modifyList(attr(h, "arguments"), l)
                                     
  # add distribution of power-transformed values
  if (distributions)
    attr(h, "distributions")$bc <- bc.dist
  # add quantiles of power-transformed values
  attr(h, "quantiles")$bc <- bc.qs
  h
}
  
  



