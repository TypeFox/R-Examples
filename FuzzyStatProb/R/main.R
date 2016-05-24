##' Fuzzy stationary probabilities of Markov chains from observations
##'
##' Computation of LR fuzzy numbers representing fuzzy stationary probabilities of an unknown Markov chain from which a sequence of observations has been drawn.
##' The fuzzy Markov chain considered during the processing follows the approach proposed by J. Buckley (see the reference section).
##' @param data This argument can be: (a) an array of either strings or natural numbers representing the observed states of the chain at consecutive time points. 
##'   The function first coerces the elements to a factor integer. (b) a 2D square matrix of strings representing fuzzy transition probabilities directly given by the user.
##'   Each string should be contained in \code{names(fuzzynumbers)} and refers to the corresponding \code{FuzzyNumber} object in the \code{fuzzynumbers} vector.
##'   When the transition probability from state i to j is 0 (in the crisp sense), then entry (i,j) must be NA. The \code{colnames} and \code{rownames}
##'   of the \code{data} matrix should have been set before calling this function.
##' @param options A tagged list containing the following parameters:
##'  \itemize{
##'    \item \code{verbose}: boolean, set to TRUE if progress information should be printed during the process. It is set to FALSE if this option is not specified. 
##'    \item \code{states}: an array of strings indicating the states for which the stationary distribution should be computed. 
##'			The values should match those specified in the \code{data} argument.
##'      If this option is not specified, the fuzzy stationary probabilities are computed for every state of the chain.
##'    \item \code{regression}: a string with the type of the regression to be applied at the end of the algorithm for fitting the membership functions of the fuzzy stationary 
##'			probabilities. 
##'      Possible values are \sQuote{linear}, \sQuote{quadratic}, \sQuote{cubic}, \sQuote{gaussian}, \sQuote{spline} and \sQuote{piecewise} (piecewise linear interpolation). 
##'      In all cases (including the gaussian), a different curve is fitted for each side of tue fuzzy number.
##'      The \code{gaussian} option fits curves of the form \eqn{\mu(x)} = exp \eqn{( -1/2 |(x-c)/s|^m)}.
##'      The \code{spline} option performs interpolation by a monotone cubic spline according to the Hyman method (see \code{splinefun} documentation) while \code{piecewise} computes a
##'      piecewise linear membership function by connecting consecutive points of the \eqn{\alpha}-cuts with straight lines, using the built-in \code{PiecewiseLinearFuzzyNumber} 
##'			subclass of the 
##'      \pkg{FuzzyNumbers} package. If this option is not specified, quadratic regression is carried out by default.
##'    \item \code{acutsonly}: boolean, set to TRUE if no regression should be done after computing the \eqn{\alpha}-cuts. This option is set to FALSE if not specified.
##'    \item \code{ncores}: positive integer representing the maximum number of cores that can be used when running in parallel. If set to more than 1, 
##'      then each processor takes care of all the computations involving one of the values of \eqn{\alpha} that have to be sampled, via 
##'			 \code{\link{parallel}} package. Defaults to 1 (sequential) if not specified. If \code{ncores} is greater than the actual number of cores in the computer,
##'      all available cores are used.
##'    \item{fuzzynumbers}: a tagged list with all the different FuzzyNumber objects that appear in \code{data} when \code{data} is a matrix of labels; 
##'     ignored otherwise. Every element of the list must have a name, referenced in at least one entry of \code{data}.
##' }
##' @param step Step size for sampling \eqn{\alpha} when computing the \eqn{\alpha}-cuts. The smallest \eqn{alpha} that is always present equals 0.001, 
##' and the rest of values are calculated as \eqn{\alpha = k } \code{step} for \eqn{k \ge 1}. 
##' The greatest sampled value that is always present as well is \eqn{\alpha = 0.999}.
##' Defaults to 0.05 when not specified.
##' @param ... Further arguments to be passed to \code{\link{DEoptim.control}} to customize the algorithm that finds the lower and upper bounds of the \eqn{\alpha}-cuts
##' by solving a minimization and a maximization problem.
##' @return An object of the new S3 class \code{FuzzyStatObj}, which is a tagged list with the following components:
##'   \item{fuzzyStatProb}{A list of \code{FuzzyNumber} objects. 
##'   The length of the list equals that of the \code{states} tag of the \code{options} argument. The object at a given position \code{i} corresponds to the 
##'   fuzzy stationary probability of the state indicated at position i of the \code{states} vector.
##'   If any of the states indicated in the \code{states} option is not found in the \code{data} input vector, the corresponding position in \code{fuzzyStatProb} 
##'   will be NA. If the function was called with \code{acutsonly} set to TRUE, then the returned object will not have a \code{fuzzyStatProb} tag.
##'   }
##' \item{acuts}{
##'   A list of data frame objects containing the \eqn{\alpha}-cuts of every fuzzy stationary probability, represented as bidimensional points
##'   (lowerBound,\eqn{\alpha}) and (upperBound,\eqn{\alpha}) where \eqn{\tilde{\pi}(\alpha) = [lowerBound, upperBound]} is an \eqn{\alpha}-cut of the fuzzy number
##'   \eqn{\tilde{\pi}}. The length of the list also equals that of the \code{states} tag of the \code{options} argument. Again, object at position \code{i} corresponds
##'   to \eqn{\alpha}-cuts of the state indicated at position \code{i} of the \code{states} vector of the option list.
##'   If any of the states indicated in the \code{states} option is not found in the \code{data} input vector, the corresponding position in \code{acuts} will be NA.
##'  }
##' @details Given a sequence of consecutive observations of the state of the chain, a fuzzy transition matrix is constructed according to the approach proposed in
##' J. Buckley's \emph{Fuzzy Probabilities} book. Fuzzy transition probabilities are constructed as the superposition of intervals (\eqn{\alpha}-cuts), 
##' which in this case represent simultaneous confidence intervals for multinomial proportions, and are computed using the input sequence of observations 
##' drawn from the chain. For each value of \eqn{\alpha}, the \eqn{\alpha}-cuts of such fuzzy transition probabilities define a matrix space
##' where we seek for the the matrices producing respectively the minimum and maximum possible stationary probability for each state of the chain,
##' using heuristic optimization tools (Differential Evolution). 
##' Both points define a closed real interval that is indeed an \eqn{\alpha} cut of the output fuzzy number representing the fuzzy stationary probability for that state.
##' Solving these problems for different \eqn{\alpha} allows to reconstruct the fuzzy stationary probabilities from their \eqn{\alpha}-cuts, 
##' applying the decomposition theorem. Regression is applied at the final stage to compute the membership functions of the stationary probabilities.
##' @references Buckley, J.J. Fuzzy Probabilities: New Approach and Applications, 2nd edition, volume
##'   115 of Studies in Fuzziness and Soft Computing. Springer, 2005.
##' @references Glaz, J. and C.P. Sison. Simultaneous confidence intervals for multinomial proportions. 
##' Journal of Statistical Planning and Inference 82:251-262 (1999). 
##' @references May, W.L. and W.D. Johnson. Constructing two-sided simultaneous confidence intervals for 
##'   multinomial proportions for small counts in a large number of cells. Journal of Statistical Software 5(6) (2000).
##'   Paper and code available at \url{http://www.jstatsoft.org/v05/i06}.
##' @references Gagolewski M. FuzzyNumbers Package: Tools to deal with fuzzy numbers in R (2012).
##' Tutorial available at http://www.ibspan.waw.pl/~gagolews/FuzzyNumbers/doc/FuzzyNumbers-Tutorial.pdf
##' @examples
##' # ----------------- CREATE DATA ----------
##' # Simulate 200 observations of a 10-state Markov chain, 
##' # and compute fuzzy stationary probability of state 1
##' if(require("markovchain")){ # for simulating from a known crisp Markov chain
##' 	# Transition matrix taken from Fig. 1 of Amigoni et al. (see references)
##' 	mcPatrol <- new("markovchain", states = robotStates, byrow = TRUE,
##' 	transitionMatrix = transRobot, name = "Patrolling")
##' 	set.seed(666)
##' 	simulatedData <- rmarkovchain(n = 200, object = mcPatrol, t0 = 
##'   sample(robotStates, 1))
##' 	mcfit = markovchainFit(simulatedData) # Fit with markovchain package
##' 	vsteady = steadyStates(mcfit$estimate) # 1 x n matrix of stat. probs
##' 	# ---------------------------------------
##' 	# Simplest case: compute only alpha-cuts for alpha=0.001 and alpha=0.999
##' 	# Set itermax to 30 (too few) just for a fast example (not good results)
##' 	linear = fuzzyStationaryProb(simulatedData,list(verbose=TRUE, states="01", 
##'   	regression="piecewise"), step=1, itermax = 30) 
##' 	summary(linear)
##' 	linear$fuzzyStatProb[["01"]]
##' 	plot(linear$fuzzyStatProb[["01"]])
##' 	points(linear$acuts[["01"]])
##'	}
##' \dontrun{
##' # A more accurate approximation, with steps of 0.1 (takes much longer!)
##' # Run the previous code to create mcPatrol, vsteady and simlatedData
##' quadratic = fuzzyStationaryProb(data = simulatedData,list(verbose=TRUE, 
##'   ncores = 2, regression="quadratic"), step=0.1)
##' m <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,11),nrow = 4,ncol = 3,byrow = TRUE)
##' layout(mat = m,heights = c(0.25,0.25,0.25,0.25))
##' for (i in robotStates){
##' par(mar = c(4,4,2,1))
##'     plot(quadratic$fuzzyStatProb[[i]],col="red",main=paste("State ",i), 
##'       cex.lab = 1.1,lwd=2);    
##'     points(quadratic$acuts[[i]]);
##'     abline(v = vsteady[1,i], lty = "dashed");
##' }
##' plot(1, type = "n", axes=FALSE, xlab="", ylab="")
##' plot_colors <- c("red")
##' legend(x = "top",inset = 0, legend = c("Quadratic"), col=plot_colors, 
##'   bty = "n", lwd=2, cex=1, horiz = FALSE)
##' 
##' # Now departing from user-specified fuzzy transition probabilities
##' library(FuzzyNumbers)
##' EU = TrapezoidalFuzzyNumber(0,0,0.02,0.07); # Extremely unlikely 
##' VLC = TrapezoidalFuzzyNumber(0.04,0.1,0.18,0.23); # Very low chance
##' SC = TrapezoidalFuzzyNumber(0.17,0.22,0.36,0.42); # Small chance
##' IM = TrapezoidalFuzzyNumber(0.32,0.41,0.58,0.65); # It may
##' MC = TrapezoidalFuzzyNumber(0.58,0.63,0.8,0.86); # Meaningful chance
##' ML = TrapezoidalFuzzyNumber(0.72,0.78,0.92,0.97); # Most likely
##' EL = TrapezoidalFuzzyNumber(0.93,0.98,1,1); # Extremely likely
##' allnumbers = c(EU,VLC,SC,IM,MC,ML,EL);
##' names(allnumbers) = c("EU","VLC","SC","IM","MC","ML","EL");
##' rownames(linguisticTransitions) = robotStates; # see the package data
##' colnames(linguisticTransitions) = robotStates;
##' 
##' # Simplest case: compute only alpha-cuts for alpha=0.001 and alpha=0.999
##' # linguisticTransitions is a matrix of strings defined in the package data
##' linear = fuzzyStationaryProb(linguisticTransitions,list(verbose=TRUE, 
##'   regression="linear", ncores = 4, fuzzynumbers = allnumbers),step=0.2)
##' summary(linear)
##' }
##' @references Amigoni, F., Basilico, N., Gatti, N. Finding the Optimal Strategies for Robotic Patrolling
##' with Adversaries in Topologically-Represented Eenvironments. In Proc. of ICRA 2009, pp. 819-824.
##' @seealso \code{\link[markovchain]{markovchainFit}}
fuzzyStationaryProb<-function(data,options,step=0.05, ...){

  t1 = proc.time();

## -----------------------------------------------------------------
##                          CHECK ARGUMENTS 
## -----------------------------------------------------------------

  r = names(options);
  present = r %in% c("verbose", "regression", "states", "acutsonly", "ncores", "fuzzynumbers");
	if(sum(!present) > 0){
		stop("ERROR: unrecognized option(s) : ", r[!present],"\n");
	}

	pointlistleft = NULL;
	pointlistright = NULL;
	originalstates = NULL;

## ---------------------------------------------------------------------
##   COMPUTE ALPHA-CUTS OF STATIONARY PROBABILITIES FOR DESIRED STATES 
## ---------------------------------------------------------------------
  temp = .stationaryAcuts(data, options, step, ...);
  pointlistleft = temp[[1]]; # lower bounds of the stationary alpha-cuts
  pointlistright = temp[[2]]; # upper bounds of the stationary alpha-cuts
  originalstates = temp$originalstates; # names of the states for which stationary alpha-cuts were computed
  nresults = temp$nresults; # number of states for which the stationary probabilities must be computed
  listconf = temp$listconf; # levels of alpha to be sampled, for which alpha-cuts will be computed
  options = temp$options;
  nstates = temp$nstates;
  ncores = temp$ncores;
  
  t2 = proc.time();
  if(options$verbose) cat("...finished successfully (elapsed: ",t2[3]-t1[3],"s.)\n");

## ---------------------------------------------------------------------
## USE REGRESSION TO OBTAIN THE MEMBERSHIP FUNCTIONS FITTING THE ALPHA-CUTS
## ---------------------------------------------------------------------
  
  if(options$verbose){ 
    cat("Applying",options$regression,"regression to obtain the membership functions of fuzzy stationary probabilities...\n");
    cat("Fitting memb.func. for state: "); flush.console();
  }
  # Now do the regression and finally build the fuzzy numbers defined by the alpha-cuts computed above
  fuzzylist = vector("list", nresults);
  myframes = vector("list", nresults);
  
  for(i in 1:nresults){    
    if(options$verbose){ cat(originalstates[i]," ");    flush.console(); }
    
    myframeleft = data.frame(pointlistleft[[i]],listconf);
    myframeright = data.frame(pointlistright[[i]],listconf);    
    
    if(options$acutsonly){
        names(myframeleft) = c("x","y");
        names(myframeright) = c("x","y");
        unframe = rbind(myframeleft, myframeright);
        myframes[[i]] = unframe;
    }
    else{
      linguistic = is.matrix(data);
      mylist = .membFnRegression(myframeleft, myframeright, options$regression, linguistic);

      if(options$verbose && i==nresults) cat("\n");    
      if(!is.null(mylist)){
        fuzzylist[[i]] = mylist[[1]];
        myframeleft = mylist[[2]];
        myframeright = mylist[[3]];
        names(myframeleft) = c("x","y");
        names(myframeright) = c("x","y");
        unframe = rbind(myframeleft, myframeright);
        myframes[[i]] = unframe;
      }else{
        fuzzylist[[i]] = NA;
        myframes[[i]] = NA;
      }
    }
  }
  
  names(myframes)  = originalstates;
  names(fuzzylist) = originalstates;
    
## ---------------------------------------------------------------------  
  
  t3 = proc.time();
  
  msg = .composeOutputMessage(nstates, originalstates, step, data, ncores, options, t1, t2, t3);
  
  if(options$acutsonly){ result = list(acuts=myframes, summary=msg); }
  else{ result = list(fuzzyStatProb=fuzzylist, acuts=myframes, summary=msg); }
  
  class(result)<-"FuzzyStatObj";
  return(result);
}

## _____________________________________________________________________________________________________________

.composeOutputMessage <- function(nstates, originalstates, step, data, ncores, options, t1, t2, t3){
  
  msg = paste(     ". Fuzzy stationary probabilities of a Markov chain with",nstates,"states\n");
  msg = paste(msg, ". Probabilities have been computed for states:", paste(originalstates, collapse=" "),"\n");
  if(!is.matrix(data)){  
    msg = paste(msg, ". Number of input observations:",length(data),"\n");}
  else{
    msg = paste(msg, ". Fuzzy transition probabilities specified by the user\n");
  }  
  msg = paste(msg, ". Parameters:\n");
  msg = paste(msg, "       Step size:",step,"\n");
  if(options$ncores > 1){
    msg = paste(msg, "       Execution was done in parallel (",ncores,"cores used )\n");
  }
  else{
    msg = paste(msg, "       Execution was done sequentially\n");
  }  
  if(!options$acutsonly){
    msg = paste(msg, "       Regression curve for output membership functions:",options$regression);  
  }
  else{
    msg = paste(msg, "       (a-cuts only; no membership function was computed)\n");
  }
  msg = paste(msg, ". To retrieve the results use $fuzzyStatProb and $acuts with this object\n");
  msg = paste(msg, ". Computation of alpha-cuts took", format(round(t2[3] - t1[3],2),nsmall=2),"seconds\n");
  if(!options$acutsonly) msg = paste(msg, ". Membership functions regression took", format(round(t3[3] - t2[3],2),nsmall=2),"seconds\n");
  
  return(msg);
}

## _____________________________________________________________________________________________________________

.stationaryAcuts <- function(data, options, step, ...){

  originalstates = NULL;
  nstates = NULL;
  
  if(!is.matrix(data)){
    levelsdata = NULL;
    if(!is.factor(data)){
      data = as.factor(data);
    }
    levelsdata = levels(data);
    data = as.integer(data);
    originalstates = options$states;
    
    if(!is.null(options$states)){
      options$states = match(options$states, levelsdata);
      if(sum(is.na(options$states)) > 0){
        stop(paste0("ERROR: states not found in the data: ",originalstates[is.na(options$states)],"\n"));
      }
    }
    else{
      options$states = 1:length(levelsdata);
      originalstates = levelsdata;
    }
    
    stopifnot(min(data) == 1);
    nstates = length(levelsdata);
  }
  else{
    # The user does not have observations data but has provided a matrix of 
    # fuzzy numbers objects directly with the fuzzy transition probabilities
    originalstates = colnames(data);
    if(is.null(originalstates)){ originalstates = rownames(data); }
    if(is.null(rownames(data)) && is.null(colnames(data))){
      stop("Both the row names and the column names of the fuzzy transition matrix are null.\nPlease give names to either the rows or the columns.");
    }
    if(is.null(options$states)){ options$states = originalstates;    }
    options$states = match(options$states, originalstates);
    nstates = length(originalstates);
    
    # Check feasibility of the fuzzy transition matrix: every row is feasible
    checkFuzzyProbabilityFeasibility(data, options$fuzzynumbers);
  }
  
  nresults = length(options$states);

  if(is.null(options$regression)){
    options$regression = "quadratic";
  }
  else{    
    options$regression = tolower(options$regression);
    reg = options$regression;
    if(!(reg == "linear" || reg == "quadratic" || reg == "cubic" || reg == "gaussian" || reg == "spline" || reg == "piecewise")){
      stop("ERROR: regression must be one of: linear | quadratic | cubic | gaussian | spline | piecewise\n");
    }
  }
  
  if(is.null(options$ncores)){
    options$ncores = 1;
  }  
  else{
      if(options$ncores <= 0){
        stop("ERROR: the value of options$ncores must be a positive integer");
      }
  }
  
  if(is.null(options$verbose)){
    options$verbose = FALSE;
  }
  
  if(is.null(options$acutsonly)){
    options$acutsonly = FALSE;    
  }
  else{
    if(!(options$acutsonly == TRUE || options$acutsonly == FALSE)){
      stop("ERROR: the value of options$acutsonly must be either TRUE or FALSE\n");
    }
  }
  

## ---------------------------------------------------------------------
##          COMPUTE ALPHA-CUTS FOR DESIRED STATES 
## ---------------------------------------------------------------------

  iterations = 1+1.0/step;  
  listaconf = rep(NA,iterations);  
  pointlistleft = vector("list", length = nresults);
  pointlistright = vector("list", length = nresults);

  for(j in 1:nresults){
    pointlistleft[[j]] = rep(NA,iterations);
    pointlistright[[j]] = rep(NA,iterations);
  }

  for(i in 1:iterations){
    if(i == 1){ confidence = 0.001;   }
    else{ 
      if(i == iterations){ confidence = 0.999; }
      else{ confidence = step * (i-1); }
    }    
    listaconf[i] = confidence;
  }
  ncores = 1;
 
  ## PARALLEL -------------------------------------------
  if(options$ncores > 1){
    ncores = min(detectCores(), options$ncores);
    cl = makeCluster(ncores);
    if(options$verbose){ cat("Parallel computation of a-cuts in progress"); flush.console(); }
    res = parLapply(cl, X = listaconf, .computeFuzzyStationary, data = data, d = nstates, states = options$states, fuzzynumbers = options$fuzzynumbers, ...);
    stopCluster(cl);
    for(k in 1:iterations){
      for(j in 1:nresults){
        pointlistleft[[j]][k] = res[[k]]$left[j];
        pointlistright[[j]][k] = res[[k]]$right[j];
      }
    }
  }
  ## SEQUENTIAL -----------------------------------------
  else{
    if(options$verbose){ cat("Computing a-cuts for a = "); flush.console();}
    for(i in 1:iterations){
      confidence = listaconf[i];
      if(options$verbose){ cat(confidence," "); flush.console(); }
      res = .computeFuzzyStationary(confidence,data,nstates,options$states, fuzzynumbers = options$fuzzynumbers, ...);
      for(j in 1:nresults){
        pointlistleft[[j]][i] = res$left[j];
        pointlistright[[j]][i] = res$right[j];
      }
    }
  }
  
  return(list(pointlistleft, pointlistright, originalstates = originalstates, nresults = nresults, 
         listconf = listaconf, options = options, nstates = nstates, ncores = ncores));
}
## _____________________________________________________________________________________________________________

checkFuzzyProbabilityFeasibility <- function(data, fuzzynumbers){

  data[is.na(data)] = "impossible.at.all";
  impossible = TrapezoidalFuzzyNumber(0,0,0,0);
  fuzzynumbers = c(fuzzynumbers, impossible);
  names(fuzzynumbers)[length(fuzzynumbers)] = "impossible.at.all";
  for(i in 1:nrow(data)){    
    ith.row = fuzzynumbers[data[i,]];
    cores = t(sapply(ith.row, core)); # matrix of dimensions |nstates| x 2
    limits = colSums(cores);
    if(limits[1] > 1){
      stop("Row ", i, " is not a well-formed prob.distrib: lower bound of the support of the sum cannot be greater than 1");
    }
    if(limits[2] < 1){
      stop("Row ", i, " is not a well-formed prob.distrib: upper bound of the support of the sum cannot be smaller than 1");
    }
  }
}

## _____________________________________________________________________________________________________________

# Returns a square matrix m such that m[i,j] is the number of times a direct transition from state i to j has been observed
.countTransitions<-function(x,nstates){
    result = matrix(0,nstates,nstates);
    n = length(x);
    for(i in 1:n-1){
      result[x[i],x[i+1]] = result[x[i],x[i+1]] + 1;
    }
    return(result);
}

## _____________________________________________________________________________________________________________

# Computes the matrix of intervals representing alpha-cuts of the fuzzy probabilities. Returns a matrix of
# intervals lower bounds, a matrix of intervals upper bounds, a list of row indices and a list of column indices (negative indicates last non-null element in a row)
# indicating the entries of the matrix where the probabilities are greater than 0, and a crisp punctual estimation of the transition probabilities according to observed data.
# Argument 'a' is the alpha-level required for the alpha-cuts
.computeIntervalMatrices<-function(x,nstates,a){
    
  countsMatrix = .countTransitions(x,nstates);
  lowerBoundsMatrix = matrix(0,nstates,nstates);
  upperBoundsMatrix = matrix(0,nstates,nstates);
  punctualEstimatesMatrix = matrix(0,nstates,nstates);
  
  listaposx = {}
  listaposy = {}
  
  iniciales = {};
  
  lower = {};
  upper = {};
    
  for(i in 1:nstates){
    if(i > 1){
      listaposy[length(listaposy)] = (-1)*listaposy[length(listaposy)];          
      # Remove the last variable of this probability distribution (it is determined as 1 - sum of the rest)
      iniciales = iniciales[-length(iniciales)];
      upper = upper[-length(upper)];
      lower = lower[-length(lower)];
    }
    
    suma = sum(countsMatrix[i,]);
    if(suma > 0){
      temp = multinomialCI(countsMatrix[i,],a,verbose=FALSE);
    }
    else{
      temp = matrix(0, nstates, 2);
    }
        
    # Traverse all the row
    for(j in 1:nstates){ 
      if(countsMatrix[i,j] > 0){
        lowerBoundsMatrix[i,j] = temp[j,1];
        upperBoundsMatrix[i,j] = temp[j,2];
        listaposx = c(listaposx,i);
        listaposy = c(listaposy,j);
        iniciales = c(iniciales, (lowerBoundsMatrix[i,j] + upperBoundsMatrix[i,j])/2); # IGNORADO
        if(lowerBoundsMatrix[i,j] == upperBoundsMatrix[i,j] && lowerBoundsMatrix[i,j] != 0){
          lower = c(lower, lowerBoundsMatrix[i,j]-0.01);
        }else{
          lower = c(lower, lowerBoundsMatrix[i,j]);
        }
        upper = c(upper, upperBoundsMatrix[i,j]);      
      }
    }

    if(suma > 0){      
      punctualEstimatesMatrix[i,] = countsMatrix[i,]/suma;
    }
    else{
      punctualEstimatesMatrix[i,] = rep(0, nstates);
    }
  }  
  # Last non-zero element of the last row
  listaposy[length(listaposy)] = (-1)*listaposy[length(listaposy)];    
  iniciales = iniciales[-length(iniciales)];
  upper = upper[-length(upper)];
  lower = lower[-length(lower)];
    
  # take as initial solution the observed transition matrix (frequencies) instead 
  # of the middle point of each interval
  antiguo = length(iniciales);
  iniciales = {};
  for(indice in 1:length(listaposy)){
    if(listaposy[indice] >= 0){
      iniciales = c(iniciales, punctualEstimatesMatrix[listaposx[indice],listaposy[indice]]);
    }
  }

  # Stationary probabilities using punctual estimation of the transition matrix
  identidad<-diag(nstates);
  ceros<-matrix(0,nstates,nstates);
  mat2=.unos(ceros) %*% solve(.unos(punctualEstimatesMatrix-identidad));
  estac=mat2[1,];
  
  return(list(left=lowerBoundsMatrix,right=upperBoundsMatrix,listaposx=listaposx,listaposy=listaposy,lower=lower,upper=upper,iniciales=iniciales,puntuales=estac));
}

## _____________________________________________________________________________________________________________

# Computes the matrix of intervals representing alpha-cuts of the fuzzy probabilities. Returns a matrix of
# intervals lower bounds, a matrix of intervals upper bounds, a list of row indices and a list of column indices (negative indicates last non-null element in a row)
# indicating the entries of the matrix where the probabilities are greater than 0, and a crisp punctual estimation of the transition probabilities according to observed data.
# Argument 'a' is the alpha-level required for the alpha-cuts. Argument 'x' is a matrix of labels (strings) that the names of the elements of the 'fuzzynumbers' argument. 
# Argument 'fuzzynumbers' is a tagged list of FuzzyNumber objects. 
.getIntervalMatricesFromLabels<-function(x,nstates,a, fuzzynumbers){
      
  lowerBoundsMatrix = matrix(0,nstates,nstates);
  upperBoundsMatrix = matrix(0,nstates,nstates);
  punctualEstimatesMatrix = matrix(0,nstates,nstates);
  
  impossible = TrapezoidalFuzzyNumber(0,0,0,0);
  fuzzynumbers = c(fuzzynumbers, impossible);
  names(fuzzynumbers)[length(fuzzynumbers)] = "impossible.at.all";
  
  x[is.na(x)] = "impossible.at.all";
  
  listaposx = {}
  listaposy = {}
  
  iniciales = {};
  
  lower = {};
  upper = {};
    
  for(i in 1:nstates){
    if(i > 1){
      listaposy[length(listaposy)] = (-1)*listaposy[length(listaposy)];          
      # Remove the last variable of this probability distribution (it is determined as 1 - sum of the rest)
      iniciales = iniciales[-length(iniciales)];
      upper = upper[-length(upper)];
      lower = lower[-length(lower)];
    }
    
    suma = sum(x[i,] != "impossible.at.all"); # NA means there is 0 chance of transition from one state to another
    if(suma > 0){
      ith.row = fuzzynumbers[x[i,]]; # ith.row is a tagged list of FuzzyNumber objects
      temp = t(sapply(ith.row, alphacut, a)); # get the alpha-cuts of all the FuzzyNumbers of this row
    }
    else{
      temp = matrix(0, nstates, 2);
    }
        
    # Traverse all the row
    for(j in 1:nstates){ 
      if(x[i,j] != "impossible.at.all"){
        lowerBoundsMatrix[i,j] = temp[j,1];
        upperBoundsMatrix[i,j] = temp[j,2];
        listaposx = c(listaposx,i);
        listaposy = c(listaposy,j);
        iniciales = c(iniciales, (lowerBoundsMatrix[i,j] + upperBoundsMatrix[i,j])/2); # IGNORADO
        if(lowerBoundsMatrix[i,j] == upperBoundsMatrix[i,j] && lowerBoundsMatrix[i,j] != 0){
          lower = c(lower, lowerBoundsMatrix[i,j]-0.01);
        }else{
          lower = c(lower, lowerBoundsMatrix[i,j]);
        }
        upper = c(upper, upperBoundsMatrix[i,j]);      
      }
    }

    if(suma > 0){
      core = sapply(ith.row, alphacut, 1.0);      
      punctualEstimatesMatrix[i,] = (core[1,]+core[2,])/2; # central points of the TrFNs, which do not form a prob.distrib.
    }
    else{
      punctualEstimatesMatrix[i,] = rep(0, nstates);
    }
  }  
  # Last non-zero element of the last row
  listaposy[length(listaposy)] = (-1)*listaposy[length(listaposy)];    
  iniciales = iniciales[-length(iniciales)];
  upper = upper[-length(upper)];
  lower = lower[-length(lower)];
    
  # take as initial solution the observed transition matrix (frequencies) instead 
  # of the middle point of each interval
  antiguo = length(iniciales);
  iniciales = {};
  for(indice in 1:length(listaposy)){
    if(listaposy[indice] >= 0){
      iniciales = c(iniciales, punctualEstimatesMatrix[listaposx[indice],listaposy[indice]]);
    }
  }

  # Stationary probabilities using punctual estimation of the transition matrix
  identidad<-diag(nstates);
  ceros<-matrix(0,nstates,nstates);
  mat2=.unos(ceros) %*% solve(.unos(punctualEstimatesMatrix-identidad));
  estac=mat2[1,];
  
  return(list(left=lowerBoundsMatrix,right=upperBoundsMatrix,listaposx=listaposx,listaposy=listaposy,lower=lower,upper=upper,iniciales=iniciales,puntuales=estac));
}

## _____________________________________________________________________________________________________________

## Computes the left and right alpha-cut bounds (with the specified alpha argument) of the stationary probabilities for the
## states indicated in states. Returns a list with the left bounds of the alpha-cut(s) and another list with the right bounds. 
## The size of both lists matches the length of argument states. Argument d is the total number of states of the chain.
## Argument ... is a tagged list passed to DEoptim.control for customizing the DE algorithm
.computeFuzzyStationary<-function(alfa,data, d, states, fuzzynumbers, ...){

  # Compute the interval matrices that are restrictions for the optimization process
  # The intervals are computed either (a) using multinomial confidence intervals if data is an observation vector,
  # or (b) directly taking the alpha-cuts of the user-specified fuzzy numbers objects of fuzzy transition probabilities
  res = NULL;
  if(is.matrix(data)){
    # Get the interval matrix simply taking the alpha-cuts of the fuzzy numbers provided by the user
    res=.getIntervalMatricesFromLabels(data,d,alfa, fuzzynumbers);
  }
  else{  
    # Compute the interval matrix using multinomial confidence intervals from the observations provided by the user
    res=.computeIntervalMatrices(data,d,alfa);
  }
  stationaryLeft = {};
  stationaryRight = {};
  
  if((alfa < 0.999 && !is.matrix(data)) || is.matrix(data)){ 
    # When dealing with Fuzzy Numbers provided by the user, we assume that the 0.999-cut always exists
    for(j in 1:length(states)){    
       state = states[j];
       pob1 = .generateInitialPopulation(length(res$lower), res$left, res$right);
       pob2 = .generateInitialPopulation(length(res$lower), res$left, res$right);
       
       mycontrol = list();
       mycontrol$trace=FALSE; mycontrol$itermax=200; mycontrol$CR=0.8; 
       mycontrol$reltol=1E-3; mycontrol$steptol=20;
       
       testVoid <- function(...){ if(length(list(...))){ FALSE } else { TRUE } }
 
       if(!testVoid(...)){
         replacement = list(...);
         mycontrol[names(replacement)] = replacement;
       }

       mycontrol$initialpop=pob1;
       control1 = do.call(DEoptim.control, mycontrol);
       
       mycontrol$initialpop=pob2;
       control2 = do.call(DEoptim.control, mycontrol);
       
       opt = DEoptim(.stationary, res$lower, res$upper, control = control1, d, res$listaposx, res$listaposy,component=state,res$left,res$right,FALSE);
       stationaryLeft = c(stationaryLeft, opt$optim$bestval);
       opt = DEoptim(.stationary, res$lower, res$upper, control = control2, d, res$listaposx, res$listaposy,component=state,res$left,res$right,TRUE);
       stationaryRight = c(stationaryRight, (-1)*opt$optim$bestval);
    }
  }
  else{ # Punctual estimations for both lower and upper bounds (no need to launch DE)
    for(j in 1:length(states)){
      stationaryLeft = c(stationaryLeft, res$puntuales[states[j]]);
      stationaryRight = c(stationaryRight, res$puntuales[states[j]]);
    }
  }

  return(list(left=stationaryLeft, right=stationaryRight, puntuales=res$puntuales));
}

## _____________________________________________________________________________________________________________
