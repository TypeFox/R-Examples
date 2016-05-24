#' Calculate power values
#' 
#' Calculates local power values, expected number of rejections, the power to
#' reject at least one hypothesis and the power to reject all hypotheses.
#' 
#' 
#' @param x A matrix containing the rejected hypothesis, as produces by the
#' graphTest function.
#' @param f List of user defined power functions. If one is interested in the
#' power to reject hypotheses 1 and 3 one could specify \code{function(x) {x[1]
#' && x[3]}}. If f is a named list, the result will contain corresponding items 
#' with the same names (among the default elements described in the following).
#' @return A list containg at least the following four elements and
#' an element for each element in the parameter \code{f}.
#' \describe{
#' \item{\code{LocPower}}{A numeric giving the local powers for the hypotheses}
#' \item{\code{ExpNrRej}}{The expected number of rejections}
#' \item{\code{PowAtlst1}}{The power to reject at least one hypothesis}
#' \item{\code{RejectAll}}{The power to reject all hypotheses}
#' }
#' @keywords htest
#' @export extractPower
extractPower <- function(x, f=list()) {
  pow <- colMeans(x)
  avgPow <- sum(x)/nrow(x)
  atleast1 <- mean(rowSums(x)>0)
  allPow <- mean(rowSums(x)==dim(x)[2])
  result <- list(LocalPower = pow, ExpRejections = avgPow,
		  PowAtlst1 = atleast1, RejectAll = allPow)
  if (is.function(f)) {f <- list(f)}
  if (length(f)>0) {
    n <- names(f)
    if (is.null(n) || all(is.na(n))) n <- paste("func", 1:length(f), sep="")
    n[n=="" | is.na(n)] <- paste("func", 1:sum(n==""), sep="")
    names(f) <- n
  }
  for (fn in names(f)) {
    suppressWarnings(# We don't want to see "coercing argument of type 'double' to logical" warnings.
      result[[fn]] <- sum(apply(x,1, f[[fn]]))/dim(x)[1]
    )
  }
  result
}

#' Calculate power values
#' 
#' Given the distribution under the alternative (assumed to be multivariate
#' normal), this function calculates the power to reject at least one
#' hypothesis, the local power for the hypotheses as well as the expected
#' number of rejections.
#' 
#' @param weights Initial weight levels for the test procedure (see graphTest
#' function). Alternatively a \code{\link{graphMCP}} object can be given as parameter \code{graph}.
#' @param alpha Overall alpha level of the procedure, see graphTest function.
#' (For entangled graphs \code{alpha} should be a numeric vector of length 
#' equal to the number of graphs, each element specifying the partial alpha 
#' for the respective graph.
#' The overall alpha level equals \code{sum(alpha)}.)
#' @param G Matrix determining the graph underlying the test procedure. Note
#' that the diagonal need to contain only 0s, while the rows need to sum to 1.
#' When multiple graphs should be used this needs to be a list containing the
#' different graphs as elements. Alternatively a \code{\link{graphMCP}} object can be given as parameter \code{graph}.
#' @param graph A graph of class \code{\link{graphMCP}}.
#' @param mean Mean under the alternative
#' @param corr.sim Covariance matrix under the alternative.
#' @param corr.test Correlation matrix that should be used for the parametric test.
#' If \code{corr.test==NULL} the Bonferroni based test procedure is used. Can contain
#' NAs.
#' @param type What type of random numbers to use. \code{quasirandom} uses a
#' randomized Lattice rule, and should be more efficient than
#' \code{pseudorandom} that uses ordinary (pseudo) random numbers.
#' @param n.sim Monte Carlo sample size. If type = "quasirandom" this number is
#' rounded up to the next power of 2, e.g. 1000 is rounded up to
#' \eqn{1024=2^10}{1024=2^10} and at least 1024.
#' @param f List of user defined power functions (or just a single power
#' function).  If one is interested in the power to reject hypotheses 1 and 3
#' one could specify: \cr\code{f=function(x) {x[1] && x[3]}}.\cr If the power
#' of rejecting hypotheses 1 and 2 is also of interest one would use a
#' (optionally named) list: \cr 
#' \code{f=list(power1and3=function(x) {x[1] && x[3]},}\cr
#' \code{power1and2=function(x) {x[1] && x[2]})}.
#' If the list has no names, the functions will be referenced 
#' to as "func1", "func2", etc. in the output.
#' @param test In the parametric case there is more than one way to handle
#' subgraphs with less than the full alpha. If the parameter \code{test} is
#' missing, the tests are performed as described by Bretz et al. (2011), i.e.
#' tests of intersection null hypotheses always exhaust the full alpha level
#' even if the sum of weights is strictly smaller than one. If
#' \code{test="simple-parametric"} the tests are performed as defined in
#' Equation (3) of Bretz et al. (2011).
#' @param upscale Logical. If \code{upscale=FALSE} then for each intersection 
#' of hypotheses (i.e. each subgraph) a weighted test is performed at the 
#' possibly reduced level alpha of sum(w)*alpha, 
#' where sum(w) is the sum of all node weights in this subset.
#' If \code{upscale=TRUE} all weights are upscaled, so that sum(w)=1.
#' @param ... For backwards compatibility. For example up to version 0.8-7
#' the parameters \code{corr.model} and \code{corr.test} were called \code{sigma}
#' and \code{cr}. Also instead of supplying a graph object one could 
#' supply a parameter \code{weights} and a transition matrix \code{G}.
#' @return A list containg three elements
#' \describe{
#' \item{\code{LocalPower}}{A numeric giving the local powers for the hypotheses}
#' \item{\code{ExpRejections}}{The expected number of rejections}
#' \item{\code{PowAtlst1}}{The power to reject at least one hypothesis}
#' }
#' @references
#' 
#' Bretz, F., Maurer, W., Brannath, W. and Posch, M. (2009) A graphical
#' approach to sequentially rejective multiple test procedures. Statistics in
#' Medicine, 28, 586--604
#' 
#' Bretz, F., Maurer, W. and Hommel, G. (2010) Test and power considerations
#' for multiple endpoint analyses using sequentially rejective graphical
#' procedures, to appear in Statistics in Medicine
#' @keywords htest
#' @examples
#' 
#' ## reproduce example from Stat Med paper (Bretz et al. 2010, Table I)
#' ## first only consider line 2 of Table I
#' ## significance levels
#' graph <- simpleSuccessiveII()
#' ## alternative (mvn distribution)
#' corMat <- rbind(c(1, 0.5, 0.5, 0.5/2),
#'                 c(0.5,1,0.5/2,0.5),
#'                 c(0.5,0.5/2,1,0.5),
#'                 c(0.5/2,0.5,0.5,1))
#' theta <- c(3, 0, 0, 0)
#' calcPower(graph=graph, alpha=0.025, mean=theta, corr.sim=corMat, n.sim= 100000)
#' 
#' 
#' ## now reproduce all 14 simulation scenarios
#' ## different graphs
#' weights1 <- c(rep(1/2, 12), 1, 1)
#' weights2 <- c(rep(1/2, 12), 0, 0)
#' eps <- 0.01
#' gam1 <- c(rep(0.5, 10), 1-eps, 0, 0, 0)
#' gam2 <- gam1
#' ## different multivariate normal alternatives
#' rho <- c(rep(0.5, 8), 0, 0.99, rep(0.5,4))
#' th1 <- c(0, 3, 3, 3, 2, 1, rep(3, 7), 0)
#' th2 <- c(rep(0, 6), 3, 3, 3, 3, 0, 0, 0, 3)
#' th3 <- c(0, 0, 3, 3, 3, 3, 0, 2, 2, 2, 3, 3, 3, 3)
#' th4 <- c(0,0,0,3,3,3,0,2,2,2,0,0,0,0)
#' 
#' ## function that calculates power values for one scenario
#' simfunc <- function(nSim, a1, a2, g1, g2, rh, t1, t2, t3, t4, Gr){
#'   al <- c(a1, a2, 0, 0)
#'   G <- rbind(c(0, g1, 1-g1, 0), c(g2, 0, 0, 1-g2), c(0, 1, 0, 0), c(1, 0, 0, 0))
#'   corMat <- rbind(c(1, 0.5, rh, rh/2), c(0.5,1,rh/2,rh), c(rh,rh/2,1,0.5), c(rh/2,rh,0.5,1))
#'   mean <- c(t1, t2, t3, t4)
#'   calcPower(weights=al, alpha=0.025, G=G, mean=mean, corr.sim=corMat, n.sim = nSim)
#' }
#' 
#' ## calculate power for all 14 scenarios
#' outList <- list()
#' for(i in 1:14){
#'   outList[[i]] <- simfunc(10000, weights1[i], weights2[i], 
#'                     gam1[i], gam2[i], rho[i], th1[i], th2[i], th3[i], th4[i])
#' }
#' 
#' ## summarize data as in Stat Med paper Table I 
#' atlst1 <- as.numeric(lapply(outList, function(x) x$PowAtlst1))
#' locpow <- do.call("rbind", lapply(outList, function(x) x$LocalPower))
#' 
#' round(cbind(atlst1, locpow), 5)
#' 
#' @export calcPower
calcPower <- function(weights, alpha, G, mean = rep(0, nrow(corr.sim)),
                      corr.sim = diag(length(mean)), corr.test = NULL,
                      n.sim = 10000, type = c("quasirandom", "pseudorandom"),
					  f=list(), upscale=FALSE, graph, ...) {
  if (!is.null(list(...)[["nSim"]]) && missing(n.sim)) { n.sim <- list(...)[["nSim"]] }
  if (!is.null(list(...)[["sigma"]]) && missing(corr.sim)) { corr.sim <- list(...)[["sigma"]] }
  if (!is.null(list(...)[["cr"]]) && missing(corr.test)) { corr.test <- list(...)[["cr"]] }  
  if (!missing(graph)) {
      G <- graph@m
      weights <- graph@weights
  } 
	type <- match.arg(type)
	if (any(is.na(corr.sim))) stop("While parameter 'corr.test' can contain NAs, this does not make sense for 'corr.sim'.")
	#print(G)
	if (is.list(mean)) {
	  result <- list()
	  for (m in mean) {
		  sims <- rqmvnorm(n.sim, mean = m, sigma = corr.sim, type = type)
		  pvals <- pnorm(sims, lower.tail = FALSE)
		  out <- graphTest(pvalues=pvals, weights=weights, alpha=alpha, G=G, cr=corr.test, upscale=upscale)
		  out <- extractPower(out, f)
		  label <- attr(m, "label")		  
		  if (!is.null(label)) {
			  attr(out, "label") <- label 
		  }
		  result[[length(result)+1]] <- out 
	  }
	  return(result)
  } else {
	  sims <- rqmvnorm(n.sim, mean = mean, sigma = corr.sim, type = type)
	  pvals <- pnorm(sims, lower.tail = FALSE)
	  out <- graphTest(pvalues=pvals, weights=weights, alpha=alpha, G=G, cr=corr.test, upscale=upscale)
	  extractPower(out, f)
  }
}

calcMultiPower <- function(weights, alpha, G, ncpL, muL, sigmaL, nL,
		corr.sim = diag(length(ncpL[[1]])), corr.test = NULL,
		n.sim = 10000, type = c("quasirandom", "pseudorandom"),
		f=list(), digits=4, variables=NULL, test, upscale=FALSE, graph, ...) {
  if (!is.null(list(...)[["nSim"]]) && missing(n.sim)) { n.sim <- list(...)[["nSim"]] }
  if (!is.null(list(...)[["sigma"]]) && missing(corr.sim)) { corr.sim <- list(...)[["sigma"]] }
  if (!is.null(list(...)[["cr"]]) && missing(corr.test)) { corr.test <- list(...)[["cr"]] }  
  if (!missing(graph)) {
    G <- graph@m
    weights <- graph@weights
  } 
  ### END OF API COMPATIBILITY CODE  
  if (!is.null(list(...)[["sigma"]]) && missing(corr.sim)) {
    corr.sim <- list(...)[["sigma"]]
  }
  if (!is.null(list(...)[["cr"]]) && missing(corr.test)) {
    corr.test <- list(...)[["cr"]]
  }
  if (!missing(ncpL) && (!missing(muL)||!missing(sigmaL)||!missing(nL))) {
    warning("Only parameter 'ncpL' will be used, not 'muL', 'sigmaL' or 'nL'.")
  }
  if (missing(ncpL)) {
    ncpL <- list()
    for (mu in muL) {
      for (s in sigmaL) {
        for (n in nL) {
          newSetting <- mu*sqrt(n)/s
          attr(newSetting, "label") <- paste("mu: ",paste(mu,collapse=","),", sigma: ",paste(s,collapse=","),", n: ",paste(n,collapse=","),sep="")
          ncpL[[length(ncpL)+1]] <- newSetting 
        }
      }
    }
  } else {
    for (i in 1:length(ncpL)) {
      attr(ncpL[[i]], "label") <- names(ncpL)[i]
    }
  }
  if (length(f)>0) {
    n <- names(f)
    if (is.null(n) || all(is.na(n))) n <- paste("func", 1:length(f), sep="")
    n[n=="" | is.na(n)] <- paste("func", 1:sum(n==""), sep="")
    names(f) <- n
  }
  result <- data.frame(Scenario=character(0))
  vnames <- c()
  if (!is.null(variables)) {
    vnames <- names(variables)
    result <- cbind(result, as.data.frame(setNames(replicate(length(vnames), numeric(0), simplify = F), vnames)))
  }
  probs <- c(paste("LocalPower", names(weights)), "ExpRejections", "PowAtlst1", "RejectAll", names(f))
  result <- cbind(result, as.data.frame(setNames(replicate(length(probs), numeric(0), simplify = F), probs)))
  
	sResult <- ""
	g <- matrix2graph(G)
	g <- setWeights(g, weights)
	if (is.null(variables)) {
		sResult <- paste(sResult, "Graph:",paste(capture.output(print(g)), collapse="\n"), sep="\n")
		resultL <- calcPower(graph=g, alpha=alpha, mean = ncpL, corr.sim=corr.sim, corr.test=corr.test, n.sim=n.sim, type=type, f=f, upscale=upscale)
		result <- addResult2DF(result, resultL, digits=digits)
		sResult <- paste(sResult, resultL2Text(resultL, digits), sep="\n")
	} else {
	  graphs <- replaceVariables(graph, variables)
	  if (!is.list(graphs)) graphs <- list(graphs)
		for (GII in graphs) {
			additionalLabel <- "" #paste(",", attr(GII, "label"))
			variables <- attr(GII, "variables")
			resultL <- calcPower(graph=GII, alpha=alpha, mean = ncpL, corr.sim=corr.sim, corr.test=corr.test, n.sim=n.sim, type=type, f=f, upscale=upscale)
			result <- addResult2DF(result, resultL, additionalLabel=additionalLabel, digits=digits, variables=variables)
			sResult <- paste(sResult, resultL2Text(resultL, digits, additionalLabel=additionalLabel), sep="\n")
		}		
	}	
	#return(sResult)
	return(list(result, c("Scenario", vnames, probs)))
	#return(paste(capture.output(print(result)), collapse="\n"))
}

addResult2DF <- function(resultM, resultL, additionalLabel="", digits, variables=NULL) {
  resultRow <-  cbind( data.frame(Scenario=" "), as.data.frame(setNames(replicate(dim(resultM)[2]-1, 0, simplify = F), colnames(resultM)[-1])) )
  
  for(result in resultL) {
    resultRow[1] <- paste(attr(result, "label"), additionalLabel, sep="")
    skip <- 1
    if (!is.null(variables)) {
      for (i in 1:length(variables)) { 
        resultRow[1+i] <- variables[i]
      }
      skip <- i+1
    }
    for (i in 1:length(result$LocalPower)) { 
      resultRow[skip+i] <- result$LocalPower[i]
    }
    for (j in 2:length(result)) {
      resultRow[skip+i+j-1] <- result[[j]]
    }
    resultRow[2:dim(resultRow)[2]] <- round(resultRow[2:dim(resultRow)[2]], digits)
    
    resultM <- rbind(resultM, resultRow)
  }
  return(resultM)
}

resultL2Text <- function(resultL, digits, additionalLabel="") {
	sResult <- ""
	for(result in resultL) {
		label <- attr(result, "label")
		title <- paste("Setting: ", label, additionalLabel, sep="")		
		sResult <- paste(sResult, title, paste(rep("=", nchar(title)),collapse=""), sep="\n")			
		sResult <- paste(sResult, "Local Power:",paste(capture.output(print(round(result$LocalPower, digits))), collapse="\n"), sep="\n")
		sResult <- paste(sResult, "\nExpected number of rejections:", round(result$ExpRejections, digits), sep="\n")
		sResult <- paste(sResult, "Prob. to reject at least one hyp.:", round(result$PowAtlst1, digits), sep="\n")
		sResult <- paste(sResult, "Prob. to reject all hypotheses:", round(result$RejectAll, digits), sep="\n")
		if (length(result)>4) {
			for (i in 5:length(result)) {
				#TODO pF <- attr(result, "label")
				pF <- attr(result[i], "label")
				if (is.null(pF)) pF <- names(result)[i]
				sResult <- paste(sResult, paste(pF, ":", sep=""), result[i], sep="\n")
			}
		}
		sResult <- paste(sResult, "\n", sep="\n")		
	}
	return(sResult)
}

createCalcPowerCall <- function(alpha, ncpL, corr.sim = diag(length(ncpL[[1]])), corr.test = NULL,
                                n.sim = 10000, type = c("quasirandom", "pseudorandom"),
                                f="", digits=4, variables="", test, upscale=FALSE, graph, loop=TRUE, seed=1234) {	
  command <- dputGraph(graph, "graph")
  command <- paste(command, "\n", "ncpL <- ", ncpL,"\n", sep="")
  command <- paste(command, "\n", "f <- ", f,"\n", sep="")
  if (!missing(variables)) {
    command <- paste(command, "\n", "variables <- ", variables,"\n", sep="")
  }
  if (!missing(corr.test)) {
    command <- paste(command, "\n", dputMatrix(corr.test, name="corr.test", indent=TRUE),"\n", sep="")
  }
  if (!missing(corr.sim)) {
    command <- paste(command, "\n", dputMatrix(corr.sim, name="corr.sim", indent=TRUE),"\n", sep="")
  }
  command <- paste(command, "set.seed(", seed ,")\n\n", sep="")
  if (loop) {
    command <- paste(command, "df <- c()\n", sep="")
    if (!missing(variables)) {
      command <- paste(command, "for (g in replaceVariables(graph, variables, list=TRUE)) {\n", sep="")
    }
    command <- paste(command, "  result <- calcPower(graph=g, mean=ncpL, f=f", sep="")
  } else {
    command <- paste(command, "gMCP:::calcMultiPower(graph=graph, ncpL=ncpL, f=f", sep="")
    if (!missing(variables)) {
      command <- paste(command, ", variables=variables", sep="")
    }
  }
  if (!missing(test)) {
    command <- paste(command, ", test=\"",test,"\"", sep="")
  }
  if (!missing(type)) {
    command <- paste(command, ", type=\"",type,"\"", sep="")
  }
  if (upscale) {
    command <- paste(command, ", upscale=TRUE", sep="")
  }	
  if (!missing(corr.test)) {
    command <- paste(command, ", corr.test=corr.test", sep="")
  }
  if (!missing(corr.sim)) {
    command <- paste(command, ", corr.sim=corr.sim", sep="")
  }
  
  command <- paste(command, ", alpha=",dput2(alpha), sep="")
  command <- paste(command, ", n.sim=",dput2(n.sim), sep="")
  command <- paste(command, ")\n", sep="")
  if (loop) {
    command <- paste(command, "  df <- rbind(df, matrix(unlist(result), nrow=3, byrow=TRUE))\n", sep="")
    if (!missing(variables)) {
    command <- paste(command, "}\n", sep="")
    }
    command <- paste(command, "print(round(df,",digits,"))\n", sep="")
  }
  return(command)
}

#x <- calcMultiPower(weights=BonferroniHolm(3)@weights, alpha=0.05, G=BonferroniHolm(3)@m, muL=list(c(0,0,0),c(10,10,10),c(10,20,30)), sigmaL=list(c(1,1,1)), nL=list(c(10,10,10),c(20,20,20)), f=list(p1=function(x){x[1]&&x[2]}))
#cat(x)
