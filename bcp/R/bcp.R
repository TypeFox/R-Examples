#
# bcp: an R package for performing a Bayesian analysis
# of change point problems.
#
# Copyright (C) 2011 Chandra Erdman and John W. Emerson
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, a copy is available at
# http://www.r-project.org/Licenses/
#
#-------------------
# FILE: bcp.R

#' @title Performs a Bayesian analysis of change point problems
#' 
#' @description
#' \code{bcp()} implements the Bayesian change point analysis methods given in Wang and Emerson (2015), of which the Barry and Hartigan (1993) product
#' partition model for the normal errors change point problem is a specific case. 1. Multivariate (or univariate) Bayesian change point analysis: We assume there exists an unknown partition of a data series y
#' into blocks such that the mean is constant within each block.  In the multivariate
#' case, a common change point structure
#' is assumed; means are constant within each block of each sequence, but may differ across sequences
#' within a given block. Conditional on the partition, the model assumes that observations are independent, identically distributed normal, with constant means within blocks and
#' constant variance throughout each sequence.  
#' 2. Linear regression Bayesian change point analysis: As with the previous model, we assume the observations (x,y), where x may be multivariate, are partitioned into blocks, and that linear models are appropriate within each block. 
#' 
#' If an adjacency structure is provided, the data are assumed to reside on nodes of a graph with the given adjacency structure; additional parameters are used in this graph change point model. If no adjacency structure is provided, the data are assumed to be sequential and the blocks are forced to be contiguous. 
#'
#' @details 
#' The primary result is an estimate of the posterior mean (or its distribution if
#' \code{return.mcmc} is \code{TRUE}) at every location.  Unlike a frequentist or
#' algorithmic approach to the problem, these estimates will not be constant within
#' regions, and no single partition is identified as best.  Estimates of the
#' probability of a change point at any given location are provided, however.
#' 
#' The user may set \code{.Random.seed} to control the MCMC iterations.
#' 
#' The functions \code{\link{summary.bcp}}, \code{\link{print.bcp}}, and \code{\link{plot.bcp}} are
#' used to obtain summaries of the results; \code{\link{legacyplot}} is included
#' from package versions prior to 3.0.0 and will only work for univariate change
#' point analyses.
#' 
#' @param y a vector or matrix of numerical data (with no missing values). For
#' the multivariate change point problems, each column corresponds to a series.
#' @param x (optional) a matrix of numerical data (with no missing values) representing the predicting variables for a linear regression. 
#' @param id (optional) a vector of integers specifying the location ID for each observation in \code{y}, starting from location 1.
#' @param adj (optional) an adjacency list. Indexing the observations from 1 to \eqn{n}, the \eqn{i}-th element of the list is a vector of indices (offset by 1) for nodes that share an edge with node \eqn{i}.
#' @param w0 (optional) a single numeric value in the multivariate case or a vector of values in the regression case; in both, the value(s), between 0 and 1, is/are the parameter(s) in the uniform prior(s) on the signal-to-noise
#'          ratio(s).  If no value is specified, the default value of 0.2 is used, as    
#'          recommended by Barry and Hartigan (1993).
#' @param p0 (optional) a value between 0 and 1. For sequential data, it is the parameter of the prior on change point probabilities, \eqn{U(0,}\code{ p0}\eqn{)}, on the probability
#'          of a change point at each location in the sequence; for data on a graph, it is the parameter in the partition prior, \code{p0}\eqn{^{l(\rho)}}, where \eqn{l(\rho)} is the boundary length of the partition. 
#' @param d (optional) a positive number only used for linear regression change point models. Lower \code{d} means higher chance of fitting the full linear model (instead of the intercept-only model); see prior for \eqn{\tau_S} in Wang and Emerson (2015).
#' @param burnin the number of burnin iterations.
#' @param mcmc the number of iterations after burnin. 
#' @param return.mcmc logical. If set to \code{TRUE}, the posterior means and the partitions
#'                   in each iteration are returned.
#' @param boundaryType (optional) only applicable for graph change point analysis. Values can be ``node'' (default) if we count nodes in the boundary length calculation, or ``edge'' if we count edges in the boundary length calculation. See Wang and Emerson (2015) for details.
#' @param p1 (optional) only applicable for graph change point analysis. The proportion of Active Pixel Passes run that are the actual Active Pixel Passes specified in Barry and Hartigan (1994). \code{p1 = 0} corresponds to exclusively using the pseudo-Active Pixel Passes given in Wang and Emerson (2015). 
#' @param freqAPP (optional) only applicable for graph change point analysis. A positive integer for the number of Active Pixel Passes run in each step of the MCMC algorithm. 
#'
#' @return Returns a list containing the following components:
#' \describe{
#'   \item{data}{a copy of the data.}
#' \item{return.mcmc}{\code{TRUE} or \code{FALSE} as specified by the user; see the arguments, above.}
#' \item{mcmc.means}{if \code{return.mcmc=TRUE}, \code{mcmc.means} contains the means for each iteration conditional on the state of the partition.}
#' \item{mcmc.rhos}{if \code{return.mcmc=TRUE}, \code{mcmc.rhos} contains the partitions after each iteration. A value of 1 indicates the end of a block.}
#' \item{blocks}{a vector of the number of blocks after each iteration.}
#' \item{posterior.mean}{a vector or matrix of the estimated posterior means. In the regression case, the matrix includes posterior means for the response variable.}
#' \item{posterior.var}{a vector or matrix of the estimated posterior variances. In the regression case, the estimated posterior variances of the response are provided.}
#' \item{posterior.prob}{a vector of the estimated posterior probabilities of changes at each location.}
#' \item{burnin}{the number of burnin iterations.}
#' \item{mcmc}{the number of iterations after burnin.}
#' \item{w0}{see the arguments, above.}
#' \item{p0}{see the arguments, above.}
#' }
#' @references 
#' \enumerate{
#' \item J. Bai and P. Perron (2003), Computation and Analysis of Multiple Structural Change Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22. \url{http://qed.econ.queensu.ca/jae/2003-v18.1/bai-perron/}.
#' 
#' \item Daniel Barry and J. A. Hartigan (1993), A Bayesian Analysis for Change Point Problems, \emph{Journal of The American Statistical Association}, \bold{88}, 309-19.
#' 
#' \item Daniel Barry and J. A. Hartigan (1994), A Product Partition Model for Image Restoration, \emph{New Directions in Statistical Data Analysis and Robustness}, (Monte Verita : Proceedings of the Cento Stefano Franscini Ascona), Birkhauser, 9-23.
#' 
#' \item Chandra Erdman and John W. Emerson (2008), A Fast Bayesian Change Point Analysis for the Segmentation of Microarray Data, \emph{Bioinformatics}, 24(19), 2143-2148. \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btn404}.
#' 
#' \item Chandra Erdman and John W. Emerson (2007), bcp: An R Package for Performing a Bayesian Analysis of Change Point Problems. \emph{Journal of Statistical Software}, 23(3), 1-13. \url{http://www.jstatsoft.org/v23/i03/}.
#' 
#' \item A. B. Olshen, E. S. Venkatraman, R. Lucito, M. Wigler (2004), Circular binary segmentation for the analysis of array-based DNA copy number data, \emph{Biostatistics}, \bold{5}, 557-572.  \url{http://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html}.
#' 
#' \item Snijders \emph{et al.} (2001), Assembly of microarrays for genome-wide measurement of DNA copy number, \emph{Nature Genetics}, \bold{29}, 263-264. 
#' 
#' \item Xiaofei Wang and John W. Emerson (2015). Bayesian Change Point Analysis of Linear Models on General Graphs, \emph{Working Paper}.
#' 
#' \item Achim Zeileis, Friedrich Leisch, Kurt Hornik, Christian Kleiber (2002), strucchange: An R Package for Testing for Structural Change in Linear Regression Models, \emph{Journal of Statistical Software}, \bold{7}(2), 1-38. \url{http://www.jstatsoft.org/v07/i02/}. 
#' }
#' @author Xiaofei Wang, Chandra Erdman, and John W. Emerson
#' @seealso \code{\link{plot.bcp}}, \code{\link{summary.bcp}}, and \code{\link{print.bcp}} for summaries of the results.
#' @import graphics
#' @examples
#' 
#' ##### univariate sequential data #####
#' # an easy problem with 2 true change points
#' set.seed(5)
#' x <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
#' bcp.1a <- bcp(x)
#' plot(bcp.1a, main="Univariate Change Point Example")
#' legacyplot(bcp.1a)
#' 
#' # a hard problem with 1 true change point
#' set.seed(5)
#' x <- rep(c(0,1), each=50)
#' y <- x + rnorm(50, sd=1)
#' bcp.1b <- bcp(y)
#' plot(bcp.1b, main="Univariate Change Point Example")
#' 
#' ##### multivariate sequential data #####
#' # an easy problem in k=3 dimensions
#' set.seed(5)
#' x <- rnorm(6, sd=3)
#' y <- rbind(cbind(rnorm(50, x[1]), rnorm(50, x[2]), rnorm(50, x[3])),
#'            cbind(rnorm(50, x[4]), rnorm(50, x[5]), rnorm(50, x[6])))
#' bcp.2a <- bcp(y)
#' plot(bcp.2a, main="Multivariate (k=3) Change Point Example")
#' plot(bcp.2a, separated=TRUE, main="Multivariate (k=3) Change Point Example")
#' 
#' # a harder problem in k=5 dimensions
#' set.seed(5)
#' means1 <- rep(0, 5)
#' means2 <- rep(1, 5)
#' x <- rbind(matrix(rep(means1, each=50), nrow=50),
#'            matrix(rep(means2, each=50), nrow=50))
#' y <- x + rnorm(length(x), sd=1)
#' bcp.2b <- bcp(cbind(y))
#' plot(bcp.2b, main="Multivariate (k=5) Change Point Example")
#' 
#' ##### linear models with sequential data #####
#' # 1 true change point at location 50; the predicting variable x is not related to location
#' x <- rnorm(100)
#' b <- rep(c(3,-3), each=50)
#' y <- b*x + rnorm(100)
#' bcp.3a <- bcp(y, x)
#' # in the two plots that follow, the location IDs are used as the plot characters
#' par(mfrow=c(1,2))
#' plot(y ~ x, type="n", main="Linear Regression: Raw Data")
#' text(x, y, as.character(1:100), col=(b/3)+2)
#' plot(y ~ x, type="n", main="Linear Regression: Posterior Means")
#' text(x, bcp.3a$posterior.mean[,1], as.character(1:100), col=(b/3)+2)
#' plot(bcp.3a, main="Linear Regression Change Point Example")
#' 
#' # 1 true change point at location 50; the predicting variable x is equal to location
#' x <- 1:100
#' b <- rep(c(3,-3), each=50)
#' y <- b*x + rnorm(100, sd=50)
#' bcp.3b <- bcp(y, x)
#' plot(bcp.3b, main="Linear Regression Change Point Example")
#' 
#' ##### univariate data on a grid #####  
#' set.seed(5)
#' adj <- makeAdjGrid(20)
#' z <- rep(c(0, 2), each=200)
#' y <- z + rnorm(400, sd=1)
#' out <- bcp(y, adj=adj, burnin=500, mcmc=500)
#' 
#' if (require("ggplot2")) {
#'   df <- data.frame(mean=z, data = y, post.means = out$posterior.mean[,1], 
#'                    post.probs = out$posterior.prob, 
#'                    i = rep(1:20, each=20), j = rep(1:20, times=20))
#' 
#'   # visualize the means
#'   g <- ggplot(df, aes(i,j)) + 
#'          geom_tile(aes(fill = mean), color='white') +
#'          scale_fill_gradientn(limits=range(y), colours=c('white', 'steelblue'))+
#'          ggtitle("True Means")
#'   print(g)
#' 
#'   # visualize the data
#'   g <- ggplot(df, aes(i,j)) + 
#'          geom_tile(aes(fill = data), color='white') +
#'          scale_fill_gradientn(limits=range(y), colours=c('white', 'steelblue'))+
#'          ggtitle("Observed Data")
#'   print(g)
#' 
#'   # visualize the posterior means/probs
#'   g <- ggplot(df, aes(i,j)) + 
#'          geom_tile(aes(fill = post.means), color='white') +
#'          scale_fill_gradientn(limits=range(y), colours=c('white', 'steelblue'))+
#'          ggtitle("Posterior Means")
#'   print(g)
#' 
#'   g <- ggplot(df, aes(i,j)) + 
#'          geom_tile(aes(fill = post.probs), color='white') +
#'          scale_fill_gradientn(limits=c(0, 1), colours=c('white', 'steelblue'))+
#'          ggtitle("Posterior Boundary Probabilities")
#'   print(g)
#' }
#' 
#' 
#' \dontrun{
#' ##### multivariate data on a grid #####
#' set.seed(5)
#' x <- rnorm(6, sd=3)
#' y <- rbind(cbind(rnorm(50, x[1]), rnorm(50, x[2]), rnorm(50, x[3])),
#'            cbind(rnorm(50, x[4]), rnorm(50, x[5]), rnorm(50, x[6])))
#' adj <- makeAdjGrid(10)
#' a <- bcp(y, adj=adj, p0=0.4, burnin=500, mcmc=500)
#' 
#' ##### linear models on a grid #####
#' set.seed(5)
#' x <- rnorm(100)
#' b <- rep(c(3,-3), each=50)
#' y <- b*x + rnorm(100)
#' adj <- makeAdjGrid(10)
#' a <- bcp(y,x,adj=adj, p0=0.4, burnin=500, mcmc=500)
#' 
#' ##### linear models on a grid using pseudo-APPs #####
#' x <- rnorm(100)
#' b <- rep(c(3,-3), each=50)
#' y <- b*x + rnorm(100)
#' adj <- makeAdjGrid(10)
#' a <- bcp(y,x,adj=adj, p0=0.4, burnin=500, mcmc=500, p1 = 0)}
#' 
#' ##### univariate data on a graph ##### 
#' \dontrun{ 
#'   demo(bcpgraph)
#' }
#' 
#' ###### Real Data Examples ######
#' \dontrun{
#' # Coriell chromosome 11: univariate sequential data
#' demo(coriell)
#' 
#' # U.S. ex-post interest rate: univariate sequential data
#' demo(RealInt)
#' 
#' # Lombard: univariate sequential data (with and without linear models)
#' demo(Lombard)
#' 
#' # Quebec rivers: multivariate sequential data
#' demo(QuebecRivers)
#' 
#' # New Haven housing: linear models on a graph
#' demo(NewHaven)
#' } 
#' 
#' @keywords datasets
#' @export
"bcp" <- function(y, x=NULL, id = NULL, adj=NULL, w0=NULL, p0=0.2, d = 10, 
                  burnin=50, mcmc=500, return.mcmc=FALSE, 
                  boundaryType = "node", p1 = 1, freqAPP = 20) {
  
  ######################################################
  ########################### BEGIN THE WORKER FUNCTION:
  ######################################################
  
  "worker.bcp" <- function(mcmc, y, x, id, w0, p0, d, burnin, return.mcmc, membinit,
                           boundaryType, adj, p1, freqAPP) {
    
    require(bcp)
    
    # INITIALIZATION
    #  if (is.data.frame(x)) x <- matrix(as.double(x), nrow=nrow(x), ncol=ncol(x))
    #  if (is.vector(x)) x <- matrix(as.double(x), ncol=1)
    #  if (!is.matrix(x)) stop("x must be a vector, matrix or a data frame")
    #  if (nrow(x)==1) {
    #    warning("coercing data to a single series")
    #    x <- matrix(as.vector(x), ncol=1)
    #  }
    if (p0 > 1 | p0 < 0 | p1 > 1 | p1 < 0) stop ("p0 and p1 must each be between 0 and 1.")
    if (is.null(id)) {
      if (is.matrix(y)) id <- 1:nrow(y)
      else id <- 1:length(y)
    }
    if (min(id) == 1) id <- id - 1
    if (is.null(x)) {
      # doing multivariate bcp
      if (is.null(w0)) w0 <- 0.2
      if (w0 > 1 | w0 < 0) stop("w0 must be between 0 and 1.")
      if (is.vector(y)) y <- cbind(y)
      # Do the work in C:
      if (is.null(adj)) {
        out <- rcpp_bcpM(
          y, as.integer(id),
          as.integer(return.mcmc),
          as.integer(burnin), 
          as.integer(mcmc),     
          as.double(p0),
          as.double(w0))      
        attr(out, "structure") <- "series"
      } else {
        if (freqAPP < 0) stop("freqAPP must be a positive integer.")
        out <- rcpp_ppm(y, as.integer(id),
                        adj,
                        as.integer(return.mcmc),
                        as.integer(burnin), 
                        as.integer(mcmc),     
                        as.double(p0),
                        as.double(w0),
                        membinit,
                        as.integer(boundaryType),
                        as.double(p1),
                        as.integer(freqAPP))
        attr(out, "structure") <- "graph"
      }
      attr(out, "model") <- "multivariate"
      out$data <- cbind(id+1, y) 
      if (is.null(colnames(x))) {
        colnames(out$posterior.mean) <- paste(rep("X", ncol(out$data)-1), 
                                              1:(ncol(out$data)-1), sep="")
      } else
        colnames(out$posterior.mean) <- colnames(x)
    } else {
      # doing regression bcp
      if (is.vector(x)) x <- cbind(x)
      if (sum(x[,1]==1) != nrow(x)) x <- cbind(1, x) 
      if (!is.double(x)) x <- matrix(as.double(x), nrow(x), ncol(x))
      if (is.null(w0)) {
        w0 <- rep(0.2, ncol(x))
      } else {
        if (any(w0 > 1) | any(w0 < 0)) stop("Each element in w0 must be between 0 and 1.")
        if (length(w0) != ncol(x)) {
          if (length(w0) == 1) {
            w0 <- rep(w0, ncol(x))
            print("Incorrect length for w0. I'll assume you wanted each error-to-signal ratio to be iid from U(0, w0).")
          } else stop("Incorrect length for w0.")
        }
      }  
        
        indmat <- t(sapply(unique(id), function(y) id == y)*1)   
        
        # Do the work in C:
        if (is.null(adj)) {
          out <- rcpp_bcpR( 
            as.double(y),
            x,
            indmat,
            as.integer(id),
            as.integer(return.mcmc),
            as.integer(burnin), 
            as.integer(mcmc),     
            as.double(p0),
            as.double(w0),
            d
          )
          out$data = cbind(id+1, y,x[,-1])
          attr(out, "structure") <- "series"
        } else {
          if (freqAPP < 0) stop("freqAPP must be a positive integer.")
          out <- rcpp_ppmR(y, x,                    
                           indmat,
                           as.integer(id), 
                           adj,
                           as.integer(return.mcmc),
                           as.integer(burnin), 
                           as.integer(mcmc),     
                           as.double(p0),
                           as.double(w0),
                           membinit,
                           as.integer(boundaryType),
                           d,
                           as.double(p1),
                           as.integer(freqAPP))
          out$data = cbind(id+1, y,x)
          attr(out, "structure") <- "graph"
        }
        attr(out, "model") <- "regression"
        out$posterior.mean <- cbind(as.numeric(out$posterior.mean)/table(id))
        rownames(out$posterior.mean) <- NULL
        colnames(out$posterior.mean) <- "y"
      }
      if (attr(out, "structure") == "series") 
        out$posterior.prob[length(out$posterior.prob)] <- NA     # Fix up the last position, always NA
      
      # RETURN RESULTS
      z <- list(data=out$data,
                return.mcmc=return.mcmc,
                mcmc.means=out$mcmc.means,
                mcmc.rhos=out$mcmc.rhos,
                blocks=out$blocks,
                posterior.mean=out$posterior.mean,
                posterior.var=out$posterior.var,
                posterior.prob=out$posterior.prob,
                burnin=burnin,            
                mcmc=mcmc,     
                p0=p0,     
                w0=w0)     
      attr(z, "model") <- attr(out, "model")
      attr(z, "structure") <- attr(out, "structure")
      class(z) <- "bcp"
      return(z)
    }
    
    ###################################################
    ########################### END THE WORKER FUNCTION
    ########################### BEGIN THE MAIN SECTION:
    ###################################################
    
    # Function header and foreach setup, from above:
    #
    #"bcp" <- function(x, w0=0.2, p0=0.2, burnin=50, mcmc=500, return.mcmc=FALSE) {
    #
    
    if (!is.null(adj)) {
      if (boundaryType == "node") {
        boundaryType <- 1
      } else {
        boundaryType <- 2
      }
      
      if (is.vector(y)) {
        dataToRank <- y
      } else if (is.matrix(y)) {
        dataToRank <- y[,1]
      }
      if (!is.null(id)) { # varying num of obs per loc, we'll sample one obs per loc to rank
        inds <- sapply(1:max(id), function(g) {
          inds <- which(id==g)
          if (length(inds)==1) return(inds)
          return(sample(inds, 1))
        })
        dataToRank <- dataToRank[inds]
      }
      numNodes <- length(dataToRank)
      # if (is.null(membinit)) {
      Minit <- ceiling(sqrt(numNodes))
      o <- rank(dataToRank, ties.method="first")    
      membinit <- ceiling(o/Minit)
      membinit <- pmin(membinit, Minit)-1
      # } else if (length(membinit) == 1) {
      #   Minit <- ceiling(numNodes/membinit)
      #   o <- rank(dataToRank, ties.method="first")    
      #   membinit <- ceiling(o/Minit)
      #   membinit <- pmin(membinit, Minit)-1
      # } else {
      #   nComponents <- max(membinit)
      #   if (length(setdiff(unique(membinit), 1:nComponents))>0) {
      #     stop("Error in membinit")
      #   } else {
      #     membinit <- membinit-1
      #   }
      # }
      relabelMap <- order(unique(membinit))
      membinit <- relabelMap[membinit+1]-1
    }
    
    ans <- worker.bcp(mcmc, y=y, x=x, id=id, w0=w0, p0=p0, d=d, 
                      burnin=burnin, return.mcmc=return.mcmc,
                      membinit=membinit, boundaryType=boundaryType,
                      adj=adj, p1=p1, freqAPP=freqAPP)
    #  ==================================
    #  ===  Reformat the mcmc.means   ===
    #  ==================================
    if (return.mcmc) {
      if (ncol(ans$mcmc.means) > 1) {
        mcmc.means <- vector('list', ncol(ans$mcmc.means)) 
        for (i in 1:length(mcmc.means)) {
          mcmc.means[[i]] <- matrix(ans$mcmc.means[,i], nrow=burnin+mcmc, byrow=TRUE)
        }
      } else {
        mcmc.means <- matrix(ans$mcmc.means, nrow=burnin+mcmc, byrow=TRUE)
      }
      ans$mcmc.means <- mcmc.means
    } 
    return(ans)
  }
  
  #' @title Creating the adjacency structure for grid graphs
  #' 
  #' @description
  #' \code{makeAdjGrid()} produces a sparse representation of the adjacency structure for grid graphs, useful as the \code{adj} argument in \code{bcp()}. 
  #'
  #' @param n the number of rows of vertices in the graph data.
  #' @param m (optional) the number of column of vertices in the graph data. If not given, we assume \code{m = n}.
  #' @param k (optional) the number of neighbors assumed for a typical vertex (see details below), either 4 or 8. Default number of neighbors is assumed to be 8.
  #' 
  #' @author Xiaofei Wang
  #' @details 
  #' \code{makeAdjGrid()} produces a list representation of the adjacency structure for grid graphs. The \eqn{i}-th entry in the list gives a vector of neighbor ids for the \eqn{i}-th node. Note that neighbor ids are offset by 1 because indexing starts at 0 in C++.
  #' If \code{k = 8}, then we assume each node is joined via edges to its 8 neighbors in the (top left, top middle, top right, left, right, bottom left, bottom middle, and bottom right) directions, where applicable. If \code{k = 4}, then we assume each node is joined via edges to its 4 neighbors in the (top, right, bottom, left) directions, where applicable.
  #' @seealso \code{\link{bcp}} for performing Bayesian change point analysis.
  #' @examples 
  #' # generates an adjacency list for a 10 node by 5 node grid, assuming a maximum of 8 neighbors
  #' adj <- makeAdjGrid(10, 5) 
  #' 
  #' # generates an adjacency list for a 10 node by 5 node grid, assuming a maximum of 4 neighbors
  #' adj4 <- makeAdjGrid(10, 5, 4)
  #' 
  #' 
  #' ### show a grid example
  #' \dontrun{
  #' set.seed(5)
  #' adj <- makeAdjGrid(20)
  #' z <- rep(c(0, 2), each=200)
  #' y <- z + rnorm(400, sd=1)
  #' out <- bcp(y, adj=adj, burnin=500, mcmc=500)
  #' 
  #' if (require("ggplot2")) {
  #'   df <- data.frame(mean=z, data = y, post.means = out$posterior.mean[,1], 
  #'                    post.probs = out$posterior.prob, 
  #'                    i = rep(1:20, each=20), j = rep(1:20, times=20))
  #' 
  #'   # visualize the data
  #'   g <- ggplot(df, aes(i,j)) + 
  #'          geom_tile(aes(fill = data), color='white') +
  #'          scale_fill_gradientn(limits=range(y), colours=c('white', 'steelblue'))+
  #'          ggtitle("Observed Data")
  #'   print(g)
  #' 
  #'   # visualize the means
  #'   g <- ggplot(df, aes(i,j)) + 
  #'          geom_tile(aes(fill = mean), color='white') +
  #'          scale_fill_gradientn(limits=range(y), colours=c('white', 'steelblue'))+
  #'          ggtitle("True Means")
  #'   print(g)
  #' 
  #'   # visualize the posterior means/probs
  #'   g <- ggplot(df, aes(i,j)) + 
  #'          geom_tile(aes(fill = post.means), color='white') +
  #'          scale_fill_gradientn(limits=range(y), colours=c('white', 'steelblue'))+
  #'          ggtitle("Posterior Means")
  #'   print(g)
  #' 
  #'   g <- ggplot(df, aes(i,j)) + 
  #'          geom_tile(aes(fill = post.probs), color='white') +
  #'          scale_fill_gradientn(limits=c(0, 1), colours=c('white', 'steelblue'))+
  #'          ggtitle("Posterior Boundary Probabilities")
  #'   print(g)
  #' }
  #' }
  #' 
  #' @keywords datasets
  #' @export
  makeAdjGrid <- function(n,m=NULL, k=8) {
    if (is.null(m)) m <- n
    adj <- vector('list', n*m)
    
    if (k == 8) {
      for (i in 2:(n-1)) {
        for (j in 2:(m-1)) {
          adj[[(j-1)*n+i]] <- c((j-2)*n+(i-1):(i+1),(j-1)*n+i-c(1,-1), j*n+(i-1):(i+1))-1
        }
      }
      i <- 1
      for (j in 2:(m-1)) adj[[(j-1)*n+i]] <- c((j-2)*n+1:2, (j-1)*n+2, (j)*n+1:2)-1
      i <- n
      for (j in 2:(m-1))  adj[[(j-1)*n+i]] <- c((j-1)*n-1:0, j*n-1, (j+1)*n-1:0)-1
      j <- 1
      for (i in 2:(n-1)) adj[[(j-1)*n+i]] <- c(i-1,i+1, n+(i-1):(i+1))-1
      j <- m
      for (i in 2:(n-1)) adj[[(j-1)*n+i]] <- c((m-2)*n+(i-1):(i+1), (m-1)*n+i-1, (m-1)*n+i+1)-1
      adj[[1]] <- c(2, n+1:2)-1
      adj[[n]] <- c((1:2)*n-1,2*n)-1
      adj[[(m-1)*n+1]] <-  c((m-2)*n+1:2, (m-1)*n+2)-1
      adj[[n*m]] <- c((m-1)*n-(1:0), n*m-1)-1
    } else if (k == 4) {
      for (i in 2:(n-1)) {
        for (j in 2:(m-1)) {
          adj[[(j-1)*n+i]] <- c((j-2)*n+i,(j-1)*n+i-c(1,-1), j*n+i)-1
        }
      }
      i <- 1
      for (j in 2:(m-1)) adj[[(j-1)*n+i]] <- c((j-2)*n+1, (j-1)*n+2, (j)*n+1)-1
      i <- n
      for (j in 2:(m-1))  adj[[(j-1)*n+i]] <- c((j-1)*n, j*n-1, (j+1)*n)-1
      j <- 1
      for (i in 2:(n-1)) adj[[(j-1)*n+i]] <- c(i-1,i+1, n+i)-1
      j <- m
      for (i in 2:(n-1)) adj[[(j-1)*n+i]] <- c((n-2)*m+i, (n-1)*m+i-1, (n-1)*m+i+1)-1
      adj[[1]] <- c(2, n+1)-1
      adj[[n]] <- c(n-1,2*n)-1
      adj[[(m-1)*n+1]] <-  c((m-2)*n+1, (m-1)*n+2)-1
      adj[[n*m]] <- c((m-1)*n, n*m-1)-1
    } else {
      stop("Error: k must be 4 or 8.")
    }
    return(adj)
  }
  
  
  
