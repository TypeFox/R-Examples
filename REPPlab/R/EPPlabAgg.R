# Helpfunction (probably not required anymore, you can delete it in the next version)
  B2P_fixed <- function(x){tcrossprod(x)/sum(x^2)  }

  
#' Function to Aggregate Directions From epplab Objects
#'
#' Function that automatically aggregates the projection directions from one or more \code{epplab} objects.
#' Three options are available on how to choose the final projection which can have a rank larger than one.
#' The parameter \code{x} can either be a single object or a list of epplab objects.
#' Options for \code{method} are \code{inverse}, \code{sq.inverse} and \code{cumulative}.
#'
#' @param x An object of class \code{epplab} or a list of \code{epplab} objects.
#' @param method The type of method, see details. Options are \code{inverse}, \code{sq.inverse} and \code{cumulative}.
#' @param percentage Threshold for the relative eigenvalue sum to retain, see details.
#' @return A list with class 'epplabagg' containing the following components:
#' \item{P}{The estimated average orthogonal projection matrix.}
#' \item{O}{An orthogonal matrix on which P is based upon.}
#' \item{k}{The rank of the average orthogonal projection matrix.}
#' \item{eigenvalues}{The relevant eigenvalues, see details. Only given if \code{method="cumulative"}.}
#' @details Denote \eqn{p_i, i=1,...,m}, the projection vectors contained in the list of \code{epplab} objects 
#' and \eqn{P_i, i=1,..,m}, the corresponding orthogonal projection matrices (each having rank one). 
#' The method \code{cumulative} is based on the eigenvalue decomposition of 
#'   \eqn{P_w=\frac 1m \sum_{i=1}^m P_i}{P_w=1/m sum(P_i)}
#'  and transforms as \code{O} the eigenvectors such that the corresponding
#'  relative eigenvalues sum is at least \code{percentage}.
#'  The number of eigenvectors retained corresponds to the rank \code{k} and
#'  \code{P} is the corresponding orthogonal projection matrix.
#'  The methods \code{inverse} and \code{sq.inverse} are automatic rules to
#'  choose the number of eigenvectors to retain as implemented by the function
#'  \code{\link[LDRTools]{AOP}}.
#' @author Daniel Fischer, Klaus Nordhausen, Anne Ruiz-Gazen
#' @seealso \code{\link{EPPlab}}, \code{\link[LDRTools]{AOP}}
#' @references \cite{Liski, E., Nordhausen, K., Oja, H. and Ruiz-Gazen, A. (201?), Combining Linear Dimension Reduction Estimates, to appear in the Proceedings of \emph{ICORS 2015}, pp. ??-??.}
#' @keywords multivariate  
#' @examples
#' 
#'  library(tourr)
#'  data(olive)
#'  # To keep the runtime short, maxiter and n.simu were chosen very 
#'  # small for demonstration purposes, real life applications would
#'  # rather choose larger values, e.g. n.simu=100, maxiter=200
#'  olivePP.kurt.max <-
#'    EPPlab(olive[,3:10],PPalg="PSO",PPindex="KurtosisMax",n.simu=10, maxiter=20)
#'  
#'  olivePP.fried <-
#'    EPPlab(olive[,3:10],PPalg="PSO",PPindex="Friedman",n.simu=10, maxiter=20)
#'  
#'  olivePPs <- list(olivePP.kurt.max, olivePP.fried)
#'  
#'  EPPlabAgg(olivePP.kurt.max)$k
#'  EPPlabAgg(olivePPs, "cum", 0.99)$k
#'  
#'  pairs(olivePP.kurt.max$x %*% EPPlabAgg(olivePPs, "cum", 0.99)$O,
#'        col=olive[,2], pch=olive[,1])
#'  
#'  
#'  olivAOP.sq <- EPPlabAgg(olivePPs, "inv")
#'  oliveProj <- olivePP.kurt.max$x %*% olivAOP.sq$O
#'  plot(density(oliveProj))
#'  rug(oliveProj[olive$region==1],col=1)
#'  rug(oliveProj[olive$region==2],col=2)
#'  rug(oliveProj[olive$region==3],col=3)
#' 
#' @export EPPlabAgg
EPPlabAgg <- function(x, method="cumulative", percentage=0.85){
  
      # Input checks
        method <- match.arg(method, c("inverse", "sq.inverse", "cumulative"))
        if(class(x)=="epplab") x <- list(x)  
      
      # Now combine the results
        switch(method,
               cumulative={
                # Store the the averages, initialize with 0
                  avgMatrix <- matrix(0, nrow=dim(x[[1]]$x)[2], ncol=dim(x[[1]]$x)[2])
                # Go through each REPPlab result object (using this option, length(x) is probably always equal to 1)
                  for(i in 1:length(x)){
                  # Now go through all directions 
                    for(dirRun in 1:dim(x[[i]]$PPdir)[2]){
                  # Now sum them up
                      avgMatrix <- avgMatrix + tcrossprod(x[[i]]$PPdir[, dirRun])
                    }
                  }
                # Divide by the summands to get the average
                  avgMatrix <- avgMatrix / (dim(x[[i]]$PPdir)[2] * length(x))
                # Calculate the eigenelements
                  eigmave<-eigen(avgMatrix) 
                  lmave<-eigmave$values
                  umave<-eigmave$vectors
                # eliminate the directions that are associated with less than '1-percentage'% of the information
                  takeThese <- 1:min(sum(((cumsum(lmave)/sum(lmave))<percentage))+1, length(lmave))
                  keepmave<-umave[,takeThese,drop=FALSE] 
                # project the data on the directions we keep
                  coord <- x[[i]]$x %*% keepmave           
                # Write out the results
                  res <- list(P=O2P(umave), O=keepmave, k= ncol(keepmave), eigen=lmave)      
                },
                {
                  B2P.output <-  list()
                  lresB2P <- list()
                  lresB2P.all <- c()
                  for(i in 1:length(x)){
                    B2P.output[[i]] <- apply(coef(x[[i]]),2 , B2P)
                    lresB2P[[i]] <- tapply(B2P.output[[i]], gl(ncol(B2P.output[[i]]), nrow(B2P.output[[i]])), matrix, nrow=dim(x[[1]]$x)[2], ncol=dim(x[[1]]$x)[2])
                    lresB2P.all <- c(lresB2P.all, lresB2P[[i]])
                  }
                  res <- AOP(lresB2P.all, weights = method)  
                })
      # Return the results
        class(res) <- "epplabagg"
        res
}
