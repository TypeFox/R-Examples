#' Make CRF features
#' 
#' Make the data structure of CRF features
#' 
#' This function makes the data structure of features need for modeling and training CRF.
#'
#' The parameters \code{n.nf} and \code{n.ef} specify the number of node and edge features,
#' respectively.
#' 
#' The objects \code{node.par} and \code{edge.par} define the corresponding
#' parameters used with each feature. \code{node.par} is a 3-dimensional arrays, 
#' and element \code{node.par[n,i,f]} is the index of parameter associated with the 
#' corresponding node potential \code{node.pot[n,i]} and node feature \code{f}.
#' \code{edge.par} is a list of 3-dimensional arrays, and element 
#' \code{edge.par[[e]][i,j,f]} is the index of parameter associated with the 
#' corresponding edge potential \code{edge.pot[[e]][i,j]} and edge feature \code{f}.
#' The value 0 is used to indicate the corresponding node or edge potential 
#' does not depend on that feature.
#' 
#' For detail of calculation of node and edge potentials from features and parameters,
#' please see \code{\link{crf.update}}.
#' 
#' @param crf The CRF
#' @param n.nf The number of node features
#' @param n.ef The number of edge features
#' @return This function will directly modify the CRF and return the same CRF.
#' 
#' @seealso \code{\link{crf.update}}, \code{\link{make.par}}, \code{\link{make.crf}}
#' 
#' @export
make.features <- function(crf, n.nf = 1, n.ef = 1)
{
	crf$n.nf <- n.nf
	crf$n.ef <- n.ef
	crf$node.par <- array(0, dim=c(crf$n.nodes, crf$max.state, n.nf))
	crf$edge.par <- lapply(1:crf$n.edges, function(i) array(0, dim=c(crf$n.states[crf$edges[i,1]], crf$n.states[crf$edges[i,2]], n.ef)))
	crf
}



#' Make CRF parameters
#' 
#' Make the data structure of CRF parameters
#' 
#' This function makes the data structure of parameters need for modeling and training CRF.
#' The parameters are stored in \code{par}, which is a numeric vector of length \code{n.par}.
#' 
#' @param crf The CRF
#' @param n.par The number of parameters
#' @return This function will directly modify the CRF and return the same CRF.
#'
#' @seealso \code{\link{crf.update}}, \code{\link{make.features}}, \code{\link{make.crf}}
#' 
#' @export
make.par <- function(crf, n.par = 1)
{
	crf$n.par <- n.par
	crf$par <- numeric(crf$n.par)
	crf$nll <- numeric(1)
	crf$gradient <- numeric(crf$n.par)
	crf
}



#' Update MRF potentials
#' 
#' Update node and edge potentials of MRF model 
#' 
#' The function updates \code{node.pot} and \code{edge.pot} of MRF model.
#' 
#' @param crf The CRF
#' @return This function will directly modify the CRF and return the same CRF.
#' 
#' @seealso \code{\link{mrf.nll}}, \code{\link{train.mrf}}
#' 
#' @export
mrf.update <- function(crf)
  .Call(MRF_Update, crf)



#' Update CRF potentials
#' 
#' Update node and edge potentials of CRF model
#' 
#' This function updates \code{node.pot} and \code{edge.pot} of CRF model by using 
#' the current values of parameters and features.
#' 
#' There are two ways to model the relationship between parameters and features.
#' The first one exploits the special structure of features to reduce the memory
#' usage. However it may not suitable for all circumstances. The other one is more 
#' straighforward by explicitly specifying the coefficients of each parameter to
#' calculate the potentials, and may use much more memory. Two approaches can be
#' used together.
#' 
#' The first way uses the objects \code{node.par} and \code{edge.par} to define
#' the structure of features and provides the feature information in variables
#' \code{node.fea} and \code{edge.fea}. The second way directly provides the
#' feature information in variables \code{node.ext} and \code{edge.ext} without
#' any prior assumption on feature structure. \code{node.ext} is a list and
#' each element has the same structure as \code{node.pot}. \code{edge.ext} is
#' a list and each element has the same structure as \code{edge.pot}.
#' 
#' In detail, the node potential is updated as follows:
#' 
#' \deqn{
#' node.pot[n,i] = \sum_{f} par[node.par[n,i,f]] * node.fea[f,n] + \sum_{k} par[k] * node.ext[[k]][n,i]
#' }
#' 
#' and the edge potential is updated as follows:
#' 
#' \deqn{
#' edge.pot[[e]][i,j] = \sum_{f} par[edge.par[[e]][i,j,f]] * edge.fea[f,e] + \sum_{k} par[k] * edge.ext[[k]][[e]][i,j]
#' }
#' 
#' @param crf The CRF
#' @param node.fea The node features matrix with dimension \code{(n.nf, n.nodes)}
#' @param edge.fea The edge features matrix with dimension \code{(n.ef, n.edges)}
#' @param node.ext The extended information of node features
#' @param edge.ext The extended information of edge features
#' @return This function will directly modify the CRF and return the same CRF.
#' 
#' @seealso \code{\link{crf.nll}}, \code{\link{train.crf}}
#' 
#' @export
crf.update <- function(crf, node.fea = NULL, edge.fea = NULL, node.ext = NULL, edge.ext = NULL)
  .Call(CRF_Update, crf, node.fea, edge.fea, node.ext, edge.ext)



#' Calculate MRF sufficient statistics
#' 
#' Calculate the sufficient statistics of MRF model 
#' 
#' This function calculates the sufficient statistics of MRF model. This function
#' much be called before the first calling to \code{\link{mrf.nll}}. 
#' In the training data matrix \code{instances}, each row is an instance and 
#' each column corresponds a node in CRF.
#' 
#' @param crf The CRF
#' @param instances The training data matrix of MRF model
#' @return This function will return the value of MRF sufficient statistics.
#' 
#' @seealso \code{\link{mrf.nll}}, \code{\link{train.mrf}}
#' 
#' @export
mrf.stat <- function(crf, instances)
  .Call(MRF_Stat, crf, instances)



#' Calculate MRF negative log-likelihood
#' 
#' Calculate the negative log-likelihood of MRF model 
#' 
#' This function calculates the negative log-likelihood of MRF model as well as 
#' the gradient. This function is intended to be called by optimization algorithm
#' in training process. Before calling this function, the MRF sufficient 
#' statistics must be calculated and stored in object \code{par.stat} of CRF.
#' 
#' In the training data matrix \code{instances}, each row is an instance and 
#' each column corresponds a node in CRF.
#' 
#' @param crf The CRF
#' @param par The parameter vector of CRF
#' @param instances The training data matrix of MRF model
#' @param infer.method The inference method used to compute the likelihood
#' @param ... Other parameters need by the inference method
#' @return This function will return the value of MRF negative log-likilihood.
#' 
#' @seealso \code{\link{mrf.stat}}, \code{\link{mrf.update}}, \code{\link{train.mrf}}
#' 
#' @export
mrf.nll <- function(par, crf, instances, infer.method = infer.chain, ...)
  .Call(MRF_NLL, crf, par, instances, quote(infer.method(crf, ...)), environment())



#' Calculate CRF negative log likelihood
#' 
#' Calculate the negative log likelihood of CRF model
#' 
#' This function calculates the negative log likelihood of CRF model as well as
#' the gradient. This function is intended to be called by optimization algorithm
#' in training process.
#' 
#' In the training data matrix \code{instances}, each row is an instance and 
#' each column corresponds a node in CRF.
#' The variables \code{node.fea}, \code{edge.fea}, \code{node.ext}, \code{edge.ext}
#' are lists of length equal to the number of instances, and their elements are
#' defined as in \code{\link{crf.update}} respectively.
#' 
#' @param crf The CRF
#' @param par The parameter vector of CRF
#' @param instances The training data matrix of CRF model
#' @param node.fea The list of node features
#' @param edge.fea The list of edge features
#' @param node.ext The list of extended information of node features
#' @param edge.ext The list of extended information of edge features
#' @param infer.method The inference method used to compute the likelihood
#' @param ... Other parameters need by the inference method
#' @return This function will return the value of CRF negative log-likelihood.
#' 
#' @seealso \code{\link{crf.update}}, \code{\link{train.crf}}
#' 
#' @export
crf.nll <- function(par, crf, instances, node.fea = NULL, edge.fea = NULL, node.ext = NULL, edge.ext = NULL, infer.method = infer.chain, ...)
  .Call(CRF_NLL, crf, par, instances, node.fea, edge.fea, node.ext, edge.ext, quote(infer.method(crf, ...)), environment())



#' Train MRF model
#' 
#' Train the MRF model to estimate the parameters
#' 
#' This function trains the Markov Random Fields (MRF) model, which is a simple variant of CRF model.
#' 
#' In the training data matrix \code{instances}, each row is an instance and 
#' each column corresponds a node in CRF.
#' 
#' @param crf The CRF
#' @param instances The training data matrix of CRF model
#' @param nll The function to calculate negative log likelihood
#' @param trace Non-negative integer to control the tracing informtion of the optimization process
#' @return This function will directly modify the CRF and return the same CRF.
#' 
#' @seealso \code{\link{mrf.update}}, \code{\link{mrf.stat}}, \code{\link{mrf.nll}}, \code{\link{make.crf}}
#' 
#' @export
train.mrf <- function(crf, instances, nll = mrf.nll, trace = 0)
{
  gradient <- function(par, crf, ...) { crf$gradient }
	crf$par.stat <- mrf.stat(crf, instances)
	solution <- optim(crf$par, nll, gradient, crf, instances, method = "L-BFGS-B", control = list(trace = trace))
	crf$par <- solution$par
	mrf.update(crf)
	crf
}



#' Train CRF model
#' 
#' Train the CRF model to estimate the parameters
#' 
#' This function train the CRF model.
#' 
#' In the training data matrix \code{instances}, each row is an instance and 
#' each column corresponds a node in CRF.
#' The variables \code{node.fea}, \code{edge.fea}, \code{node.ext}, \code{edge.ext}
#' are lists of length equal to the number of instances, and their elements are
#' defined as in \code{\link{crf.update}} respectively.
#' 
#' @param crf The CRF
#' @param instances The training data matrix of CRF model
#' @param node.fea The list of node features
#' @param edge.fea The list of edge features
#' @param node.ext The list of extended information of node features
#' @param edge.ext The list of extended information of edge features
#' @param nll The function to calculate negative log likelihood
#' @param trace Non-negative integer to control the tracing informtion of the optimization process
#' @return This function will directly modify the CRF and return the same CRF.
#' 
#' @seealso \code{\link{crf.update}}, \code{\link{crf.nll}}, \code{\link{make.crf}}
#' 
#' @export
train.crf <- function(crf, instances, node.fea = NULL, edge.fea = NULL, node.ext = NULL, edge.ext = NULL, nll = crf.nll, trace = 0)
{
  gradient <- function(par, crf, ...) { crf$gradient }
  solution <- optim(crf$par, nll, gradient, crf, instances, node.fea, edge.fea, node.ext, edge.ext, method = "L-BFGS-B", control = list(trace = trace))
	crf$par <- solution$par
	crf.update(crf, node.fea[[1]], edge.fea[[1]], node.ext[[1]], edge.ext[[1]])
	crf
}
