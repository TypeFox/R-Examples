###############################################################################
###############################################################################

## BN generics

###############################################################################


#' learn a network (structure and parameters) of a \link{BN} from a \link{BNDataset}.
#' 
#' Learn a network (structure and parameters) of a \link{BN} from a \link{BNDataset} (see the \code{Details} section).
#' 
#' Learn the structure (the directed acyclic graph) of a \code{\link{BN}} object according to a \code{\link{BNDataset}}.
#' We provide three algorithms in order to learn the structure of the network, that can be chosen with the \code{algo} parameter.
#' The first is the Silander-Myllym\"aki (\code{sm})
#' exact search-and-score algorithm, that performs a complete evaluation of the search space in order to discover
#' the best network; this algorithm may take a very long time, and can be inapplicable when discovering networks
#' with more than 25--30 nodes. Even for small networks, users are strongly encouraged to provide
#' meaningful parameters such as the layering of the nodes, or the maximum number of parents -- refer to the 
#' documentation in package manual for more details on the method parameters.
#' 
#' The second algorithm (and the default one) is the Max-Min Hill-Climbing heuristic (\code{mmhc}), that performs a statistical
#' sieving of the search space followed by a greedy evaluation. It is considerably faster than the complete method,
#' at the cost of a (likely)
#' lower quality. Also note that in the case of a very dense network and lots of obsevations, the statistical evaluation
#' of the search space may take a long time. Also for this algorithm there are parameters that may need to be tuned,
#' mainly the confidence threshold of the statistical pruning.
#' 
#' The third method is the Structural Expectation-Maximization (\code{sem}) algorithm,
#' for learning a network from a dataset with missing values. It iterates a sequence of Expectation-Maximization (in order to ``fill in''
#' the holes in the dataset) and structure learning from the guessed dataset, until convergence. The structure learning used inside SEM,
#' due to computational reasons, is MMHC. Convergence of SEM can be controlled with the parameters \code{struct.threshold}
#' and \code{param.threshold}, for the structure and the parameter convergence, respectively.
#' 
#' Search-and-score methods also need a scoring function to compute an estimated measure of each configuration of nodes.
#' We provide three of the most popular scoring functions, \code{BDeu} (Bayesian-Dirichlet equivalent uniform, default),
#' \code{AIC} (Akaike Information Criterion) and \code{BIC} (Bayesian Information Criterion). The scoring function
#' can be chosen using the \code{scoring.func} parameter.
#' 
#' Then, the parameters of the network are learnt using MAP (Maximum A Posteriori) estimation (if not using bootstrap).
#' 
#' See documentation for \code{\link{learn.structure}} and \code{\link{learn.params}} for more informations.
#' 
#' @name learn.network
#' @rdname learn.network
#' 
#' @param x can be a \code{\link{BN}} or a \code{\link{BNDataset}}. If \code{x} is a \code{\link{BN}},
#' then also the \code{dataset} parameter must be given.
#' @param y a \code{\link{BNDataset}} object, to be provided only if \code{x} is a \code{\link{BN}}.
#' @param algo the algorithm to use. Currently, one among:
#'        \code{sm} (Silander-Myllymaki),
#'        \code{mmhc} (Max-Min Hill Climbing, default) and
#'        \code{sem} (Structural Expectation Maximization).
#' @param scoring.func the scoring function to use. Currently, one among
#'        \code{BDeu}, \code{AIC}, \code{BIC}.
#' @param initial.network network srtructure to be used as starting point for structure search.
#'        Can take different values:
#'        a \code{BN} object, a matrix containing the adjacency matrix of the structure of the network,
#'        or the string \code{random.chain} to sample a random chain as starting point.
#' @param alpha confidence threshold (only for \code{mmhc}).
#' @param ess Equivalent Sample Size value.
#' @param bootstrap \code{TRUE} to use bootstrap samples. 
#' @param layering vector containing the layers each node belongs to.
#' @param max.fanin.layers matrix of available parents in each layer (only for \code{sm}).
#' @param max.fanin maximum number of parents for each node (only for \code{sm}).
#' @param layer.struct \code{0/1} matrix for indicating which layers can contain parent nodes
#'        for nodes in a layer (only for \code{mmhc}).
#' @param cont.nodes vector containing the index of continuous variables.
#' @param use.imputed.data \code{TRUE} to learn the structure from the imputed dataset
#' (if available, a check is performed). Default is to use raw dataset
#' @param use.cpc (when using \code{mmhc}) compute Candidate Parent-and-Children sets instead of 
#' starting the Hill Climbing from an empty graph.
#' @param ... potential further arguments for methods.
#' 
#' @return new \code{\link{BN}} object with structure (DAG) and conditional probabilities
#' as learnt from the given dataset.
#' 
#' @seealso learn.structure learn.params
#' 
#' @examples
#' \dontrun{
#' mydataset <- BNDataset("data.file", "header.file")
#' 
#' # starting from a BN
#' net <- BN(mydataset)
#' net <- learn.network(net, mydataset)
#' 
#' # start directly from the dataset
#' net <- learn.network(mydataset)
#' }
#' 
#' @exportMethod learn.network
setGeneric("learn.network", function(x, ...)#dataset, algo="mmhc", scoring.func="BDeu", alpha=0.05, ess=1, bootstrap=FALSE,
                                     #layering=c(), max.fanin.layers=NULL, max.fanin=num.variables(dataset),
                                     #layer.struct = NULL,
                                     #cont.nodes=c(), use.imputed.data=FALSE, use.cpc=TRUE, ...)
  standardGeneric("learn.network"))


#' learn the parameters of a \link{BN}.
#' 
#' Learn the parameters of a \link{BN} object according to a \link{BNDataset}
#' using MAP (Maximum A Posteriori) estimation.
#' 
#' @name learn.params
#' @rdname learn.params
#' 
#' @param bn a \code{\link{BN}} object.
#' @param dataset a \code{\link{BNDataset}} object.
#' @param ess Equivalent Sample Size value.
#' @param use.imputed.data use imputed data.
#' 
#' @return new \code{\link{BN}} object with conditional probabilities.
#' 
#' @seealso learn.network
#' 
#' @examples
#' \dontrun{
#' ## first create a BN and learn its structure from a dataset
#' dataset <- BNDataset("file.header", "file.data")
#' bn <- BN(dataset)
#' bn <- learn.structure(bn, dataset)
#' bn <- learn.params(bn, dataset, ess=1)
#' }
#' 
#' @exportMethod learn.params
setGeneric("learn.params", function(bn, dataset, ess=1, use.imputed.data=F) standardGeneric("learn.params"))


#' learn the structure of a network.
#' 
#' Learn the structure (the directed acyclic graph) of a \code{\link{BN}} object according to a \code{\link{BNDataset}}.
#'
#' We provide three algorithms in order to learn the structure of the network, that can be chosen with the \code{algo} parameter.
#' The first is the Silander-Myllym\"aki (\code{sm})
#' exact search-and-score algorithm, that performs a complete evaluation of the search space in order to discover
#' the best network; this algorithm may take a very long time, and can be inapplicable when discovering networks
#' with more than 25--30 nodes. Even for small networks, users are strongly encouraged to provide
#' meaningful parameters such as the layering of the nodes, or the maximum number of parents -- refer to the 
#' documentation in package manual for more details on the method parameters.
#' 
#' The second algorithm (and the default one) is the Max-Min Hill-Climbing heuristic (\code{mmhc}), that performs a statistical
#' sieving of the search space followed by a greedy evaluation. It is considerably faster than the complete method, at the cost of a (likely)
#' lower quality. Also note that in the case of a very dense network and lots of obsevations, the statistical evaluation
#' of the search space may take a long time. Also for this algorithm there are parameters that may need to be tuned,
#' mainly the confidence threshold of the statistical pruning.
#' 
#' The third method is the Structural Expectation-Maximization (\code{sem}) algorithm,
#' for learning a network from a dataset with missing values. It iterates a sequence of Expectation-Maximization (in order to ``fill in''
#' the holes in the dataset) and structure learning from the guessed dataset, until convergence. The structure learning used inside SEM,
#' due to computational reasons, is MMHC. Convergence of SEM can be controlled with the parameters \code{struct.threshold}
#' and \code{param.threshold}, for the structure and the parameter convergence, respectively.
#' 
#' Search-and-score methods also need a scoring function to compute an estimated measure of each configuration of nodes.
#' We provide three of the most popular scoring functions, \code{BDeu} (Bayesian-Dirichlet equivalent uniform, default),
#' \code{AIC} (Akaike Information Criterion) and \code{BIC} (Bayesian Information Criterion). The scoring function
#' can be chosen using the \code{scoring.func} parameter.
#' 
#' @name learn.structure
#' @rdname learn.structure
#' 
#' @param bn a \code{\link{BN}} object.
#' @param dataset a \code{\link{BNDataset}}.
#' @param algo the algorithm to use. Currently, one among \code{sm} (Silander-Myllymaki), \code{mmhc}
#'        (Max-Min Hill Climbing, default) and \code{sem} (Structural Expectation Maximization).
#' @param scoring.func the scoring function to use. Currently, one among \code{BDeu}, \code{AIC}, \code{BIC}.
#' @param initial.network network srtructure to be used as starting point for structure search.
#'        Can take different values:
#'        a \code{BN} object, a matrix containing the adjacency matrix of the structure of the network,
#'        or the string \code{random.chain} to sample a random chain as starting point.
#' @param alpha confidence threshold (only for \code{mmhc}).
#' @param ess Equivalent Sample Size value.
#' @param bootstrap \code{TRUE} to use bootstrap samples. 
#' @param layering vector containing the layers each node belongs to (only for \code{sm}).
#' @param max.fanin.layers matrix of available parents in each layer (only for \code{sm}).
#' @param max.fanin maximum number of parents for each node (only for \code{sm}).
#' @param layer.struct prior knowledge for layering structure (only for \code{mmhc}).
#' @param cont.nodes vector containing the index of continuous variables.
#' @param use.imputed.data \code{TRUE} to learn the structure from the imputed dataset
#' (if available, a check is performed). Default is to use raw dataset
#' @param use.cpc (when using \code{mmhc}) compute Candidate Parent-and-Children sets instead of 
#' starting the Hill Climbing from an empty graph.
#' @param ... potential further arguments for method.
#' 
#' @return new \code{\link{BN}} object with DAG.
#' 
#' @seealso learn.network
#' 
#' @examples
#' \dontrun{
#' dataset <- BNDataset("file.header", "file.data")
#' bn <- BN(dataset)
#' # use MMHC
#' bn <- learn.structure(bn, dataset, alpha=0.05, ess=1, bootstrap=FALSE)
#' 
#' # now use Silander-Myllymaki
#' layers <- layering(bn)
#' mfl <- as.matrix(read.table(header=F,
#' text='0 1 1 1 1 0 1 1 1 1 0 0 8 7 7 0 0 0 14 6 0 0 0 0 19'))
#' bn <- learn.structure(bn, dataset, algo='sm', max.fanin=3, cont.nodes=c(),
#'                       layering=layers, max.fanin.layers=mfl, use.imputed.data=FALSE)
#' }
#' 
#' @exportMethod learn.structure
setGeneric("learn.structure", function(bn, dataset, algo="mmhc", scoring.func="BDeu", initial.network=NULL,
                                       alpha=0.05, ess=1, bootstrap=FALSE,
                                       layering=c(), max.fanin.layers=NULL, max.fanin=num.variables(dataset),
                                       layer.struct = NULL,
                                       cont.nodes=c(), use.imputed.data=FALSE, use.cpc=TRUE, ...) standardGeneric("learn.structure"))


#' return the layering of the nodes.
#' 
#' Compute the topological ordering of the nodes of a network, in order to divide the network in layers.
#' 
#' @name layering
#' @rdname layering
#' 
#' @param x a \code{\link{BN}} object.
#' 
#' @examples
#' \dontrun{
#' dataset <- BNDataset("file.header", "file.data")
#' x <- BN(dataset)
#' x <- learn.network(x, dataset)
#' layering(x)
#' }
#' 
#' @return a vector containing layers the nodes can be divided into.
#' 
#' @exportMethod layering
setGeneric("layering", function(x) standardGeneric("layering"))


#' compute the most probable values to be observed.
#' 
#' Return an array containing the values that each variable of the network is more likely to take, according to the CPTS.
#' In case of ties take the first value.
#' 
#' @name get.most.probable.values
#' @rdname get.most.probable.values
#' 
#' @param x a \code{\link{BN}} or \code{\link{InferenceEngine}} object.
#' 
#' @return array containing, in each position, the most probable value for the corresponding variable.
#' 
#' @examples
#' \dontrun{
#' # try with a BN object x
#' get.most.probable.values(x)
#' 
#' # now build an InferenceEngine object
#' eng <- InferenceEngine(x)
#' get.most.probable.values(eng)
#' }
#'  
#' @exportMethod get.most.probable.values
setGeneric("get.most.probable.values", function(x) standardGeneric("get.most.probable.values"))


#' sample a row vector of values for a network.
#' 
#' @name sample.row
#' @rdname sample.row
#' 
#' @param x a \code{\link{BN}} or \code{\link{InferenceEngine}} object.
#' 
#' @return a vector of values.
#' 
#' @exportMethod sample.row
setGeneric("sample.row", function(x) standardGeneric("sample.row"))


#' sample a \code{\link{BNDataset}} from a network of an inference engine.
#' 
#' @name sample.dataset
#' @rdname sample.dataset
#' 
#' @param x a \code{\link{BN}} or \code{\link{InferenceEngine}} object.
#' @param n number of items to sample.
#' 
#' @return a \code{\link{BNDataset}}
#' 
#' @exportMethod sample.dataset
setGeneric("sample.dataset", function(x, n=100) standardGeneric("sample.dataset"))


#' compute the list of inferred marginals of a BN.
#' 
#' Given an \code{\link{InferenceEngine}}, it returns a list containing the marginals for the variables
#' in the network, according to the propagated beliefs.
#' 
#' @name marginals
#' @rdname marginals
#' 
#' @param x an \code{\link{InferenceEngine}}
#' @param ... potential further arguments of methods.
#' 
#' @return a list containing the marginals of each variable, as probability tables.
#' 
#' @examples
#' \dontrun{
#' eng <- InferenceEngine(net)
#' marginals(eng)
#' }
#' 
#' @exportMethod marginals
setGeneric("marginals", function(x, ...) standardGeneric("marginals"))


# ' query BN given observations
# ' 
# ' @name query
# ' @rdname query
# ' 
# ' @param x a BN.
# ' @param observed.vars vector of observed variables.
# ' @param observed.vals vector of observed values for corresponding variables in \code{observed.vars}.
# ' 
# ' @return most probable values given observations
# ' 
# ' @exportMethod query
# setGeneric("query", function(x, observed.vars=c(), observed.vals=c()) standardGeneric("query"))


#' save a \code{\link{BN}} picture as \code{.eps} file.
#' 
#' Save an image of a Bayesian Network as an \code{.eps} file.
#' 
#' @name save.to.eps
#' @rdname save.to.eps
#' 
#' @param x a \code{\link{BN}} object
#' @param filename name (with path, if needed) of the file to be created
#' 
#' @examples
#' \dontrun{
#' save.to.eps(x, "out.eps")
#' }
#' 
#' @seealso \code{\link{plot}}
#' 
#' @exportMethod save.to.eps
setGeneric("save.to.eps", function(x, filename) standardGeneric("save.to.eps"))


#' Read a network from a \code{.dsc} file.
#' 
#' Read a network described in a \code{.dsc}-formatted file, and
#' build a \code{\link{BN}} object.
#' 
#' The method relies on a coherent ordering of variable values and parameters in the file.
#' 
#' @name read.dsc
#' @rdname read.dsc
#' 
#' @param x the \code{.dsc} file, with absolute/relative position.
#' 
#' @return a \code{\link{BN}} object.
#' 
#' @exportMethod read.dsc
setGeneric("read.dsc", function(x) standardGeneric("read.dsc"))


#' Read a network from a \code{.bif} file.
#' 
#' Read a network described in a \code{.bif}-formatted file, and
#' build a \code{\link{BN}} object.
#' 
#' The method relies on a coherent ordering of variable values and parameters in the file.
#' 
#' @name read.bif
#' @rdname read.bif
#' 
#' @param x the \code{.bif} file, with absolute/relative position.
#' 
#' @return a \code{\link{BN}} object.
#' 
#' @exportMethod read.bif
setGeneric("read.bif", function(x) standardGeneric("read.bif"))


#' Read a network from a \code{.net} file.
#' 
#' Read a network described in a \code{.net}-formatted file, and
#' build a \code{\link{BN}} object.
#' 
#' The method relies on a coherent ordering of variable values and parameters in the file.
#' 
#' @name read.net
#' @rdname read.net
#' 
#' @param x the \code{.net} file, with absolute/relative position.
#' 
#' @return a \code{\link{BN}} object.
#' 
#' @exportMethod read.net
setGeneric("read.net", function(x) standardGeneric("read.net"))


#' Write a network saving it in a \code{.dsc} file.
#' 
#' Write a network on disk, saving it in a \code{.dsc}-formatted file.
#' 
#' @name write.dsc
#' @rdname write.dsc
#' 
#' @param x the \code{\link{BN}} object.
#' @param path the relative or absolute path of the directory of the created file.
#' 
#' @exportMethod write.dsc
setGeneric("write.dsc", function(x, path="./") standardGeneric("write.dsc"))


#' Read the scoring function used to learn the structure of a network.
#' 
#' Read the scoring function used in the \code{\link{learn.structure}} method.
#' Outcome is meaningful only if the structure of a network has been learnt.
#' 
#' @name scoring.func
#' @rdname scoring.func
#' 
#' @param x the \code{\link{BN}} object.
#' 
#' @return the scoring function used.
#' 
#' @exportMethod scoring.func
setGeneric("scoring.func", function(x) standardGeneric("scoring.func"))


#' Set the scoring function used to learn the structure of a network.
#' 
#' Set the scoring function used in the \code{\link{learn.structure}} method.
#' 
#' @name scoring.func<-
#' @rdname scoring.func-set
#' 
#' @param x the \code{\link{BN}} object.
#' @param value the scoring function used.
#' 
#' @return updated BN.
#' 
#' @exportMethod scoring.func<-
setGeneric("scoring.func<-", function(x, value) standardGeneric("scoring.func<-"))


#' Read the algorithm used to learn the structure of a network.
#' 
#' Read the algorithm used in the \code{\link{learn.structure}} method.
#' Outcome is meaningful only if the structure of a network has been learnt.
#' 
#' @name struct.algo
#' @rdname struct.algo
#' 
#' @param x the \code{\link{BN}} object.
#' 
#' @return the structure learning algorithm used.
#' 
#' @exportMethod struct.algo
setGeneric("struct.algo", function(x) standardGeneric("struct.algo"))


#' Set the algorithm used to learn the structure of a network.
#' 
#' Set the algorithm used in the \code{\link{learn.structure}} method.
#' 
#' @name struct.algo<-
#' @rdname struct.algo-set
#' 
#' @param x the \code{\link{BN}} object.
#' @param value the scoring function used.
#' 
#' @return updated BN.
#' 
#' @exportMethod struct.algo<-
setGeneric("struct.algo<-", function(x, value) standardGeneric("struct.algo<-"))


#' Initialize a WPDAG from a DAG.
#' 
#' Given a \code{\link{BN}} object with a \code{dag}, return a network
#' with its \code{wpdag} set as the CPDAG computed starting from the \code{dag}.
#' 
#' @name wpdag.from.dag
#' @rdname wpdag.from.dag
#' 
#' @param x a \code{\link{BN}} object.
#' @param layering vector containing the layers each node belongs to.
#' 
#' @return a \code{\link{BN}} object with an initialized \code{wpdag}.
#' 
#' @examples
#' \dontrun{
#' net <- learn.network(dataset, layering=layering)
#' wp.net <- wpdag.from.dag(net, layering)
#' }
#' 
#' @seealso \code{\link{dag.to.cpdag}}
#' 
#' @exportMethod wpdag.from.dag
setGeneric("wpdag.from.dag", function(x, layering=NULL) standardGeneric("wpdag.from.dag"))


###############################################################################
###############################################################################

## BNDataset generics

###############################################################################


#' check if a BNDataset contains raw data.
#' 
#' Check whether a \code{\link{BNDataset}} object actually contains raw data.
#' 
#' @name has.raw.data
#' @rdname has.raw.data
#' 
#' @param x a \code{\link{BNDataset}}.
#' 
#' @examples
#' \dontrun{
#' x <- BNDataset()
#' has.raw.data(x) # FALSE
#' 
#' x <- read.dataset(x, "file.header", "file.data")
#' has.raw.data(x) # TRUE, since read.dataset() actually reads raw data.
#' }
#' 
#' @seealso \code{\link{has.imputed.data}}, \code{\link{raw.data}}, \code{\link{imputed.data}}
#' 
#' @exportMethod has.raw.data
setGeneric("has.raw.data", function(x) standardGeneric("has.raw.data"))

#' check if a BNDataset contains impited data.
#' 
#' Check whether a \code{\link{BNDataset}} object actually contains imputed data.
#' 
#' @name has.imputed.data
#' @rdname has.imputed.data
#' 
#' @param x a \code{\link{BNDataset}}.
#' 
#' @examples
#' \dontrun{
#' x <- BNDataset()
#' has.imputed.data(x) # FALSE
#' 
#' x <- read.dataset(x, "file.header", "file.data")
#' has.imputed.data(x) # FALSE, since read.dataset() actually reads raw data.
#' 
#' x <- impute(x)
#' has.imputed.data(x) # TRUE
#' }
#' 
#' @seealso \code{\link{has.raw.data}}, \code{\link{raw.data}}, \code{\link{imputed.data}}
#' 
#' @exportMethod has.imputed.data
setGeneric("has.imputed.data", function(x) standardGeneric("has.imputed.data"))


#' get raw data of a BNDataset.
#' 
#' Return raw data contained in a \code{\link{BNDataset}} object, if any.
#' 
#' @name raw.data
#' @rdname raw.data
#' 
#' @param x a \code{\link{BNDataset}}.
#' 
#' @seealso \code{\link{has.raw.data}}, \code{\link{has.imputed.data}}
#' 
#' @exportMethod raw.data
setGeneric("raw.data", function(x) standardGeneric("raw.data"))


#' get imputed data of a BNDataset.
#' 
#' Return imputed data contained in a \code{\link{BNDataset}} object, if any.
#' 
#' @name imputed.data
#' @rdname imputed.data
#' 
#' @param x a \code{\link{BNDataset}}.
#' 
#' @seealso \code{\link{has.raw.data}}, \code{\link{has.imputed.data}}, \code{\link{raw.data}}
#' 
#' @exportMethod imputed.data
setGeneric("imputed.data", function(x) standardGeneric("imputed.data"))


#' add raw data.
#' 
#' Insert raw data in a \code{\link{BNDataset}} object.
#' 
#' @name raw.data<-
#' @rdname raw.data-set
#' 
#' @param x a \code{\link{BNDataset}}.
#' @param value a matrix of integers containing a dataset.
#' 
#' @seealso \code{\link{has.raw.data}}, \code{\link{raw.data}}, \code{\link{read.dataset}}
#' 
#' @exportMethod raw.data<-
setGeneric("raw.data<-", function(x, value) standardGeneric("raw.data<-"))


#' add imputed data.
#' 
#' Insert imputed data in a \code{\link{BNDataset}} object.
#' 
#' @name imputed.data<-
#' @rdname imputed.data-set
#' 
#' @param x a \code{\link{BNDataset}}.
#' @param value a matrix of integers containing a dataset.
#' 
#' @seealso \code{\link{has.imputed.data}}, \code{\link{imputed.data}}, \code{\link{read.dataset}}
#' 
#' @exportMethod imputed.data<-
setGeneric("imputed.data<-", function(x, value) standardGeneric("imputed.data<-"))


#' Subset a \code{\link{BNDataset}} to get only complete cases.
#' 
#' Given a \code{\link{BNDataset}}, return a copy of the original object where
#' the \code{raw.data} consists only in the observations that do not contain missing values.
#' 
#' Non-missingness can be required on a subset of variables (by default, on all variables).
#' 
#' If present, imputed data and bootstrap samples are eliminated from the
#' new \code{\link{BNDataset}}, as using this method *after* using \code{\link{impute}}
#' or \code{\link{bootstrap}}, there may likely be a loss of correspondence between
#' the subsetted \code{raw.data} and the previously generated \code{imputed.data}
#' and \code{bootstrap} samples.
#' 
#' @name complete
#' @rdname complete
#' 
#' @param x a \code{\link{BNDataset}}.
#' @param complete.vars vector containing the indices of the variables to be considered
#' for the subsetting; variables not included in the vector can still contain \code{NA}s.
#' 
#' @return a copy of the original \code{\link{BNDataset}} containing only complete observations.
#' 
#' @exportMethod complete
setGeneric("complete", function(x, complete.vars=seq_len(num.variables(x))) standardGeneric("complete"))


#' Read a dataset from file.
#' 
#' There are two ways to build a BNDataset: using two files containing respectively header informations
#' and data, and manually providing the data table and the related header informations
#' (variable names, cardinality and discreteness).
#' 
#' The key informations needed are:
#' 1. the data;
#' 2. the state of variables (discrete or continuous);
#' 3. the names of the variables;
#' 4. the cardinalities of the variables (if discrete), or the number of levels they have to be quantized into
#' (if continuous). 
#' Names and cardinalities/leves can be guessed by looking at the data, but it is strongly advised to provide
#' _all_ of the informations, in order to avoid problems later on during the execution.
#' 
#' Data can be provided in form of data.frame or matrix. It can contain NAs. By default, NAs are indicated with '?';
#' to specify a different character for NAs, it is possible to provide also the \code{na.string.symbol} parameter.
#' The values contained in the data have to be numeric (real for continuous variables, integer for discrete ones).
#' The default range of values for a discrete variable \code{X} is \code{[1,|X|]}, with \code{|X|} being
#' the cardinality of \code{X}. The same applies for the levels of quantization for continuous variables.
#' If the value ranges for the data are different from the expected ones, it is possible to specify a different
#' starting value (for the whole dataset) with the \code{starts.from} parameter. E.g. by \code{starts.from=0}
#' we assume that the values of the variables in the dataset have range \code{[0,|X|-1]}.
#' Please keep in mind that the internal representation of bnstruct starts from 1,
#' and the original starting values are then lost. 
#' 
#' It is possible to use two files, one for the data and one for the metadata,
#' instead of providing manually all of the info. 
#' bnstruct requires the data files to be in a format subsequently described.
#' The actual data has to be in (a text file containing data in) tabular format, one tuple per row,
#' with the values for each variable separated by a space or a tab. Values for each variable have to be
#' numbers, starting from \code{1} in case of discrete variables.
#' Data files can have a first row containing the names of the corresponding variables.
#' 
#' In addition to the data file, a header file containing additional informations can also be provided.
#' An header file has to be composed by three rows of tab-delimited values:
#' 1. list of names of the variables, in the same order of the data file;
#' 2. a list of integers representing the cardinality of the variables, in case of discrete variables,
#'   or the number of levels each variable has to be quantized in, in case of continuous variables;
#' 3. a list that indicates, for each variable, if the variable is continuous
#'   (\code{c} or \code{C}), and thus has to be quantized before learning,
#'   or discrete (\code{d} or \code{D}).
#' 
#' @name read.dataset
#' @rdname read.dataset
#' 
#' @param object the \code{\link{BNDataset}} object.
#' @param data.file the \code{data} file.
#' @param header.file the \code{header} file.
#' @param data.with.header \code{TRUE} if the first row of \code{dataset} file is an header (e.g. it contains the variable names).
#' @param na.string.symbol character that denotes \code{NA} in the dataset.
#' @param sep.symbol separator among values in the dataset.
#' @param starts.from starting value for entries in the dataset (observed values, default is 1).
#' 
#' @seealso BNDataset
#' 
#' @examples
#' \dontrun{
#' dataset <- BNDataset()
#' dataset <- read.dataset(dataset, "file.data", "file.header")
#' }
#' 
#' @exportMethod read.dataset
setGeneric("read.dataset", function(object, data.file, header.file, data.with.header = FALSE,
                                    na.string.symbol = '?', sep.symbol = '', starts.from = 1)
                            standardGeneric("read.dataset"))


#' Impute a \code{\link{BNDataset}} raw data with missing values.
#'
#' @name impute
#' @rdname impute
#' 
#' @param object the \code{\link{BNDataset}} object.
#' @param k.impute number of neighbours to be used; for discrete variables we use mode,
#' for continuous variables the median value is instead taken.
#' 
#' @examples
#' \dontrun{
#' dataset <- BNDataset("file.data", "file.header")
#' dataset <- impute(dataset)
#' }
#' 
#' @exportMethod impute
setGeneric("impute", function(object, k.impute=10) standardGeneric("impute"))


#' Perform bootstrap.
#' 
#' Create a list of \code{num.boots} samples of the original dataset.
#' 
#' @name bootstrap
#' @rdname bootstrap
#' 
#' @param object the \code{\link{BNDataset}} object.
#' @param num.boots number of sampled datasets for bootstrap.
#' @param seed random seed.
#' @param imputation \code{TRUE} if imputation has to be performed. Default is \code{FALSE}.
#' @param k.impute number of neighbours to be used; for discrete variables we use mode, for continuous variables the median value is instead taken (useful only if imputation == TRUE).
#' 
#' @examples
#' \dontrun{
#' dataset <- BNDataset("file.data", "file.header")
#' dataset <- bootstrap(dataset, num.boots = 1000)
#' }
#' 
#' @exportMethod bootstrap
setGeneric("bootstrap", function(object, num.boots = 100, seed = 0, imputation = FALSE, k.impute = 10)
                             standardGeneric("bootstrap"))


#' get selected element of bootstrap list.
#' 
#' Given a \code{\link{BNDataset}}, return the sample corresponding to given index.
#' 
#' @name boot
#' @rdname boot
#' 
#' @param dataset a \code{\link{BNDataset}} object.
#' @param index the index of the requested sample.
#' @param use.imputed.data \code{TRUE} if samples from imputed dataset are to be used. Default if \code{FALSE}.
#' 
#' @seealso bootstrap
#' 
#' @examples
#' \dontrun{
#' dataset <- BNDataset("file.data", "file.header")
#' dataset <- bootstrap(dataset, num.boots = 1000)
#' 
#' for (i in 1:num.boots(dataset))
#'    print(boot(dataset, i))
#' }
#' 
#' @seealso \code{\link{bootstrap}}
#' 
#' @exportMethod boot
setGeneric("boot", function(dataset, index, use.imputed.data = FALSE) standardGeneric("boot"))


###############################################################################
###############################################################################

## InferenceEngine generics

###############################################################################

#' build a JunctionTree.
#' 
#' Starting from the adjacency matrix of the directed acyclic graph of the network
#' contained in an InferenceEngine, build a JunctionTree for the network
#' and store it into an InferenceEngine.
#' 
#' @name build.junction.tree
#' @rdname build.junction.tree
#' 
#' @param object an \code{\link{InferenceEngine}} object.
#' @param ... potential further arguments for methods.
#' 
#' @seealso InferenceEngine
#' 
#' @examples
#' \dontrun{
#' dataset <- BNDataset("file.header", "file.data")
#' net <- BN(dataset)
#' eng <- InferenceEngine()
#' eng <- build.junction.tree(eng)
#' }
#' 
#' @exportMethod build.junction.tree
setGeneric("build.junction.tree", function(object, ...) standardGeneric("build.junction.tree"))


#' perform belief propagation.
#' 
#' Perform belief propagation for the network of an InferenceEngine, given a set of observations when present.
#' In the current version of \code{bnstruct}, belief propagation can be computed only over a junction tree.
#' 
#' @name belief.propagation
#' @rdname belief.propagation
#' 
#' @param ie an \code{\link{InferenceEngine}} object.
#' @param observations list of observations, consisting in two vector, \code{observed.vars} for the observed variables,
#' and \code{observed.vals} for the values taken by variables listed in \code{observed.vars}. If no observations
#' are provided, the \code{InferenceEngine} will use the ones it already contains.
#' @param return.potentials if TRUE only the potentials are returned, instead of the default \code{\link{BN}}.
#' 
#' @return updated \code{\link{InferenceEngine}} object.
#' 
#' @examples
#' \dontrun{
#' dataset <- BNDataset("file.header", "file.data")
#' bn <- BN(dataset)
#' ie <- InferenceEngine(bn)
#' ie <- belief.propagation(ie)
#' 
#' observations(ie) <- list("observed.vars"=("A","G","X"), "observed.vals"=c(1,2,1))
#' belief.propagation(ie)
#' }
#' 
#' @exportMethod belief.propagation
setGeneric("belief.propagation", function(ie, observations = NULL,
                                          return.potentials = FALSE) standardGeneric("belief.propagation"))


#' check if an updated \code{\link{BN}} is present in an \code{\link{InferenceEngine}}.
#' 
#' Check if an InferenceEngine actually contains an updated network, in order to provide the chance of
#' a fallback and use the original network if no belief propagation has been performed.
#' 
#' @name test.updated.bn
#' @rdname test.updated.bn
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' 
#' @return \code{TRUE} if an updated network is contained in the InferenceEngine, \code{FALSE} otherwise.
#' 
#' @examples
#' \dontrun{
#' dataset <- BNDataset("file.header", "file.data")
#' bn <- BN(dataset)
#' ie <- InferenceEngine(bn)
#' test.updated.bn(ie) # FALSE
#' 
#' observations(ie) <- list("observed.vars"=("A","G","X"), "observed.vals"=c(1,2,1))
#' ie <- belief.propagation(ie)
#' test.updated.bn(ie) # TRUE
#' }
#' 
#' @exportMethod test.updated.bn
setGeneric("test.updated.bn", function(x) standardGeneric("test.updated.bn"))


#' expectation-maximization algorithm.
#' 
#' Learn parameters of a network using the Expectation-Maximization algorithm.
#' 
#' @name em
#' @rdname em
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' @param dataset observed dataset with missing values for the Bayesian Network of \code{x}.
#' @param threshold threshold for convergence, used as stopping criterion.
#' @param max.em.iterations maximum number of iterations to run in case of no convergence.
#' @param ess Equivalent Sample Size value.
#' 
#' @return a list containing: an \code{\link{InferenceEngine}} with a new updated network (\code{"InferenceEngine"}),
#'         and the imputed dataset (\code{"BNDataset"}).
#' 
#' @examples
#' \dontrun{
#' em(x, dataset)
#' }
#' 
#' @exportMethod em
setGeneric("em", function(x, dataset, threshold = 0.001,
                          max.em.iterations = 10, ess = 1) standardGeneric("em"))


# ' Structural Expectation-Maximization algorithm.
# ' 
# ' Learn structure and parameters of a network with the Structural EM algorithm.
# ' 
# ' @name sem
# ' @rdname sem
# ' 
# ' @param x a \code{\link{BN}} object.
# ' @param dataset observed dataset with missing values for the Bayesian Network of \code{x}.
# ' @param struct.threshold threshold for convergence of the structure learning step, used as stopping criterion.
# ' @param param.threshold threshold for convergence of the parameter learning step, used as stopping criterion.
# ' @param max.em.iterations maximum number of iterations of EM to run in case of no convergence.
# ' @param max.sem.iterations maximum number of iterations of SEM to run in case of no convergence.
# ' @param scoring.func the scoring function to use. Currently, one among \code{AIC} and \code{BIC}
# ' (default - \code{BDeu} supported  as linear approximation).
# ' @param alpha confidence threshold (only for \code{mmhc}).
# ' @param ess Equivalent Sample Size value.
# ' @param bootstrap \code{TRUE} to use bootstrap samples. 
# ' @param layering vector containing the layers each node belongs to (only for \code{sm}).
# ' @param max.fanin.layers matrix of available parents in each layer (only for \code{sm}).
# ' @param max.fanin maximum number of parents for each node (only for \code{sm}).
# ' @param cont.nodes vector containing the index of continuous variables.
# ' @param use.imputed.data \code{TRUE} to learn the structure from the imputed dataset
# ' (if available, a check is performed). Default is to use raw dataset
# ' @param use.cpc (when using \code{mmhc}) compute Candidate Parent-and-Children sets instead of 
# ' starting the Hill Climbing from an empty graph.
# ' @param ... further potential arguments for method.
# ' 
# ' @return a (\code{"BN"}) network with the new structure.
# ' 
# exportMethod sem
setGeneric("sem", function(x, dataset, struct.threshold = 0, param.threshold = 0, max.sem.iterations = 25,
                           max.em.iterations = 10, scoring.func = "BDeu", initial.network = NULL,
                           alpha = 0.05, ess = 1, bootstrap = FALSE,
                           layering = c(), max.fanin.layers = NULL,
                           max.fanin = num.variables(dataset), cont.nodes = c(), use.imputed.data = FALSE,
                           use.cpc = T, ...) standardGeneric("sem"))


###############################################################################
###############################################################################

## Accessors and mutators

###############################################################################

#' get name of an object.
#' 
#' Return the name of an object, of class \code{\link{BN}} or \code{\link{BNDataset}}.
#' 
#' @name name
#' @rdname name
#' 
#' @param x an object.
#' 
#' @return name of the object.
#' 
#' @exportMethod name
setGeneric("name", function(x) standardGeneric("name"))


#' get number of nodes of an object.
#' 
#' Return the name of an object, of class \code{\link{BN}} or \code{\link{InferenceEngine}}.
#' 
#' @name num.nodes
#' @rdname num.nodes
#' 
#' @param x an object.
#' 
#' @return number of nodes of the desired object.
#' 
#' @exportMethod num.nodes
setGeneric("num.nodes", function(x) standardGeneric("num.nodes"))


#' get variables of an object.
#' 
#' Get the list of variables (with their names) of a \code{\link{BN}} or \code{\link{BNDataset}}.
#' 
#' @name variables
#' @rdname variables
#' 
#' @param x an object.
#' 
#' @return vector of the variables names of the desired object.
#' 
#' @exportMethod variables
setGeneric("variables", function(x) standardGeneric("variables"))


#' get status (discrete or continuous) of the variables of an object.
#' 
#' Get a vector representing the status of the variables (with their names) of a \code{\link{BN}} or \code{\link{BNDataset}}.
#' Elements of the vector are \code{c} if the variable is continue, and \code{d} if the variable is discrete.
#' 
#' @name discreteness
#' @rdname discreteness
#'
#' @param x an object.
#' 
#' @return vector contaning, for each variable of the desired object,
#'         \code{c} if the variable is continue, and \code{d} if the variable is discrete.
#' 
#' @exportMethod discreteness
setGeneric("discreteness", function(x) standardGeneric("discreteness"))


#' get size of the variables of an object.
#' 
#' Return a list containing the size of the variables of an object. It is the actual cardinality
#' of discrete variables, and the cardinality of the discretized variable for continuous variables.
#' 
#' @name node.sizes
#' @rdname node.sizes
#' 
#' @param x an object.
#' 
#' @return vector contaning the size of each variable of the desired object.
#' 
#' @exportMethod node.sizes
setGeneric("node.sizes", function(x) standardGeneric("node.sizes"))


#' get the list of conditional probability tables of a \code{\link{BN}}.
#' 
#' Return the list of conditional probability tables of the variables of a \code{\link{BN}} object.
#' Each probability table is associated to the corresponding variable, and its dimensions are named according
#' to the variable they represent.
#' 
#' Each conditional probability table is represented as a multidimensional array. 
#' The ordering of the dimensions of each variable is not guaranteed to follow the actual conditional distribution.
#' E.g. dimensions for conditional probability \code{P(C|A,B)} can be either \code{(C,A,B)} or \code{(A,B,C)}, depending on
#' if some operations have been performed, or how the probability table has been computed.
#' Users should not rely on dimension numbers, but should instead select the dimensions using their names.
#' 
#' @name cpts
#' @rdname cpts
#' 
#' @param x an object.
#' 
#' @return list of the conditional probability tables of the desired object.
#' 
#' @exportMethod cpts
setGeneric("cpts", function(x) standardGeneric("cpts"))


#' get adjacency matrix of a network.
#' 
#' Return the adjacency matrix of the directed acyclic graph representing the structure of a network.
#' 
#' @name dag
#' @rdname dag
#'
#' @param x an object.
#' 
#' @return matrix containing the adjacency matrix of the directed acyclic graph representing
#'         the structure of the object.
#' 
#' @exportMethod dag
setGeneric("dag", function(x) standardGeneric("dag"))


#' get the WPDAG of an object.
#' 
#' Return the weighted partially directed acyclic graph of a network, when available (e.g. when bootstrap on dataset is performed).
#' 
#' @name wpdag
#' @rdname wpdag
#'
#' @param x an object.
#' 
#' @return matrix contaning the WPDAG of the object.
#' 
#' @exportMethod wpdag
setGeneric("wpdag", function(x) standardGeneric("wpdag"))


#' get header file of a \code{\link{BNDataset}}.
#' 
#' Return the header filename of a dataset (with the path to its position, as given by the user),
#' present if the dataset has been read from a file and not manually inserted.
#' The header file contains three rows:
#' \enumerate{
#' \item list of names of the variables, in the same order as in the data file;
#' \item list of cardinalities of the variables, if discrete, or levels for quantization if continuous;
#' \item list of status of the variables: \code{c} for continuous variables, \code{d} for discrete ones.
#' }
#' 
#' @name header.file
#' @rdname header.file
#' 
#' @param x a \code{\link{BNDataset}}.
#' 
#' @return header filename of the dataset.
#' 
#' @seealso \code{\link{data.file}}
#' 
#' @exportMethod header.file
setGeneric("header.file", function(x) standardGeneric("header.file"))


#' get data file of a \code{\link{BNDataset}}.
#' 
#' Return the data filename of a dataset (with the path to its position, as given by the user).
#' The data filename may contain a header in the first row, containing the list of names of the variables,
#' in the same order as in the header file.
#' After the header, if present, the file contains a data.frame with the observations, one item per row.
#' 
#' @name data.file
#' @rdname data.file
#'
#' @param x a \code{\link{BNDataset}}.
#' 
#' @return data filename of the dataset.
#' 
#' @seealso \code{\link{data.file}}
#' 
#' @exportMethod data.file
setGeneric("data.file", function(x) standardGeneric("data.file"))


#' get number of variables of a \code{\link{BNDataset}}.
#' 
#' Return the number of the variables contained in a dataset. This value corresponds to the value
#' of \code{\link{num.nodes}} of a network built upon the same dataset.
#' 
#' @name num.variables
#' @rdname num.variables
#' 
#' @param x a \code{\link{BNDataset}} object.
#' 
#' @return number of variables of the desired dataset.
#' 
#' @seealso \code{\link{num.nodes}}
#' 
#' @exportMethod num.variables
setGeneric("num.variables", function(x) standardGeneric("num.variables"))


#' get number of items of a \code{\link{BNDataset}}.
#' 
#' Return the number of items in a dataset, that is, the number of rows in its data slot.
#' 
#' @name num.items
#' @rdname num.items
#' 
#' @param x a \code{\link{BNDataset}} object.
#' 
#' @return number of items of the desired dataset.
#' 
#' @exportMethod num.items
setGeneric("num.items", function(x) standardGeneric("num.items"))


#' check whether a \code{\link{BNDataset}} has bootstrap samples or not.
#' 
#' Return \code{TRUE} if the given dataset contains samples for bootstrap, \code{FALSE} otherwise.
#' 
#' @name has.boots
#' @rdname has.boots
#'
#' @param x a \code{\link{BNDataset}} object.
#' 
#' @return \code{TRUE} if dataset has bootstrap samples.
#' 
#' @seealso \code{\link{has.imputed.boots}}, \code{\link{boots}}, \code{\link{imp.boots}}
#' 
#' @exportMethod has.boots
setGeneric("has.boots", function(x) standardGeneric("has.boots"))


#' check whether a \code{\link{BNDataset}} has bootstrap samples from imputed data or not.
#' 
#' Return \code{TRUE} if the given dataset contains samples for bootstrap from inputed dataset, \code{FALSE} otherwise.
#' 
#' @name has.imputed.boots
#' @name has.imputed.boots
#'
#' @param x a \code{\link{BNDataset}} object.
#' 
#' @return \code{TRUE} if dataset has bootstrap samples from imputed data.
#' 
#' @seealso \code{\link{has.boots}}, \code{\link{boots}}, \code{\link{imp.boots}}
#' 
#' @exportMethod has.imputed.boots
setGeneric("has.imputed.boots", function(x) standardGeneric("has.imputed.boots"))


#' get list of bootstrap samples of a \code{\link{BNDataset}}.
#' 
#' Return the list of samples computed from raw data of a dataset.
#' 
#' @name boots
#' @rdname boots
#'
#' @param x a \code{\link{BNDataset}} object.
#' 
#' @return the list of bootstrap samples.
#' 
#' @seealso \code{\link{has.boots}}, \code{\link{has.imputed.boots}}, \code{\link{imp.boots}}
#' 
#' @exportMethod boots
setGeneric("boots", function(x) standardGeneric("boots"))


#' get list of bootstrap samples from imputed data of a \code{\link{BNDataset}}.
#' 
#' Return the list of samples computed from raw data of a dataset.
#' 
#' @name imp.boots
#' @rdname imp.boots
#'
#' @param x a \code{\link{BNDataset}} object.
#' 
#' @return the list of bootstrap samples from imputed data.
#' 
#' @seealso \code{\link{has.boots}}, \code{\link{has.imputed.boots}}, \code{\link{boots}}
#' 
#' @exportMethod imp.boots
setGeneric("imp.boots", function(x) standardGeneric("imp.boots"))


#' get number of bootstrap samples of a \code{\link{BNDataset}}.
#' 
#' Return the number of bootstrap samples computed from a dataset.
#' 
#' @name num.boots
#' @rdname num.boots
#' 
#' @param x a \code{\link{BNDataset}} object.
#' 
#' @return the number of bootstrap samples.
#' 
#' @exportMethod num.boots
setGeneric("num.boots", function(x) standardGeneric("num.boots"))


#' get the junction tree of an \code{\link{InferenceEngine}}.
#' 
#' Return the adjacency matrix representing the junction tree computed for a network.
#' 
#' Rows and columns are named after the (variables in the) cliques that each node of the junction tree represent.
#' 
#' @name junction.tree
#' @rdname junction.tree
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' 
#' @return the junction tree contained in the \code{\link{InferenceEngine}}.
#' 
#' @seealso \code{\link{build.junction.tree}}
#' 
#' @exportMethod junction.tree
setGeneric("junction.tree", function(x) standardGeneric("junction.tree"))


#' get the list of cliques of the junction tree of an \code{\link{InferenceEngine}}.
#' 
#' Return the list of cliques containing the variables associated to each node of a junction tree.
#' 
#' @name jt.cliques
#' @rdname jt.cliques
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' 
#' @return the list of cliques of the junction tree contained in the \code{\link{InferenceEngine}}.
#' 
#' @exportMethod jt.cliques
setGeneric("jt.cliques", function(x) standardGeneric("jt.cliques"))


#' get the list of joint probability tables compiled by an \code{\link{InferenceEngine}}.
#' 
#' Return the list of joint probability tables for the cliques of the junction tree 
#' obtained after belief propagation has been performed.
#' 
#' Each joint probability table is represented as a multidimensional array. 
#' To retrieve single dimensions (e.g. to compute marginals), users should not rely on dimension numbers,
#' but should instead select the dimensions using their names.
#' 
#' @name jpts
#' @rdname jpts
#'
#' @param x an \code{\link{InferenceEngine}}.
#' 
#' @return the list of joint probability tables compiled by the \code{\link{InferenceEngine}}.
#' 
#' @exportMethod jpts
setGeneric("jpts", function(x) standardGeneric("jpts"))


#' get the \code{\link{BN}} object contained in an \code{\link{InferenceEngine}}.
#' 
#' Return a network contained in an InferenceEngine.
#' 
#' @name bn
#' @rdname bn-method
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' 
#' @return the \code{\link{BN}} object contained in an \code{\link{InferenceEngine}}.
#'         
#' @exportMethod bn
setGeneric("bn", function(x) standardGeneric("bn"))


#' get the updated \code{\link{BN}} object contained in an \code{\link{InferenceEngine}}.
#' 
#' Return an updated network contained in an InferenceEngine.
#' 
#' @name updated.bn
#' @rdname updated.bn-method
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' 
#' @return the updated \code{\link{BN}} object contained in an \code{\link{InferenceEngine}}.
#'         
#' @exportMethod updated.bn
setGeneric("updated.bn", function(x) standardGeneric("updated.bn"))

#' get the list of observations of an \code{\link{InferenceEngine}}.
#' 
#' Return the list of observations added to an InferenceEngine.
#' 
#' Output is a list in the following format:
#' \itemize{
#' \item{\code{observed.vars}}{vector of observed variables;}
#' \item{\code{observed.vals}}{vector of values observed for the variables in \code{observed.vars} in the corresponding position.}
#' }
#' 
#' @name observations
#' @rdname observations
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' 
#' @return the list of observations of the \code{\link{InferenceEngine}}.
#' 
#' @exportMethod observations
setGeneric("observations", function(x) standardGeneric("observations"))


###############################################################################


#' set name of an object.
#' 
#' Set the \code{name} slot of an object of type \code{\link{BN}} or \code{\link{BNDataset}}.
#' 
#' @name name<-
#' @rdname name-set
#' 
#' @param x an object.
#' @param value the new name of the object.
#' 
#' @exportMethod name<-
setGeneric("name<-", function(x, value) standardGeneric("name<-"))


#' set number of nodes of an object.
#' 
#' Set the number of nodes of an object of type \code{\link{BN}} (number of nodes of the network)
#' or \code{\link{InferenceEngine}} (where parameter contains the number of nodes of the junction tree).
#' 
#' @name num.nodes<-
#' @rdname num.nodes-set
#' 
#' @param x an object.
#' @param value the number of nodes in the object.
#' 
#' @exportMethod num.nodes<-
setGeneric("num.nodes<-", function(x, value) standardGeneric("num.nodes<-"))


#' set variables of an object.
#' 
#' Set the list of variable names in a \code{\link{BN}} or \code{\link{BNDataset}} object.
#' 
#' @name variables<-
#' @rdname variables-set
#' 
#' @param x an object.
#' @param value vector containing the variable names of the object.
#'        Overwrites \code{num.nodes} slot if non-matching.
#' 
#' @exportMethod variables<-
setGeneric("variables<-", function(x, value) standardGeneric("variables<-"))


#' set status (discrete or continuous) of the variables of an object.
#' 
#' Set the list of variable status for the variables in a network or a dataset.
#' 
#' @name discreteness<-
#' @rdname discreteness-set
#' 
#' @param x an object.
#' @param value a vector of elements in \{\code{c},\code{d}\} for continuous and discrete variables (respectively).
#' 
#' @exportMethod discreteness<-
setGeneric("discreteness<-", function(x, value) standardGeneric("discreteness<-"))


#' set the size of variables of an object.
#' 
#' Set the size of the variables of a BN or BNDataset object. It represents the actual cardinality
#' of discrete variables, and the cardinality of the discretized variable for continuous variables.
#' 
#' @name node.sizes<-
#' @rdname node.sizes-set
#' 
#' @param x an object.
#' @param value vector contaning the size of each variable of the object.
#' 
#' @exportMethod node.sizes<-
setGeneric("node.sizes<-", function(x, value) standardGeneric("node.sizes<-"))


#' set the list of conditional probability tables of a network.
#' 
#' Set the list of conditional probability tables of a \code{\link{BN}} object.
#' 
#' Each conditional probability table is represented as a multidimensional array. 
#' To retrieve single dimensions (e.g. to compute marginals), users should provide dimensions names.
#' 
#' @name cpts<-
#' @rdname cpts-set
#' 
#' @param x an object.
#' @param value list of the conditional probability tables of the object.
#' 
#' @exportMethod cpts<-
setGeneric("cpts<-", function(x, value) standardGeneric("cpts<-"))


#' set adjacency matrix of an object.
#' 
#' Set the adjacency matrix of the directed acyclic graph representing the structure of a network.
#' 
#' @name dag<-
#' @rdname dag-set
#' 
#' @param x an object.
#' @param value matrix containing the adjacency matrix of the directed acyclic graph representing
#'        the structure of the object.
#'         
#' @exportMethod dag<-
setGeneric("dag<-", function(x, value) standardGeneric("dag<-"))


#' set WPDAG of the object.
#' 
#' Set the weighted partially directed acyclic graph of a network (e.g. in case bootstrap on dataset is performed).
#' 
#' @name wpdag<-
#' @rdname wpdag-set
#' 
#' @param x an object.
#' @param value matrix contaning the WPDAG of the object.
#' 
#' @exportMethod wpdag<-
setGeneric("wpdag<-", function(x, value) standardGeneric("wpdag<-"))


#' set header file of a \code{\link{BNDataset}}.
#' 
#' Set the header filename of a dataset (with the path to its position, as given by the user).
#' The header file has to contain three rows:
#' \enumerate{
#' \item list of names of the variables, in the same order as in the data file;
#' \item list of cardinalities of the variables, if discrete, or levels for quantization if continuous;
#' \item list of status of the variables: \code{c} for continuous variables, \code{d} for discrete ones.
#' }
#' Further rows are ignored.
#' 
#' @name header.file<-
#' @rdname header.file-set
#' 
#' @param x a \code{\link{BNDataset}}.
#' @param value header filename.
#' 
#' @seealso \code{\link{data.file<-}}
#' 
#' @exportMethod header.file<-
setGeneric("header.file<-", function(x, value) standardGeneric("header.file<-"))


#' set data file of a \code{\link{BNDataset}}.
#' 
#' Set the data filename of a dataset (with the path to its position, as given by the user).
#' The data filename may contain a header in the first row, containing the list of names of the variables,
#' in the same order as in the header file.
#' After the header, if present, the file contains a data.frame with the observations, one item per row.
#' 
#' @name data.file<-
#' @rdname data.file-set
#' 
#' @param x a \code{\link{BNDataset}}.
#' @param value data filename.
#' 
#' @seealso \code{\link{header.file<-}}
#' 
#' @exportMethod data.file<-
setGeneric("data.file<-", function(x, value) standardGeneric("data.file<-"))


#' set number of variables of a \code{\link{BNDataset}}.
#' 
#' Set the number of variables observed in a dataset.
#' 
#' @name num.variables<-
#' @rdname num.variables-set
#' 
#' @param x a \code{\link{BNDataset}} object.
#' @param value number of variables of the dataset.
#' 
#' @exportMethod num.variables<-
setGeneric("num.variables<-", function(x, value) standardGeneric("num.variables<-"))


#' set number of items of a \code{\link{BNDataset}}.
#' 
#' Set the number of observed items (rows) in a dataset.
#' 
#' @name num.items<-
#' @rdname num.items-set
#' 
#' @param x a \code{\link{BNDataset}} object.
#' @param value number of items of the desired dataset.
#' 
#' @exportMethod num.items<-
setGeneric("num.items<-", function(x, value) standardGeneric("num.items<-"))


#' set list of bootstrap samples of a \code{\link{BNDataset}}.
#' 
#' Add to a dataset a list of samples from raw data computed using bootstrap.
#' 
#' @name boots<-
#' @rdname boots-set
#' 
#' @param x a \code{\link{BNDataset}} object.
#' @param value the list of bootstrap samples.
#' 
#' @exportMethod boots<-
setGeneric("boots<-", function(x, value) standardGeneric("boots<-"))


#' set number of bootstrap samples of a \code{\link{BNDataset}}.
#' 
#' Set the length of the list of samples of a dataset computed using bootstrap.
#' 
#' @name num.boots<-
#' @rdname num.boots-set
#' 
#' @param x a \code{\link{BNDataset}} object.
#' @param value the number of bootstrap samples.
#' 
#' @exportMethod num.boots<-
setGeneric("num.boots<-", function(x, value) standardGeneric("num.boots<-"))


#' set list of bootstrap samples from imputed data of a \code{\link{BNDataset}}.
#' 
#' Add to a dataset a list of samples from imputed data computed using bootstrap.
#' 
#' @name imp.boots<-
#' @rdname imp.boots-set
#' 
#' @param x a \code{\link{BNDataset}} object.
#' @param value the list of bootstrap samples from imputed data.
#' 
#' @exportMethod imp.boots<-
setGeneric("imp.boots<-", function(x, value) standardGeneric("imp.boots<-"))


#' set the junction tree of an \code{\link{InferenceEngine}}.
#' 
#' Set the adjacency matrix of the junction tree computed for a network.
#' 
#' @name junction.tree<-
#' @rdname junction.tree-set
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' @param value the junction tree to be inserted in the \code{\link{InferenceEngine}}.
#' 
#' @exportMethod junction.tree<-
setGeneric("junction.tree<-", function(x, value) standardGeneric("junction.tree<-"))


#' set the list of cliques of the junction tree of an \code{\link{InferenceEngine}}.
#' 
#' Add to the InferenceEngine a list containing the cliques of variables composing the nodes of the junction tree.
#' 
#' @name jt.cliques<-
#' @rdname jt.cliques-set
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' @param value the list of cliques of the junction tree contained in the \code{\link{InferenceEngine}}.
#'
#' @exportMethod jt.cliques<-
setGeneric("jt.cliques<-", function(x, value) standardGeneric("jt.cliques<-"))


#' set the list of joint probability tables compiled by an \code{\link{InferenceEngine}}.
#' 
#' Add a list of joint probability tables for the cliques of the junction tree.
#' 
#' Each joint probability table is represented as a multidimensional array. 
#' To retrieve single dimensions (e.g. to compute marginals), users should provide dimension names.
#' 
#' @name jpts<-
#' @rdname jpts-set
#'
#' @param x an \code{\link{InferenceEngine}}.
#' @param value the list of joint probability tables compiled by the \code{\link{InferenceEngine}}.
#' 
#' @exportMethod jpts<-
setGeneric("jpts<-", function(x, value) standardGeneric("jpts<-"))


#' set the original \code{\link{BN}} object contained in an \code{\link{InferenceEngine}}.
#' 
#' Add an original network to an InferenceEngine.
#' 
#' @name bn<-
#' @rdname bn-set
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' @param value the \code{\link{BN}} object contained in an \code{\link{InferenceEngine}}.
#'         
#' @exportMethod bn<-
setGeneric("bn<-", function(x, value) standardGeneric("bn<-"))


#' set the updated \code{\link{BN}} object contained in an \code{\link{InferenceEngine}}.
#' 
#' Add an updated network to an InferenceEngine.
#' 
#' @name updated.bn<-
#' @rdname updated.bn-set
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' @param value the updated \code{\link{BN}} object contained in an \code{\link{InferenceEngine}}.
#'         
#' @exportMethod updated.bn<-
setGeneric("updated.bn<-", function(x, value) standardGeneric("updated.bn<-"))



#' set the list of observations of an \code{\link{InferenceEngine}}.
#' 
#' Add a list of observations to an InferenceEngine, using a list of observations composed by the two following vectors:
#' \itemize{
#' \item{\code{observed.vars}}{vector of observed variables;}
#' \item{\code{observed.vals}}{vector of values observed for the variables in \code{observed.vars} in the corresponding position.}
#' }
#' 
#' Replace previous list of observations, if present. In order to add evidence, and not just replace it,
#' one must use the \code{\link{add.observations<-}} method.
#' 
#' In case of multiple observations of the same variable, the last observation is the one used, as the most recent.
#' 
#' @name observations<-
#' @rdname observations-set
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' @param value the list of observations of the \code{\link{InferenceEngine}}.
#' 
#' @seealso \code{\link{add.observations<-}}
#' 
#' @exportMethod observations<-
setGeneric("observations<-", function(x, value) standardGeneric("observations<-"))


#' add further evidence to an existing list of observations of an \code{\link{InferenceEngine}}.
#' 
#' Add a list of observations to an InferenceEngine that already has observations,
#' using a list composed by the two following vectors:
#' \itemize{
#' \item{\code{observed.vars}}{vector of observed variables;}
#' \item{\code{observed.vals}}{vector of values observed for the variables in \code{observed.vars} in the corresponding position.}
#' }
#' 
#' In case of multiple observations of the same variable, the last observation is the one used, as the most recent.
#' 
#' @name add.observations<-
#' @rdname add.observations-set
#' 
#' @param x an \code{\link{InferenceEngine}}.
#' @param value the list of observations of the \code{\link{InferenceEngine}}.
#' 
#' @seealso \code{\link{observations<-}}
#' 
#' @exportMethod add.observations<-
setGeneric("add.observations<-", function(x, value) standardGeneric("add.observations<-"))


###############################################################################
###############################################################################

# other common generics

###############################################################################


#' print an object to \code{stdout}.
#' 
#' @name print
#' 
#' @param x an object.
#' @param show.raw.data when \code{x} is a \code{\link{BNDataset}}, print also raw dataset, if available.
#' @param show.imputed.data when \code{x} is a \code{\link{BNDataset}}, print also imputed dataset, if available.
#' @param engine when \code{x} is an \code{\link{InferenceEngine}}, specify the inference engine to be shown.
#'        Currently only \code{engine = 'jt'} is supported.
#' @param ... potential other arguments.
#' 
#' @exportMethod print
setGeneric("print", function(x, ...) standardGeneric("print"))
