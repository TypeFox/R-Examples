setClassUnion("missingOrNumeric", c("missing","numeric"))
setClassUnion("missingOrLogical", c("missing","logical"))
setClassUnion("missingOrInteger", c("missing","integer"))
setClassUnion("missingOrCharacter", c("missing","character"))
setClassUnion("missingOrNULL", c("missing", "NULL"))
setClassUnion("vectorOrNULL", c("vector", "NULL"))

###############################################################################
#
# Generic Bayesian Network class
#
###############################################################################

#' BN class definition.
#' 
#' @section Slots:
#' \describe{
#'   \item{\code{name}:}{name of the network}
#'   \item{\code{num.nodes}:}{number of nodes in the network}
#'   \item{\code{variables}:}{names of the variables in the network}
#'   \item{\code{discreteness}:}{\code{TRUE} if variable is discrete, \code{FALSE} if variable is continue}
#'   \item{\code{node.sizes}:}{if variable \code{i} is discrete, \code{node.sizes[i]} contains the cardinality of \code{i},
#'      if \code{i} is instead discrete the value is the number of states variable \code{i} takes when discretized}
#'   \item{\code{cpts}:}{list of conditional probability tables of the network}
#'   \item{\code{dag}:}{adjacency matrix of the network}
#'   \item{\code{wpdag}:}{weighted partially dag}
#'   \item{\code{scoring.func}:}{scoring function used in structure learning (when performed)}
#'   \item{\code{struct.algo}:}{algorithm used in structure learning (when performed)}
#' }
#' 
#' @name BN-class
#' @rdname BN-class
#' @docType class
#' @aliases BN,BN-class
#' 
#' @exportClass BN
setClass("BN",
         representation(
           name         = "character",
           num.nodes    = "numeric",
           variables    = "character",
           discreteness = "logical",
           node.sizes   = "numeric",
           cpts         = "list",
           dag          = "matrix",
           wpdag        = "matrix",
           scoring.func = "character",
           struct.algo  = "character"
         ),
         prototype(
           name         = "",
           num.nodes    = 0,
           variables    = c(""),
           discreteness = c(TRUE),
           node.sizes   = c(0),
           cpts         = list(NULL),
           dag          = matrix(),
           wpdag        = matrix(),
           scoring.func = "",
           struct.algo  = ""
         )
        )


# Create new class union to allow the embedding of a BN in a slot of another class,
# allowing it to take the default value of NULL.
setClassUnion("BNOrNULL", members=c("BN", "NULL"))


###############################################################################
#
# Dataset class
#
###############################################################################

#' BNDataset class.
#' 
#' Contains the all of the data that can be extracted from a given dataset:
#' raw data, imputed data, raw and imputed data with bootstrap.
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
#' In case of need of more advanced options when reading a dataset from files, please refer to the
#' documentation of the \code{\link{read.dataset}} method. Imputation and bootstrap are also available
#' as separate routines (\code{\link{impute}} and \code{\link{bootstrap}}, respectively).
#'
#'  
#' @param data raw data.frame or path/name of the file containing the raw dataset (see 'Details').
#' @param discreteness a vector of booleans indicating if the variables are discrete or continuous
#'   (\code{TRUE} and \code{FALSE}, respectively),
#'   or path/name of the file containing header information for the dataset (discreteness, variable names, cardinality - see 'Details').
#' @param variables vector of variable names.
#' @param node.sizes vector of variable cardinalities (for discrete variables) or quantization ranges (for continuous variables).
#' @param ... further arguments for reading a dataset from files (see documentation for \code{read.dataset}).
#' 
#' @return BNDataset object.
#' 
#' @seealso read.dataset, impute, bootstrap
#' 
#' @section Slots:
#' \describe{
#'   \item{\code{name}:}{name of the dataset}
#'   \item{\code{header.file}:}{name and location of the header file}
#'   \item{\code{data.file}:}{name and location of the data file}
#'   \item{\code{variables}:}{names of the variables in the network}
#'   \item{\code{node.sizes}:}{cardinality of each variable of the network}
#'   \item{\code{num.variables}:}{number of variables (columns) in the dataset}
#'   \item{\code{discreteness}:}{\code{TRUE} if variable is discrete, \code{FALSE} if variable is continue}
#'   \item{\code{num.items}:}{number of observations (rows) in the dataset}
#'   \item{\code{has.raw.data}:}{TRUE if the dataset contains data read from a file}
#'   \item{\code{has.imputed.data}:}{TRUE if the dataset contains imputed data (computed from raw data)}
#'   \item{\code{raw.data}:}{matrix containing raw data}
#'   \item{\code{imputed.data}:}{matrix containing imputed data}
#'   \item{\code{has.boots}:}{dataset has bootstrap samples}
#'   \item{\code{boots}:}{list of bootstrap samples}
#'   \item{\code{has.imputed.boots}:}{dataset has imputed bootstrap samples}
#'   \item{\code{imp.boots}:}{list of imputed bootstrap samples}
#'   \item{\code{num.boots}:}{number of bootstrap samples}
#' }
#' 
#' @name BNDataset-class
#' @rdname BNDataset-class
#' @docType class
#' @aliases BNDataset,BNDataset-class
#' 
#' @usage BNDataset(data, discreteness, variables = NULL, node.sizes = NULL, ...)
#' 
#' @examples
#' \dontrun{
#' # create from files
#' dataset <- BNDataset("file.data", "file.header")
#' 
#' # other way: create from raw dataset and metadata
#' data <- matrix(c(1:16), nrow = 4, ncol = 4)
#' dataset <- BNDataset(data = data,
#'                      discreteness = rep('d',4),
#'                      variables = c("a", "b", "c", "d"),
#'                      node.sizes = c(4,8,12,16))
#' }
#'
#' @exportClass BNDataset
setClass("BNDataset",
         representation(
           name              = "character",
           header.file       = "character",
           data.file         = "character",
           variables         = "character",
           node.sizes        = "numeric",
           num.variables     = "numeric",
           discreteness      = "logical",
           num.items         = "numeric",
           has.raw.data      = "logical",
           has.imputed.data  = "logical",
           raw.data          = "matrix",
           imputed.data      = "matrix",
           has.boots         = "logical",
           boots             = "list",
           has.imputed.boots = "logical",
           imp.boots         = "list",
           num.boots         = "numeric"
         ),
         prototype(
           name              = "",
           header.file       = "",
           data.file         = "",
           variables         = c(""),
           node.sizes        = c(0),
           num.variables     = 0,
           discreteness      = c(TRUE),
           num.items         = 0,
           has.raw.data      = FALSE,
           has.imputed.data  = FALSE,
           raw.data          = matrix(c(0)),
           imputed.data      = matrix(c(0)),
           has.boots         = FALSE,
           boots             = list(NULL),
           has.imputed.boots = FALSE,
           imp.boots         = list(NULL),
           num.boots         = 0
         )
        )


###############################################################################
#
# InferenceEngine class
#
###############################################################################


#' InferenceEngine class.
#' 
#' @section Slots:
#' \describe{
#'   \item{\code{junction.tree}:}{junction tree adjacency matrix.}
#'   \item{\code{num.nodes}:}{number of nodes in the junction tree.}
#'   \item{\code{cliques}:}{list of cliques composing the nodes of the junction tree.}
#'   \item{\code{triangulated.graph}:}{adjacency matrix of the original triangulated graph.}
#'   \item{\code{jpts}:}{inferred joint probability tables.}
#'   \item{\code{bn}:}{original Bayesian Network (as object of class \code{\link{BN}}) as provided by the user, or learnt from a dataset.
#'          \code{NULL} if missing.}
#'   \item{\code{updated.bn}:}{Bayesian Network  (as object of class \code{\link{BN}}) as modified by a belief propagation computation. In particular,
#'          it will have different conditional probability tables with respect to its original version. \code{NULL} if missing.}
#'   \item{\code{observed.vars}:}{list of observed variables, by name or number.}
#'   \item{\code{observed.vals}:}{list of observed values for the corresponding variables in \code{observed.vars}.}
#' }
#' 
#' 
#' @name InferenceEngine-class
#' @docType class
#' @rdname InferenceEngine-class
#' @aliases InferenceEngine,InferenceEngine-class
#' 
#' @exportClass InferenceEngine
setClass("InferenceEngine",
         representation(
           junction.tree      = "matrix",
           num.nodes          = "numeric",
           cliques            = "list",
           triangulated.graph = "matrix",
           jpts               = "list",
           bn                 = "BNOrNULL",
           updated.bn         = "BNOrNULL",
           observed.vars      = "vectorOrNULL",
           observed.vals      = "vectorOrNULL"
         ),
         prototype(
           junction.tree      = matrix(),
           num.nodes          = 0,
           cliques            = list(NULL),
           triangulated.graph = matrix(),
           jpts               = list(NULL),
           bn                 = NULL,
           updated.bn         = NULL,
           observed.vars      = NULL,
           observed.vals      = NULL
         )
        )


####
setClassUnion("AllTheClasses", c("BN", "BNDataset", "InferenceEngine"))