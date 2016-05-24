#' @useDynLib sdcTable
#' @import "methods"
#' @import "Rcpp"
#' @import "Rglpk"
#' @import "stringr"
#' @import "lpSolveAPI"
#' @import "data.table"
#' @importFrom "stats" "na.omit"
#' @importFrom "utils" "combn"
#' @importFrom "utils" "tail"
#' @importFrom "utils" "read.table"
#' @importFrom "utils" "str"

setClassUnion('dataframeOrNULL', c('data.frame', 'NULL'))
setClassUnion('numericOrNULL', c('numeric', 'NULL'))
setClassUnion('characterOrNULL', c('character', 'NULL'))
setClassUnion('logicalOrNULL', c('logical', 'NULL'))
setClassUnion('matrixOrNULL', c('matrix', 'NULL'))
setClassUnion('listOrNULL', c('list', 'NULL'))

#' S4 class describing a dataObj-object
#'
#' This class models a data object containing the 'raw' data for a given problem
#' as well as information on the position of the dimensional variables, the count
#' variable, additional numerical variables, weights or sampling weights within the
#' raw data. Also slot 'isMicroData' shows if slow 'rawData' consists of microdata
#' (multiple observations for each cell are possible, isMicroData==TRUE) or if data
#' have already been aggregated (isMicroData==FALSE)
#'
#' \describe{
#' \item{slot \code{rawData}:}{list with each element being a vector of either codes of dimensional variables, counts, weights that should be used for secondary cell suppression problem, numerical variables or sampling weights. }
#' \item{slot \code{dimVarInd}:}{numeric vector (or NULL) defining the indices of the dimensional variables within slot 'rawData'}
#' \item{slot \code{freqVarInd}:}{numeric vector (or NULL) defining the indices of the frequency variables within slot 'rawData'}
#' \item{slot \code{numVarInd}:}{numeric vector (or NULL) defining the indices of the numerical variables within slot 'rawData'}
#' \item{slot \code{weightVarInd}:}{numeric vector (or NULL) defining the indices of the variables holding weights within slot 'rawData'}
#' \item{slot \code{sampWeightInd}:}{numeric vector (or NULL) defining the indices of the variables holding sampling weights within slot 'rawData'}
#' \item{slot \code{isMicroData}:}{logical vector of length 1 (or NULL) that is TRUE if slot 'rawData' are microData and FALSE otherwise }
#' }
#' @name dataObj-class
#' @rdname dataObj-class
#' @exportClass dataObj
#' @note objects of class \code{dataObj} are input for slot \code{dataObj} in class \code{sdcProblem}
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setClass(
  Class='dataObj',
  representation=representation(
    rawData='listOrNULL',
    dimVarInd='numericOrNULL',
    freqVarInd='numericOrNULL',
    numVarInd='numericOrNULL',
    weightVarInd='numericOrNULL',
    sampWeightInd='numericOrNULL',
    isMicroData='logicalOrNULL'
  ),
  prototype=prototype(
    rawData=NULL,
    dimVarInd=NULL,
    freqVarInd=NULL,
    numVarInd=NULL,
    weightVarInd=NULL,
    sampWeightInd=NULL,
    isMicroData=NULL
  ),
  validity=function(object) {
    if ( !all(g_freqvar_ind(object) %in% 1:length(g_raw_data(object))) ) {
      stop("dataObj:: check input parameter 'freqVarInd'!\n")
    }
    if ( !all(g_numvar_ind(object) %in% 1:length(g_raw_data(object))) ) {
      stop("dataObj:: check input parameter 'numVarInd'!\n")
    }
    if ( length(g_weightvar_ind(object)) > 1 ) {
      stop("dataObj:: length of parameter 'weightVarInd' must not be greater than 1!\n")
    }
    if ( length(g_sampweight_ind(object)) > 1 ) {
      stop("dataObj:: length of parameter 'sampWeightInd' must not be greater than 1!\n")
    }
    return(TRUE)
  }
)

#' S4 class describing a dimInfo-object
#'
#' An object of class \code{dimInfo} holds all necessary information about the
#' dimensional variables defining a hierarchical table that needs to be protected.
#'
#' \describe{
#' \item{slot \code{dimInfo}:}{a list (or NULL) with all list elements being objects of class \code{dimVar}}
#' \item{slot \code{strID}:}{a character vector (or NULL) defining IDs that identify each table cell. The ID's are based on (default) codes of the dimensional variables defining a cell. }
#' \item{slot \code{strInfo}:}{a list object (or NULL) with each list element being a numeric vector of length 2 defining the start and end-digit that is allocated by the i-th dimensional variable in ID-codes available in slot \code{strID}  }
#' \item{slot \code{vNames}:}{a character vector (or NULL) defining the variable names of the dimensional variables defining the table structure}
#' \item{slot \code{posIndex}:}{a numeric vector (or NULL) holding the position of the dimensional variables within slot \code{rawData} of class \code{dataObj} }
#' }
#' @name dimInfo-class
#' @rdname dimInfo-class
#' @exportClass dimInfo
#' @note objects of class \code{dimInfo} are input for slots in classes \code{sdcProblem} and \code{safeObj}
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setClass(
  Class='dimInfo',
    representation=representation(
    dimInfo='listOrNULL',
    strID='characterOrNULL',
    strInfo='listOrNULL',
    vNames='characterOrNULL',
    posIndex='numericOrNULL'
  ),
  prototype=prototype(
    dimInfo=NULL,
    strID=NULL,
    strInfo=NULL,
    vNames=NULL,
    posIndex=NULL
  ),
  validity=function(object) {
    if ( length(g_varname(object)) != length(g_pos_index(object)) ) {
      stop("dimInfo:: parameter 'vNames' and 'posIndex' differ in length!\n")
    }
    if ( length(g_varname(object)) != length(g_dim_info(object)) ) {
      stop("dimInfo:: parameter 'vNames' and 'dimInfo' differ in length!\n")
    }
    if ( length(g_str_info(object)) != length(g_varname(object)) ) {
      stop("dimInfo:: parameter 'strInfo' and 'vNames' differ in length!\n")
    }
    if ( any(sapply(g_dim_info(object), "class") != "dimVar" ) ) {
      stop("dimInfo:: elements of parameter 'dimInfo' must be of class 'dimVar'!\n")
    }
    return(TRUE)
  }
)

#' S4 class describing a dimVar-object
#'
#' An object of class \code{dimVar} holds all necessary information about a single
#' dimensional variable such as original and standardized codes, the level-structure,
#' the hierarchical structure, codes that may be (temporarily) removed from
#' building the complete hierarchy (dups) and their corresponding codes that correspond
#' to these duplicated codes.
#'
#' \describe{
#' \item{slot \code{codesOriginal}:}{a character vector (or NULL) holding original variable codes}
#' \item{slot \code{codesDefault}:}{a character vector (or NULL) holding standardized codes}
#' \item{slot \code{codesMinimal}:}{a logical vector (or NULL) defining if a code is required to build the complete hierarchy or not (then the code is a (sub)total)}
#' \item{slot \code{vName}:}{character vector of length 1 (or NULL) defining the variable name of the dimensional variable}
#' \item{slot \code{levels}:}{a numeric vector (or NULL) defining the level structure. For each code the corresponding level is listed with the grand-total always having level==1}
#' \item{slot \code{structure}:}{a numeric vector (or NULL) with length of the total number of levels. Each element shows how many digits the i-th level allocates within the standardized codes (note: level 1 always allocates exactly 1 digit in the standardized codes)}
#' \item{slot \code{dims}:}{a list (or NULL) defining the hierarchical structure of the dimensional variable. Each list-element is a character vector with elements available in slot \code{codesDefault} and the first element always being a (sub)total and the remaining elements being the codes that contribute to the (sub)total}
#' \item{slot \code{dups}:}{character vector (or NULL) having showing original codes that are duplicates in the hierarchy and can temporarily removed when building a table with this dimensional variable }
#' \item{slot \code{dupsUp}:}{character vector (or NULL) with original codes that are the corresponding upper-levels to the codes that may be removed because they are duplicates and that are listed in slot \code{dups}}
#' }
#' @name dimVar-class
#' @rdname dimVar-class
#' @exportClass dimVar
#' @note objects of class \code{dimVar} form the base for elements in slot \code{dimInfo} of class \code{dimInfo}.
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setClass(
  Class='dimVar',
  representation=representation(
    codesOriginal='characterOrNULL',
    codesDefault='characterOrNULL',
    codesMinimal='logicalOrNULL',
    vName='characterOrNULL',
    levels='numericOrNULL',
    structure='numericOrNULL',
    dims='listOrNULL',
    dups='characterOrNULL',
    dupsUp='characterOrNULL'
  ),
  prototype=prototype(
    codesOriginal=NULL,
    codesDefault=NULL,
    codesMinimal=NULL,
    vName=NULL,
    levels=NULL,
    structure=NULL,
    dims=NULL,
    dups=NULL,
    dupsUp=NULL
  ),
  validity=function(object) {
    if ( length(g_original_codes(object)) != length(g_default_codes(object)) ) {
      stop("dimVar:: length of 'codesOriginal' and 'codesDefault' differ!\n")
    }
    if ( length(g_varname(object)) != 1 ) {
      stop("dimVar:: length of 'vName' must equal 1!\n")
    }
    return(TRUE)
  }
)

#' S4 class describing a problemInstance-object
#'
#' An object of class \code{problemInstance} holds the main information that is
#' required to solve the secondary cell suppression problem.
#'
#' \describe{
#' \item{slot \code{strID}:}{a character vector (or NULL) of ID's identifying table cells}
#' \item{slot \code{Freq}:}{a numeric vector (or NULL) of counts for each table cell}
#' \item{slot \code{w}:}{a numeric vector (or NULL) of weights that should be used when solving the secondary cell suppression problem }
#' \item{slot \code{numVars}:}{a list (or NULL) with each element being a numeric vector holding values of specified numerical variables for each table cell}
#' \item{slot \code{lb}:}{numeric vector (or NULL) holding assumed lower bounds for each table cell}
#' \item{slot \code{ub}:}{numeric vector (or NULL) holding assumed upper bounds for each table cell}
#' \item{slot \code{LPL}:}{numeric vector (or NULL) holding required lower protection levels for each table cell}
#' \item{slot \code{UPL}:}{numeric vector (or NULL) holding required upper protection levels for each table cell}
#' \item{slot \code{SPL}:}{numeric vector (or NULL) holding required sliding protection levels for each table cell}
#' \item{slot \code{sdcStatus}:}{character vector (or NULL) holding the current anonymization state for each cell.
#' \itemize{
#' \item \code{z}: cell is forced to be published and must not be suppressed
#' \item \code{u}: cell has been primary suppressed
#' \item \code{x}: cell is a secondary suppression
#' \item \code{s}: cell can be published}
#' }}
#' @name problemInstance-class
#' @rdname problemInstance-class
#' @exportClass problemInstance
#' @note objects of class \code{problemInstance} are used as input for slot \code{problemInstance} in class \code{sdcProblem}
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setClass(
  Class='problemInstance',
  representation=representation(
    strID='characterOrNULL',
    Freq='numericOrNULL',
    w='numericOrNULL',
    numVars='listOrNULL',
    lb='numericOrNULL',
    ub='numericOrNULL',
    LPL='numericOrNULL',
    UPL='numericOrNULL',
    SPL='numericOrNULL',
    sdcStatus='characterOrNULL'
  ),
  prototype=prototype(
    strID=NULL,
    Freq=NULL,
    w=NULL,
    numVars=NULL,
    lb=NULL,
    ub=NULL,
    LPL=NULL,
    UPL=NULL,
    SPL=NULL,
    sdcStatus=NULL
  ),
  validity=function(object) {
    if ( !all.equal(
      length(g_strID(object)),
      length(g_freq(object)),
      length(g_lb(object)),
      length(g_ub(object)),
      length(g_SPL(object)),
      length(g_LPL(object)),
      length(g_UPL(object)),
      length(g_sdcStatus(object))) ) {
        stop("problemInstance:: slots 'strID', 'freq', 'lb', 'ub', 'SPL', 'LPL', 'SPL' and 'sdcStatus' must have the same length!\n")
    }
    if ( !is.null(g_numVars(object)) ) {
      if ( !all(sapply(g_numVars(object), length) == length(g_numVars(object)[[1]])) ) {
        stop("problemInstance:: length of vectors in slot 'numVars' differ in length!\n")
      }
    }
    if ( !is.null(g_numVars(object)[[1]]) & length(g_freq(object)) != length(g_numVars(object)[[1]]) ) {
      stop("problemInstance:: parameter 'Freq' and 'numVars' differ in length!\n")
    }
    if ( !is.null(g_w(object)) & length(g_freq(object)) != length(g_w(object)) ) {
      stop("problemInstance:: parameter 'Freq' and 'w' differ in length!\n")
    }
    if ( any(g_lb(object) <= g_freq(object) - g_LPL(object)) == FALSE ) {
      stop("problemInstance:: parameter 'lb' <= 'Freq'-'LPL' in some cases!\n")
    }
    if ( any(g_freq(object) - g_LPL(object) <= g_freq(object)) == FALSE ) {
      stop("problemInstance:: parameter 'Freq'-'LPL <= 'Freq' in some cases!\n")
    }
    if ( any(g_freq(object) <= g_freq(object) + g_UPL(object)) == FALSE ) {
      stop("problemInstance:: parameter 'Freq' <= 'Freq'+'UPL in some cases!\n")
    }
    if ( any(g_freq(object) + g_UPL(object) <= g_ub(object)) == FALSE) {
      stop("problemInstance:: parameter 'Freq'+'UPL' <= 'ub' in some cases!\n")
    }
    if ( any(g_ub(object) - g_lb(object) >= g_SPL(object)) == FALSE) {
      stop("problemInstance:: parameter 'ub'-'lb' >= 'SPL' in some cases!\n")
    }
    if ( !all(g_sdcStatus(object) %in% c('w','s','u','x','z')) ) {
      stop("problemInstance:: valid codes for sdcStatus are 'w', 'z', 's', 'x' or 'u'!\n")
    }
    return(TRUE)
  }
)

setClassUnion('dataObjOrNULL', c('dataObj', 'NULL'))
setClassUnion('dimInfoOrNULL', c('dimInfo', 'NULL'))
setClassUnion('problemInstanceOrNULL', c('problemInstance', 'NULL'))

#' S4 class describing a sdcProblem-object
#'
#' An object of class \code{sdcProblem} contains the entire information that is
#' required to protect the complete table that is given by the dimensional
#' variables. Such an object holds the data itself  (slot \code{dataObj}), the
#' entire information about the dimensional variables (slot \code{dimInfo}),
#' information on all table cells (ID's, bounds, values, anonymization state in
#' slot \code{problemInstance}), the indices on the subtables that need to be
#' considered if one wants to protect primary sensitive cells using a heuristic
#' approach (slot \code{partition}, information on which groups or rather
#' subtables have already been protected while performing a heuristic method
#' (slots \code{startI} and \code{startJ}) and the time that has been elapsed
#' (slot \code{elapsedTime}).
#'
#' \describe{
#' \item{slot \code{dataObj}:}{an object of class \code{dataObj} (or NULL) holding information on the underlying data}
#' \item{slot \code{dimInfo}:}{an object of class \code{dimInfo} (or NULL) containing information on all dimensional variables}
#' \item{slot \code{problemInstance}:}{an object of class \code{problemInstance} holding information on values, bounds, required protection levels as well as the anonymization state for all table cells}
#' \item{slot \code{partition}:}{a list object (or NULL) that is typically generated with calc.multiple(type='makePartitions',...) specifying information on the subtables and the necessary order that need to be protected when using a heuristic approach to solve the cell suppression problem}
#' \item{slot \code{startI}:}{a numeric vector of length 1 defining the group-level of the subtables in which a heuristic algorithm needs to start. All subtables having a group-index less than \code{startI} have already been protected}
#' \item{slot \code{startJ}:}{a numeric vector of length 1 defining the number of the table within the group defined by parameter \code{startI} at which a heuristic algorithm needs to start. All tables in the group having an index \code{j} smaller than \code{startJ} have already been protected}
#' \item{slot \code{indicesDealtWith}:}{a numeric vector holding indices of table cells that have protected and whose anonymization state must remain fixed}
#' \item{slot \code{elapsedTime}:}{a numeric vector of length 1 holding the time that has already been elapsed during the anonymization process}
#'  }
#' @name sdcProblem-class
#' @rdname sdcProblem-class
#' @exportClass sdcProblem
#' @note objects of class \code{sdcProblem} are typically generated by function \code{\link{makeProblem}} and are the input of functions \code{\link{primarySuppression}} and \code{\link{protectTable}}
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setClass(
  Class='sdcProblem',
  representation=representation(
    dataObj='dataObjOrNULL',
    dimInfo='dimInfoOrNULL',
    problemInstance='problemInstanceOrNULL',
    partition='listOrNULL',
    startI='numeric',
    startJ='numeric',
    indicesDealtWith='numericOrNULL',
    elapsedTime='numericOrNULL'
  ),
  prototype=prototype(
    dataObj=NULL,
    dimInfo=NULL,
    problemInstance=NULL,
    partition=NULL,
    startI=1,
    startJ=1,
    indicesDealtWith=NULL,
    elapsedTime=NULL
  ),
  validity=function(object) {
    if ( g_startI(object) > g_partition(object)$nrGroups ) {
      stop("argument 'startI' must be <=",g_partition(object)$nrGroups,"!\n")
    }
    if ( g_startJ(object) > length(g_partition(object)$indices[[g_startI(object)]]) ) {
      stop("argument 'startJ' must be <=",length(g_partition(object)$indices[[g_startI(object)]]),"!\n")
    }
    if ( length(g_startI(object)) != 1 ) {
      stop("sdcProblem:: length of argument 'startI' must equal 1!\n")
    }
    if ( length(g_startJ(object)) != 1 ) {
      stop("sdcProblem:: length of argument 'startJ' must equal 1!\n")
    }
    if ( g_startI(object) < 1 ) {
      stop("sdcProblem:: argument 'startI' must be >= 1!\n")
    }
    if ( g_startJ(object) < 1 ) {
      stop("sdcProblem:: argument 'startJ' must be >= 1!\n")
    }
    if ( length(g_elapsedTime(object)) != 1 ) {
      stop("sdcProblem:: length of argument 'elapsedTime' must equal 1!\n")
    }
    return(TRUE)
  }
)

#' S4 class describing a simpleTriplet-object
#'
#' Objects of class \code{simpleTriplet} define matrices that are stored in a
#' sparse format. Only the row- and column indices and the corresponding values
#' of non-zero cells are stored. Additionally, the dimension of the matrix given
#' by the total number of rows and columns is stored.
#'
#' \describe{
#' \item{slot \code{i}:}{a numeric vector specifying row-indices with each value being geq 1 and leq of the value in \code{nrRows}}
#' \item{slot \code{j}:}{a numeric vector specifying column-indices with each value being geq 1 and leq of the value in \code{nrCols}}
#' \item{slot \code{v}:}{a numeric vector specifying the values of the matrix in cells specified by the corresponding row- and column indices}
#' \item{slot \code{nrRows}:}{a numeric vector of length 1 holding the total number of rows of the matrix}
#' \item{slot \code{nrCols}:}{a numeric vector of length 1 holding the total number of columns of the matrix}
#' }
#' @name simpleTriplet-class
#' @rdname simpleTriplet-class
#' @exportClass simpleTriplet
#' @note objects of class \code{simpleTriplet} are input of slot \code{constraints} in class \code{\link{linProb-class}} and slot slot \code{con} in class \code{\link{cutList-class}}
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setClass(
  Class='simpleTriplet',
  representation=representation(
    i='numericOrNULL',
    j='numericOrNULL',
    v='numericOrNULL',
    nrRows='numericOrNULL',
    nrCols='numericOrNULL'
  ),
  prototype=prototype(
    i=numeric(0),
    j=numeric(0),
    v=numeric(0),
    nrRows=numeric(0),
    nrCols=numeric(0)
  ),
  validity=function(object) {
    if ( length(object@i) != length(object@j) ) {
      stop("simpleTriplet:: length of 'i' and 'j' differ!\n")
    }
    if ( length(object@i) != length(object@v) ) {
      stop("simpleTriplet:: length of 'i' and 'v' differ!\n")
    }
    if ( length(object@nrRows) + length(object@nrCols) != 2 ) {
      stop("simpleTriplet:: 'nrRows' and 'nrCols' must be a vector of length 1!\n")
    }
    return(TRUE)
  }
)
setClassUnion("simpleTripletOrNULL", c("simpleTriplet", "NULL"))

#' S4 class describing a linProb-object
#'
#' An object of class \code{linProb} defines a linear problem given by the
#' objective coefficients (slot \code{objective}), a constraint matrix (slot
#' \code{constraints}), the direction (slot \code{direction}) and the right
#' hand side (slot \code{rhs}) of the constraints. Also, allowed lower (slot
#' \code{boundsLower}) and upper (slot \code{boundsUpper}) bounds of the
#' variables as well as its types (slot \code{types}) are specified.
#'
#' \describe{
#' \item{slot \code{objective}:}{a numeric vector holding coefficients of the objective function}
#' \item{slot \code{constraints}:}{an object of class \code{\link{simpleTriplet-class}} specifying the constraint matrix of the problem}
#' \item{slot \code{direction}:}{a character vector holding the directions of the constraints, allowed values are:
#' \itemize{
#' \item \code{==}: equal
#' \item \code{<}: less
#' \item \code{>}: greater
#' \item \code{<=}: less or equal
#' \item \code{>=}: greater or equal}
#' }
#' \item{slot \code{rhs}:}{numeric vector holding right hand side values of the constraints}
#' \item{slot \code{boundsLower}:}{a numeric vector holding lower bounds of the objective variables}
#' \item{slot \code{boundsUpper}:}{a numeric vector holding upper bounds of the objective variables}
#' \item{slot \code{types}:}{a character vector specifying types of the objective variables, allowed types are:
#' \itemize{
#' \item \code{C}: binary
#' \item \code{B}: continuous
#' \item \code{I}: integer}
#' }
#' }
#' @name linProb-class
#' @rdname linProb-class
#' @exportClass linProb
#' @note when solving the problems in the procedure, minimization of the objective is performed.
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setClass(
  Class='linProb',
  representation=representation(
    objective='numericOrNULL',
    constraints='simpleTripletOrNULL',
    direction='characterOrNULL',
    rhs='numericOrNULL',
    boundsLower='listOrNULL',
    boundsUpper='listOrNULL',
    types='characterOrNULL'
  ),
  prototype=prototype(
    objective=NULL,
    constraints=NULL,
    direction=NULL,
    rhs=NULL,
    boundsLower=NULL,
    boundsUpper=NULL,
    types=NULL
  ),
  validity=function(object) {
    if ( length(object@rhs) != length(object@direction) ) {
      stop("linProb:: length of 'rhs' and 'direction' differ!\n")
    }
    nrRows.constraints <- g_nr_rows(object@constraints)
    if ( length(object@direction) != nrRows.constraints ) {
      stop("linProb:: length of 'direction' and number of rows of 'constraints' differ!\n")
    }
    nrCols.constraints <- g_nr_cols(object@constraints)
    if ( length(object@objective) != nrCols.constraints ) {
      stop("linProb:: length of 'objective' and number of columns of 'constraints' differ!\n")
    }
    if ( length(object@objective) != length(object@types) ) {
      stop("linProb:: Length of 'objective' and 'types' differ!\n")
    }
    if ( !all(object@boundsLower$indices %in% 1:length(object@direction)) ) {
      stop("linProb:: wrong indices of 'boundsLower!'\n")
    }
    if ( !all(object@boundsUpper$indices %in% 1:length(object@direction)) ) {
      stop("linProb:: wrong indices of 'boundsUpper!'\n")
    }
    if ( length(object@boundsLower$indices) != length(object@boundsLower$value) ) {
      stop("linProb:: length of indices and values in 'boundsLower' differ!\n")
    }
    if ( length(object@boundsUpper$indices) != length(object@boundsUpper$value) ) {
      stop("linProb:: length of indices and values in 'boundsUpper' differ!\n")
    }
    if ( !all(object@direction %in% c("==","<",">",">=","<=")) ) {
      stop("linProb:: illegal symbols in 'direction' differ!\n")
    }
    return(TRUE)
  }
)

#' S4 class describing a cutList-object
#'
#' An object of class \code{cutList} holds constraints that can be extracted and
#' used as for objects of class \code{\link{linProb-class}}. An object of class
#' \code{cutList} consists of a constraint matrix (slot \code{con}), a vector
#' of directions (slot \code{direction}) and a vector specifying the right hand
#' sides of the constraints (slot \code{rhs}).
#'
#' \describe{
#' \item{slot \code{con}:}{an object of class \code{\link{simpleTriplet-class}} specifying the constraint matrix of the problem}
#' \item{slot \code{direction}:}{a character vector holding the directions of the constraints, allowed values are:
#' \itemize{
#' \item \code{==}: equal
#' \item \code{<}: less
#' \item \code{>}: greater
#' \item \code{<=}: less or equal
#' \item \code{>=}: greater or equal}
#' }
#' \item{slot \code{rhs}:}{numeric vector holding right hand side values of the constraints}
#' }
#' @name cutList-class
#' @rdname cutList-class
#' @exportClass cutList
#' @note objects of class \code{cutList} are dynamically generated (and removed) during the cut and branch algorithm when solving the secondary cell suppression problem
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setClass(
  Class='cutList',
  representation=representation(
    con='simpleTriplet',
    direction='character',
    rhs='numeric'
  ),
  prototype=prototype(
    con=new("simpleTriplet"),
    direction=character(0),
    rhs=numeric(0)
  ),
  validity=function(object) {
    if ( g_nr_rows(g_constraints(object)) != length(g_direction(object)) ) {
      stop("cutList:: number of rows of 'con' and length of 'direction' differs!\n")
    }
    if ( length(g_direction(object))!= length(g_rhs(object)) ) {
      stop("cutList:: length of 'direction' and 'rhs' differ!\n")
    }
    if ( !all(g_direction(object) %in% c(">", ">=", "==", "<", "<=")) ) {
      stop("cutList:: elements of 'direction' must only contain symbols '>', '>=', '==', '<' or '<='!\n")
    }
    return(TRUE)
  }
)

#' S4 class describing a safeObj-object
#'
#' Objects of class \code{safeObj} are the final result after protection a
#' tabular structure. After a successful run of \code{\link{protectTable}}
#' an object of this class is generated and returned. Objects of class
#' \code{safeObj} contain a final, complete data set (slot \code{finalData})
#' that has a column showing the anonymization state of each cell and the
#' complete information on the dimensional variables that have defined the table
#' that has been protected (slot \code{dimInfo}). Also, the number of
#' non-duplicated table cells (slot \code{nrNonDuplicatedCells}) is returned
#' along with the number of primary (slot \code{nrPrimSupps}) and secondary
#' (slot \code{nrSecondSupps}) suppressions. Furthermore, the number of cells
#' that can be published (slot \code{nrPublishableCells}), the algorithm that
#' has  been used to protect the data (slot \code{suppMethod}) and the time that
#' was needed to protect the data structure (slot \code{elapsedTime}) is
#' returned.
#'
#' \describe{
#' \item{slot \code{finalData}:}{a data.frame (or NULL) featuring columns for each variable defining the table (with their original codes), the cell counts and values of any numerical variables and the anonymization status for each cell with
#' \itemize{
#' \item \code{s, z}: cell can be published
#' \item \code{u}: cell is a primary sensitive cell
#' \item \code{x}: cell was selected as a secondary suppression}
#' }
#' \item{slot \code{dimInfo}:}{an object of class \code{\link{dimInfo-class}} holding all information on variables defining the table }
#' \item{slot \code{nrNonDuplicatedCells}:}{numeric vector of length 1 (or NULL) showing the number of non-duplicated table cells. This value is different from 0 if any dimensional variable features duplicated codes. These codes have been re-added to the final dataset.}
#' \item{slot \code{nrPrimSupps}:}{numeric vector of length 1 (or NULL) showing the number of primary suppressed cells}
#' \item{slot \code{nrSecondSupps}:}{numeric vector of length 1 (or NULL) showing the number of secondary suppressions}
#' \item{slot \code{nrPublishableCells}:}{numeric vector of length 1 (or NULL) showing the number of cells that may be published}
#' \item{slot \code{suppMethod}:}{character vector of length 1 holding information on the protection method}
#' \item{slot \code{elapsedTime}:}{numeric vector of length 1 holding the time that was required to protect the table}
#' }
#' @name safeObj-class
#' @rdname safeObj-class
#' @exportClass safeObj
#' @note objects of class \code{safeObj} are returned after the function \code{\link{protectTable}} has finished.
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setClass(
  Class='safeObj',
  representation=representation(
    finalData='dataframeOrNULL',
    dimInfo='dimInfoOrNULL',
    nrNonDuplicatedCells='numericOrNULL',
    nrPrimSupps='numericOrNULL',
    nrSecondSupps='numericOrNULL',
    nrPublishableCells='numericOrNULL',
    suppMethod='characterOrNULL',
    elapsedTime='numericOrNULL'
  ),
  prototype=prototype(
    finalData=NULL,
    dimInfo=NULL,
    nrNonDuplicatedCells=NULL,
    nrPrimSupps=NULL,
    nrSecondSupps=NULL,
    nrPublishableCells=NULL,
    suppMethod=NULL,
    elapsedTime=NULL
  ),
  validity=function(object) {
    if ( length(g_nrPrimSupps(object)) != 1 ) {
      stop("safeObj:: length of 'nrPrimSupps' must equal 1!\n")
    }
    if ( length(g_nrNonDuplicatedCells(object)) != 1 ) {
      stop("safeObj:: length of 'nrNonDuplicatedCells' must equal 1!\n")
    }
    if ( length(g_nrSecondSupps(object)) != 1 ) {
      stop("safeObj:: length of 'nrSecondSupps' must equal 1!\n")
    }
    if ( length(g_nrPublishableCells(object)) != 1 ) {
      stop("safeObj:: length of 'nrPublishableCells' must equal 1!\n")
    }
    if ( length(g_suppMethod(object)) != 1 ) {
      stop("safeObj:: length of 'suppMethod' must equal 1!\n")
    }
    if ( !g_suppMethod(object) %in% c('SIMPLEHEURISTIC', 'HITAS', 'OPT', 'HYPERCUBE') ) {
      stop("safeObj:: 'suppMethod' must bei either 'SIMPLEHEURISTIC', 'HITAS', 'HYPERCUBE' or 'OPT'!\n")
    }
    if ( length(g_elapsedTime(object)) != 1 ) {
      stop("safeObj:: length of 'elapsedTime' must equal 1!\n")
    }
    return(TRUE)
  }
)
