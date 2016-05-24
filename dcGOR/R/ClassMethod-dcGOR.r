require(Matrix)
################################################################################
################################################################################
#' @title Definition for S4 class InfoDataFrame
#' @description \code{InfoDataFrame} has two slots: data and dimLabels.
#' @return Class InfoDataFrame
#' @slot data A data.frame containing terms (rows) and measured variables (columns).
#' @slot dimLabels A character descripting labels for rows and columns.
#' @section Creation:
#' An object of this class can be created via: \code{new("InfoDataFrame", data, dimLabels)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{dim()}: }{retrieve the dimension in the object}
#' \item{\code{nrow()}: }{retrieve number of rows in the object}
#' \item{\code{ncol()}: }{retrieve number of columns in the object}
#' \item{\code{rowNames()}: }{retrieve names of rows in the object}
#' \item{\code{colNames()}: }{retrieve names of columns in the object}
#' \item{\code{dimLabels()}: }{retrieve the slot 'dimLabels', containing labels used for display of rows and columns in the object}
#' \item{\code{Data()}: }{retrieve the slot 'data' in the object}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' \item{\code{as(data.frame, "InfoDataFrame")}: }{convert a data.frame to an object of class InfoDataFrame}
#' \item{\code{[i,j]}: }{get the subset of the same class}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("InfoDataFrame")}: }{show the class definition}
#' \item{\code{showMethods(classes="InfoDataFrame")}: }{show the method definition upon this class}
#' \item{\code{getSlots("InfoDataFrame")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("InfoDataFrame")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="InfoDataFrame")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @import Matrix
#' @import igraph
#' @docType class
#' @keywords S4 classes
#' @name InfoDataFrame-class
#' @rdname InfoDataFrame-class
#' @seealso \code{\link{InfoDataFrame-method}}
#' @examples
#' # generate data on domain information on
#' data <- data.frame(x=1:10, y=I(LETTERS[1:10]), row.names=paste("Domain", 1:10, sep="_"))
#' dimLabels <- c("rowLabels", "colLabels")
#' # create an object of class InfoDataFrame
#' x <- new("InfoDataFrame", data=data, dimLabels=dimLabels)
#' x
#' # alternatively, using coerce methods
#' x <- as(data, "InfoDataFrame")
#' x
#' # look at various methods defined on class Anno
#' dimLabels(x)
#' dim(x)
#' nrow(x)
#' ncol(x)
#' rowNames(x)
#' colNames(x)
#' Data(x)
#' x[1:3,]

#' @rdname InfoDataFrame-class
#' @aliases InfoDataFrame
#' @exportClass InfoDataFrame
setClass(
    Class="InfoDataFrame",
    representation(
        data = "data.frame",
        dimLabels = "character"
    ),
    prototype = prototype(
        data = new( "data.frame" ),
        dimLabels=c("rowLabels", "colLabels")
    ),
    validity = function(object){
        if(!is.data.frame(object@data)){
            return("data is not data.frame")
        }else{
            return(TRUE)
        }
    } 
)

########################################
#' @title Methods defined for S4 class InfoDataFrame
#' @description Methods defined for class \code{InfoDataFrame}.
#' @param x an object of class \code{InfoDataFrame}
#' @param object an object of class \code{InfoDataFrame}
#' @param i an index
#' @param j an index
#' @param drop a logic for matrices and arrays. If TRUE the result is coerced to the lowest possible dimension. This only works for extracting elements, not for the replacement
#' @param ... additional parameters
#' @docType methods
#' @keywords S4 methods
#' @name InfoDataFrame-method
#' @rdname InfoDataFrame-method
#' @seealso \code{\link{InfoDataFrame-class}}

setGeneric("dimLabels", function(x) standardGeneric("dimLabels"))
#' @rdname InfoDataFrame-method
#' @aliases dimLabels
#' @export
setMethod("dimLabels", "InfoDataFrame", function(x) x@dimLabels)

#' @rdname InfoDataFrame-method
#' @aliases dim,InfoDataFrame-method
#' @export
setMethod("dim", "InfoDataFrame", function(x) dim(x@data))

#' @rdname InfoDataFrame-method
#' @aliases nrow
#' @export
setMethod("nrow", "InfoDataFrame", function(x) nrow(x@data))

#' @rdname InfoDataFrame-method
#' @aliases ncol
#' @export
setMethod("ncol", "InfoDataFrame", function(x) ncol(x@data))

setGeneric("rowNames", function(x) standardGeneric("rowNames"))
#' @rdname InfoDataFrame-method
#' @aliases rowNames
#' @export
setMethod("rowNames", "InfoDataFrame", function(x) rownames(x@data))

setGeneric("colNames", function(x) standardGeneric("colNames"))
#' @rdname InfoDataFrame-method
#' @aliases colNames
#' @export
setMethod("colNames", "InfoDataFrame", function(x) colnames(x@data))

setGeneric("Data", function(x) standardGeneric("Data"))
#' @rdname InfoDataFrame-method
#' @aliases Data
#' @export
setMethod("Data", "InfoDataFrame", function(x) x@data)

#' @rdname InfoDataFrame-method
#' @name data.frame2InfoDataFrame
setAs("data.frame", "InfoDataFrame", function(from) new("InfoDataFrame",data=from))

.wrapcat <- function(lbl, nms, total, ..., indent=2, exdent=4)
{
    lbl <- sprintf("%s:", lbl)
    txt <- paste(c(lbl,  nms), collapse=" ")
    ext <-
        if (length(nms) < total) sprintf("(%d total)", total)
        else character()
    txt <- paste(c(lbl,  nms, ext), collapse=" ")
    cat(strwrap(txt, ..., indent=indent, exdent=exdent), sep="\n")
}
.selectSomeIndex <- function(x, maxToShow=5, byrow=TRUE, ...)
{
    len <-
        if (byrow) dim(x)[[1]]
        else dim(x)[[2]]
    if (maxToShow < 3) maxToShow <- 3
    if (len > maxToShow) {
        maxToShow <- maxToShow - 1
        bot <- ceiling(maxToShow/2)
        top <- len-(maxToShow-bot-1)
        list(1:bot, "...", top:len)
    } else if (len >= 1) {
        list(1:len, NULL, NULL)
    }else{
        list(NULL, NULL, NULL)
    }
}
.showInfoDataFrame <- function(object, labels=list(0)) 
{
    lbls <- list(
        object=class(object),
        termNames=dimLabels(object)[[1]],
        varLabels="varLabels"
    )
    lbls[names(labels)] <- labels
    ## create a simplified object for extracting names
    idx <- .selectSomeIndex(Data(object), maxToShow=6)
    idy <- .selectSomeIndex(Data(object), byrow=FALSE, maxToShow=6)
    Data <- Data(object)[c(idx[[1]], idx[[3]]), c(idy[[1]], idy[[3]]), drop=FALSE]
    rnms <- rownames(Data)
    nms <- c(rnms[idx[[1]]], idx[[2]], if (!is.null(idx[[1]])) rnms[-idx[[1]]] else NULL)
    ## for terms
    .wrapcat(lbls$termNames, nms, nrow(object))
    ## for domains
    cnms <- colnames(Data)
    vars <- c(cnms[idy[[1]]], idy[[2]], cnms[-idy[[1]]])
    .wrapcat(lbls$varLabels, vars, ncol(object))
}
#' @rdname InfoDataFrame-method
#' @export
setMethod("show", "InfoDataFrame",
    function(object) {
        cat("An object of S4 class '", class(object), "'\n", sep="")
        if(sum(dim(object@data)) !=0 ){
            .showInfoDataFrame(
                object, 
                labels=list(
                    object="InfoDataFrame",
                    termNames="rowNames",
                    varLabels="colNames"
                )
            )
        }else{
            cat("but contains no data\n", sep="")
        }
    }
)

#' @rdname InfoDataFrame-method
#' @aliases [,InfoDataFrame-method
#' @export
setMethod("[", signature(x="InfoDataFrame"),
    function(x, i, j, ..., drop = FALSE) {
        if (missing(drop)){
            drop = FALSE
        }
        
        if(missing(j)) {
            labels <- x@dimLabels
            D <- x@data[i,,drop = drop]
        } else {
            labels <- x@dimLabels[j,,drop = drop]
        }
        
        if( missing( i )){
            D <- x@data[,j,drop = drop]
        }else{
            D <- x@data[i,j,drop = drop]
        }
        
        x <- new("InfoDataFrame", data=D, dimLabels=labels)
    }
)

################################################################################
################################################################################
#' @title Definition for VIRTUAL S4 class AnnoData
#' @description \code{AnnoData} is union of other classes: either matrix or dgCMatrix (a sparse matrix in the package Matrix). It is used as a virtual class
#' @return Class AnnoData
#' @import Matrix
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name AnnoData-class
#' @rdname AnnoData-class
#' @seealso \code{\link{Anno-class}}

#' @rdname AnnoData-class
#' @aliases AnnoData
#' @exportClass AnnoData
setClassUnion("AnnoData", c("matrix", "dgCMatrix"))

################################################################################
################################################################################
#' @title Definition for S4 class Anno
#' @description \code{Anno} has 3 slots: annoData, termData and domainData
#' @return Class Anno
#' @slot annoData An object of S4 class \code{\link{AnnoData}}, containing data matrix with the column number equal to nrow(termData) and the row number equal to nrow(domainData).
#' @slot termData An object of S4 class \code{\link{InfoDataFrame}}, describing information on columns in annoData.
#' @slot domainData An object of S4 class \code{\link{InfoDataFrame}}, describing information on rows in annoData.
#' @section Creation:
#' An object of this class can be created via: \code{new("Anno", annoData, termData, domainData)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{dim()}: }{retrieve the dimension in the object}
#' \item{\code{annoData()}: }{retrieve the slot 'annoData' in the object}
#' \item{\code{termData()}: }{retrieve the slot 'termData' (as class InfoDataFrame) in the object}
#' \item{\code{domainData()}: }{retrieve the slot 'domainData' (as class InfoDataFrame) in the object}
#' \item{\code{tData()}: }{retrieve termData (as data.frame) in the object}
#' \item{\code{dData()}: }{retrieve domainData (as data.frame) in the object}
#' \item{\code{termNames()}: }{retrieve term names (ie, row names of termData) in the object}
#' \item{\code{domanNames()}: }{retrieve domain names (ie, row names of domainData) in the object}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' \item{\code{as(matrix, "Anno")}: }{convert a matrix to an object of class Anno}
#' \item{\code{as(dgCMatrix, "Anno")}: }{convert a dgCMatrix (a sparse matrix) to an object of class Anno}
#' \item{\code{[i,j]}: }{get the subset of the same class}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("Anno")}: }{show the class definition}
#' \item{\code{showMethods(classes="Anno")}: }{show the method definition upon this class}
#' \item{\code{getSlots("Anno")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("Anno")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="Anno")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name Anno-class
#' @seealso \code{\link{Anno-method}}
#' @examples
#' # create an object of class Anno, only given a matrix
#' annoData <- matrix(runif(50),nrow=10,ncol=5)
#' as(annoData, "Anno")
#'
#' # create an object of class Anno, given a matrix plus information on its columns/rows
#' # 1) create termData: an object of class InfoDataFrame
#' data <- data.frame(x=1:5, y=I(LETTERS[1:5]), row.names=paste("Term", 1:5, sep="_"))
#' termData <- new("InfoDataFrame", data=data)
#' termData
#' # 2) create domainData: an object of class InfoDataFrame
#' data <- data.frame(x=1:10, y=I(LETTERS[1:10]), row.names=paste("Domain", 1:10, sep="_"))
#' domainData <- new("InfoDataFrame", data=data)
#' domainData
#' # 3) create an object of class Anno
#' # VERY IMPORTANT: make sure having consistent names between annoData and domainData (and termData)
#' annoData <- matrix(runif(50),nrow=10,ncol=5)
#' rownames(annoData) <- rowNames(domainData)
#' colnames(annoData) <- rowNames(termData)
#' x <- new("Anno", annoData=annoData, domainData=domainData, termData=termData)
#' x
#' # 4) look at various methods defined on class Anno
#' dim(x)
#' annoData(x)
#' termData(x)
#' tData(x)
#' domainData(x)
#' dData(x)
#' termNames(x)
#' domainNames(x)
#' # 5) get the subset
#' x[1:3,1:2]

#' @rdname Anno-class
#' @aliases Anno
#' @exportClass Anno
setClass(
    Class="Anno",
    representation(
        annoData = "AnnoData",
        termData = "InfoDataFrame",
        domainData = "InfoDataFrame"
    ),
    prototype = prototype(
        annoData = matrix(),
        termData = new("InfoDataFrame",dimLabels=c("termNames", "termColumns")),
        domainData = new("InfoDataFrame",dimLabels=c("domainNames", "domainColumns"))
    ),
    validity = function(object){
        msg <- NULL
        # dimension for annoData
        adim <- dim(object)
        ## annoData and domainData
        if( dim(domainData(object))[1] != 0 ){
            if (adim[1] != dim(domainData(object))[1]){
                msg <- append(msg, "domain numbers differ between annoData and domainData")
            }
            if (!identical(domainNames(object), rowNames(domainData(object)))){
                msg <- append(msg, "domain names differ between annoData and domainData")
            }
        }
        ## annoData and termData
        if( dim(termData(object))[1] != 0 ){
            if (adim[2] != dim(termData(object))[1]){
                msg <- append(msg, "term numbers differ between annoData and termData")
            }
            if (!identical(termNames(object), rowNames(termData(object)))){
                msg <- append(msg, "term names differ between annoData and termData")
            }
        }
        if (is.null(msg)) TRUE else msg
    }
)

########################################
#' @title Methods defined for S4 class Anno
#' @description Methods defined for class \code{Anno}.
#' @param x an object of class \code{Anno}
#' @param object an object of class \code{Anno}
#' @param i an index
#' @param j an index
#' @param drop a logic for matrices and arrays. If TRUE the result is coerced to the lowest possible dimension. This only works for extracting elements, not for the replacement
#' @param ... additional parameters
#' @docType methods
#' @keywords S4 methods
#' @name Anno-method
#' @rdname Anno-method
#' @seealso \code{\link{Anno-class}}

#' @rdname Anno-method
#' @aliases dim,Anno-method
#' @export
setMethod("dim", "Anno", function(x) dim(x@annoData))

setGeneric("annoData", function(x) standardGeneric("annoData"))
#' @rdname Anno-method
#' @aliases annoData
#' @export
setMethod("annoData", "Anno", function(x) x@annoData)

setGeneric("termData", function(x) standardGeneric("termData"))
#' @rdname Anno-method
#' @aliases termData
#' @export
setMethod("termData", "Anno", function(x) x@termData)

setGeneric("domainData", function(x) standardGeneric("domainData"))
#' @rdname Anno-method
#' @aliases domainData
#' @export
setMethod("domainData", "Anno", function(x) x@domainData)

setGeneric("tData", function(object) standardGeneric("tData"))
#' @rdname Anno-method
#' @aliases tData
#' @export
setMethod("tData", signature(object="Anno"), function(object){
    data <- Data(termData(object))
    if(sum(dim(data))==0){
        cat("No data is available\n", sep="")
    }else{
        data
    }
})

setGeneric("dData", function(object) standardGeneric("dData"))
#' @rdname Anno-method
#' @aliases dData
#' @export
setMethod("dData", signature(object="Anno"), function(object){
    data <- Data(domainData(object))
    if(sum(dim(data))==0){
        cat("No data is available\n", sep="")
    }else{
        data
    }
})

setGeneric("termNames", function(object) standardGeneric("termNames"))
#' @rdname Anno-method
#' @aliases termNames
#' @export
setMethod("termNames", signature(object="Anno"), function(object) rowNames(termData(object)))

setGeneric("domainNames", function(object) standardGeneric("domainNames"))
#' @rdname Anno-method
#' @aliases domainNames
#' @export
setMethod("domainNames", signature(object="Anno"), function(object) rowNames(domainData(object)))

#' @rdname Anno-method
#' @name matrix2Anno
setAs("matrix", "Anno", function(from) {
    ## for domainData
    rn <- rownames(from)
    if(is.null(rn)) rn <- 1:nrow(from)
    domainData <- new("InfoDataFrame", data=data.frame(names=rn))
    ## for termData    
    cn <- colnames(from)
    if(is.null(cn)) cn <- 1:ncol(from)
    termData <- new("InfoDataFrame", data=data.frame(names=cn))
    ## for Anno
    new("Anno", annoData=from, domainData=domainData, termData=termData)
})

#' @rdname Anno-method
#' @name dgCMatrix2Anno
setAs("dgCMatrix", "Anno", function(from) {
    ## for domainData
    rn <- rownames(from)
    if(is.null(rn)) rn <- 1:nrow(from)
    domainData <- new("InfoDataFrame", data=data.frame(names=rn))
    ## for termData    
    cn <- colnames(from)
    if(is.null(cn)) cn <- 1:ncol(from)
    termData <- new("InfoDataFrame", data=data.frame(names=cn))
    ## for Anno
    new("Anno", annoData=from, domainData=domainData, termData=termData)
})

#' @rdname Anno-method
#' @export
setMethod("show", 
    signature=signature(object="Anno"),
    function(object) {
        cat("An object of S4 class '", class(object), "'\n", sep="")
        adim <- dim(object)
        if (length(adim)>1){
            cat("@annoData:", if (length(adim)>1) paste(adim[[1]], "domains,",adim[[2]], "terms") else NULL, "\n")
        }
        ## termData
        if( dim(termData(object))[1] != 0 ){
            cat("@termData (", class(termData(object)), ")\n", sep="")
            .showInfoDataFrame(
                termData(object), 
                labels=list(
                    object="termData",
                    termNames="termNames",
                    varLabels="tvarLabels"
                )
            )
        }else{
            cat("@termData (NULL)\n", sep="")
        }
        ## domainData
        if( dim(domainData(object))[1] != 0 ){
            cat("@domainData (", class(domainData(object)), ")\n", sep="")
            .showInfoDataFrame(
                domainData(object), 
                labels=list(
                    object="domainData",
                    termNames="domainNames",
                    varLabels="dvarLabels"
                )
            )
        }else{
            cat("@domainData (NULL)\n", sep="")
        }
    }
)

#' @rdname Anno-method
#' @aliases [,Anno-method
#' @export
setMethod("[", signature(x="Anno"), 
    function(x, i, j, ..., drop = FALSE) {
        if (missing(drop)){
            drop <- FALSE
        }
        if (missing(i) && missing(j)) {
            if (length(list(...))!=0){
                stop("specify domains or terms to subset")
            }
            return(x)
        }

        if (!missing(j)) {
            tD <- termData(x)[j,, ..., drop=drop]
        }else{
            tD <- termData(x)
        }
        
        if (!missing(i)) {
            dD <- domainData(x)[i,,..., drop=drop]
        }else{
            dD <- domainData(x)
        }
        
        if (missing(i) & !missing(j)){
            aD <- annoData(x)[,j]
        }else if (!missing(i) & missing(j)){
            aD <- annoData(x)[i,]
        }else if (!missing(i) & !missing(j)){
            aD <- annoData(x)[i,j]
        }else{
            aD <- annoData(x)
        }
        
        x <- new("Anno", annoData=aD, termData=tD, domainData=dD)
    }
)

################################################################################
################################################################################
#' @title Definition for S4 class Eoutput
#' @description \code{Eoutput} is an S4 class to store output from enrichment analysis by \code{\link{dcEnrichment}}.
#' @return Class Eoutput
#' @slot domain A character specifying the domain identity
#' @slot ontology A character specifying the ontology identity
#' @slot term_info A data.frame of nTerm X 5 containing term information, where nTerm is the number of terms in consideration, and the 5 columns are "term_id" (i.e. "Term ID"), "term_name" (i.e. "Term Name"), "namespace" (i.e. "Term Namespace"), "distance" (i.e. "Term Distance") and "IC" (i.e. "Information Content for the term based on annotation frequency by it")
#' @slot anno A list of terms, each storing annotated domains (also within the background domains). Always, terms are identified by "term_id" and domain members identified by their ids (e.g. sunids for SCOP domains)
#' @slot data A vector containing input data in \code{\link{dcEnrichment}}. It is not always the same as the input data as only those mappable are retained
#' @slot background A vector containing background in \code{\link{dcEnrichment}}. It is not always the same as the input background (if provided by the user) as only those mappable are retained
#' @slot overlap A list of terms, each storing domains overlapped between domains annotated by a term and domains in the input data (i.e. the domains of interest). Always, terms are identified by "term_id" and domain members identified by their ids (e.g. sunids for SCOP domains)
#' @slot zscore A vector of terms, containing z-scores
#' @slot pvalue A vector of terms, containing p-values
#' @slot adjp A vector of terms, containing adjusted p-values. It is the p value but after being adjusted for multiple comparisons
#' @section Creation:
#' An object of this class can be created via: \code{new("Eoutput", domain, ontology, term_info, anno, data, overlap, zscore, pvalue, adjp)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{zscore()}: }{retrieve the slot 'zscore' in the object}
#' \item{\code{pvalue()}: }{retrieve the slot 'pvalue' in the object}
#' \item{\code{adjp()}: }{retrieve the slot 'adjp' in the object}
#' \item{\code{view()}: }{retrieve an integrated data.frame used for viewing the object}
#' \item{\code{write()}: }{write the object into a local file}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("Eoutput")}: }{show the class definition}
#' \item{\code{showMethods(classes="Eoutput")}: }{show the method definition upon this class}
#' \item{\code{getSlots("Eoutput")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("Eoutput")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="Eoutput")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name Eoutput-class
#' @rdname Eoutput-class
#' @seealso \code{\link{Eoutput-method}}
#' @examples
#' \dontrun{
#' # 1) load SCOP.sf (as 'InfoDataFrame' object)
#' SCOP.sf <- dcRDataLoader('SCOP.sf')
#' # randomly select 20 domains
#' data <- sample(rowNames(SCOP.sf), 20)
#' 
#' # 2) perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="GOMF")
#' eoutput
#'
#' # 3) write into the file 'Eoutput.txt' in your local directory
#' write(eoutput, file='Eoutput.txt')
#'
#' # 4) view the top 5 significant terms
#' view(eoutput, top_num=5, sortBy="pvalue", details=TRUE)
#'
#' # 4) retrieve several slots directly
#' zscore(eoutput)[1:5]
#' pvalue(eoutput)[1:5]
#' adjp(eoutput)[1:5]
#' }

#' @rdname Eoutput-class
#' @aliases Eoutput
#' @exportClass Eoutput
setClass(
    Class="Eoutput",
    representation(
        domain      = "character",
        ontology    = "character",
        term_info   = "data.frame",
        anno        = "list",
        data        = "vector",
        background  = "vector",
        overlap     = "vector",
        zscore      = "vector",
        pvalue      = "vector",
        adjp        = "vector"
    ),
    prototype = prototype(
        domain      = "domain",
        ontology    = "ontology",
        term_info   = data.frame(),
        anno        = list(),
        data        = vector(),
        background  = vector(),
        overlap     = vector(),
        zscore      = vector(),
        pvalue      = vector(),
        adjp        = vector()
    ),
    validity = function(object){
        if(!is.data.frame(object@term_info)){
            return("term_info is not data.frame")
        }else{
            return(TRUE)
        }
    }
)

########################################
#' @title Methods defined for S4 class Eoutput
#' @description Methods defined for S4 class \code{Eoutput}.
#' @param object an object of S4 class \code{Eoutput}. Usually this is an output from \code{\link{dcEnrichment}}
#' @param x an object of S4 class \code{Eoutput}. Usually this is an output from \code{\link{dcEnrichment}}
#' @param top_num the maximum number (5, by default) of terms will be viewed. If NULL or NA, all terms will be viewed (this can be used for the subsequent saving)
#' @param sortBy which statistics will be used for sorting and viewing terms. It can be "pvalue" for p value, "adjp" for adjusted p value, "zscore" for enrichment z-score, "nAnno" for the number in domains annotated by a term, "nOverlap" for the number in overlaps, and "none" for ordering simply according to ID of terms
#' @param decreasing logical to indicate whether to sort in a decreasing order. If it is null, by default it will be true for "zscore", "nAnno" or "nOverlap"; otherwise false
#' @param details logical to indicate whether the detailed information of terms is also viewed. By default, it sets to TRUE for the inclusion
#' @param file a character specifying a file name written into. By default, it is 'Eoutput.txt'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' view(x) returns a data frame with following components:
#' \itemize{
#'  \item{\code{term_id}: term ID}
#'  \item{\code{nAnno}: number in domains annotated by a term}
#'  \item{\code{nGroup}: number in domains from the input group}
#'  \item{\code{nOverlap}: number in overlaps}
#'  \item{\code{zscore}: enrichment z-score}
#'  \item{\code{pvalue}: p value}
#'  \item{\code{adjp}: adjusted p value}
#'  \item{\code{term_name}: term name}
#'  \item{\code{term_namespace}: term namespace; optional, it is only appended when "details" is true}
#'  \item{\code{term_distance}: term distance; optional, it is only appended when "details" is true}
#'  \item{\code{members}: members (represented as domain IDs) in overlaps; optional, it is only appended when "details" is true}
#' }
#' write(x) also returns the same data frame as view(x), in addition to a specified local file.
#' @docType methods
#' @keywords S4 methods
#' @name Eoutput-method
#' @rdname Eoutput-method
#' @seealso \code{\link{Eoutput-class}}

#' @rdname Eoutput-method
#' @aliases show,Eoutput-method
#' @export
setMethod("show", "Eoutput",
    function(object) {
        cat(sprintf("An object of S4 class '%s', containing following slots:", class(object)), "\n", sep="")
        cat(sprintf("  @domain: '%s'", object@domain), "\n", sep="")
        cat(sprintf("  @ontology: '%s'", object@ontology), "\n", sep="")
        cat(sprintf("  @term_info: a data.frame of %d terms X %d information", dim(object@term_info)[1],dim(object@term_info)[2]), "\n", sep="")
        cat(sprintf("  @anno: a list of %d terms, each storing annotated domains", length(object@anno)), "\n", sep="")
        cat(sprintf("  @data: a vector containing a group of %d input domains (annotatable)", length(object@data)), "\n", sep="")
        cat(sprintf("  @background: a vector containing a group of %d background domains (annotatable)", length(object@background)), "\n", sep="")
        cat(sprintf("  @overlap: a list of %d terms, each containing domains overlapped with input domains", length(object@overlap)), "\n", sep="")
        cat(sprintf("  @zscore: a vector of %d terms, containing z-scores", length(object@zscore)), "\n", sep="")
        cat(sprintf("  @pvalue: a vector of %d terms, containing p-values", length(object@pvalue)), "\n", sep="")
        cat(sprintf("  @adjp: a vector of %d terms, containing adjusted p-values", length(object@adjp)), "\n", sep="")
        cat(sprintf("In summary, a total of %d terms ('%s') are analysed for a group of %d input domains ('%s')", dim(object@term_info)[1],object@ontology,length(object@data),object@domain), "\n", sep="")
    }
)

setGeneric("zscore", function(x) standardGeneric("zscore"))
#' @rdname Eoutput-method
#' @aliases zscore
#' @export
setMethod("zscore", "Eoutput", function(x) x@zscore)

setGeneric("pvalue", function(x) standardGeneric("pvalue"))
#' @rdname Eoutput-method
#' @aliases pvalue
#' @export
setMethod("pvalue", "Eoutput", function(x) x@pvalue)

setGeneric("adjp", function(x) standardGeneric("adjp"))
#' @rdname Eoutput-method
#' @aliases adjp
#' @export
setMethod("adjp", "Eoutput", function(x) x@adjp)

setGeneric("view", function(x, ...) standardGeneric("view"))
#' @rdname Eoutput-method
#' @aliases view
#' @export
setMethod("view", "Eoutput", 
    function(x, top_num=5, sortBy=c("pvalue","adjp","zscore","nAnno","nOverlap","none"), decreasing=NULL, details=T){
        sortBy <- match.arg(sortBy)
    
        if( is.null(top_num) || is.na(top_num) ){
            top_num <- length(x@term_info$term_id)
        }
        if ( top_num > length(x@term_info$term_id) ){
            top_num <- length(x@term_info$term_id)
        }
    
        tab <- data.frame( term_id          = x@term_info$term_id,
                           nAnno            = sapply(x@anno,length),
                           nGroup           = length(x@data),
                           nOverlap         = sapply(x@overlap,length),
                           zscore           = x@zscore,
                           pvalue           = x@pvalue,
                           adjp             = x@adjp,
                           term_name        = x@term_info$term_name,
                           term_namespace   = x@term_info$term_namespace,
                           term_distance    = x@term_info$term_distance,
                           members          = sapply(x@overlap, function(x) paste(x,collapse=','))
                          )
        # members          = ifelse(dim(x@term_info)[1]==1, paste(x@overlap,collapse=','), sapply(x@overlap, function(x) paste(x,collapse=',')))
    
        if(details == T){
            res <- tab[,c(1:11)]
        }else{
            res <- tab[,c(1:8)]
        }
    
        if(is.null(decreasing)){
            if(sortBy=="zscore" | sortBy=="nAnno" | sortBy=="nOverlap"){
                decreasing <- T
            }else{
                decreasing <- F
            }
        }
    
        switch(sortBy, 
            adjp={res <- res[order(res[,7], decreasing=decreasing)[1:top_num],]},
            pvalue={res <- res[order(res[,6], decreasing=decreasing)[1:top_num],]},
            zscore={res <- res[order(res[,5], decreasing=decreasing)[1:top_num],]},
            nAnno={res <- res[order(res[,2], decreasing=decreasing)[1:top_num],]},
            nOverlap={res <- res[order(res[,4], decreasing=decreasing)[1:top_num],]},
            none={res <- res[order(res[,1], decreasing=decreasing)[1:top_num],]}
        )
        
        res
    }
)

setGeneric("write", function(x, ...) standardGeneric("write"))
#' @rdname Eoutput-method
#' @aliases write
#' @export
setMethod("write", "Eoutput", 
    function(x, file="Eoutput.txt", verbose=T){
        if(file=='' || is.na(file) || is.null(file)){
            file <- "Eoutput.txt"
        }
        
        out <- view(x, top_num=NULL, sortBy="pvalue", details=TRUE)
        
        utils::write.table(out, file=file, col.names=T, row.names=F, sep="\t")
        
        if(verbose){
            message(sprintf("A file ('%s') has been written into your local directory ('%s')", file, getwd()), appendLF=T)
        }
        
        invisible(out)
    }
)


################################################################################
################################################################################
#' @title Definition for VIRTUAL S4 class AdjData
#' @description \code{AdjData} is union of other classes: either matrix or dgCMatrix (a sparse matrix in the package Matrix). It is used as a virtual class
#' @return Class AdjData
#' @import Matrix
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name AdjData-class
#' @rdname AdjData-class
#' @seealso \code{\link{Onto-class}}

#' @rdname AdjData-class
#' @aliases AdjData
#' @exportClass AdjData
setClassUnion("AdjData", c("matrix", "dgCMatrix"))

################################################################################
################################################################################
#' @title Definition for S4 class Onto
#' @description \code{Onto} has 2 slots: nodeInfo and adjMatrix
#' @return Class Onto
#' @slot nodeInfo An object of S4 class \code{\link{InfoDataFrame}}, describing information on nodes/terms.
#' @slot adjMatrix An object of S4 class \code{\link{AdjData}}, containing adjacency data matrix (for a direct graph), with rows for parent (arrow-outbound) and columns for children (arrow-inbound)
#' @section Creation:
#' An object of this class can be created via: \code{new("Onto", nodeInfo, adjMatrix)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{dim()}: }{retrieve the dimension in the object}
#' \item{\code{adjMatrix()}: }{retrieve the slot 'adjMatrix' in the object}
#' \item{\code{nodeInfo()}: }{retrieve the slot 'nodeInfo' (as class InfoDataFrame) in the object}
#' \item{\code{nInfo()}: }{retrieve nodeInfo (as data.frame) in the object}
#' \item{\code{nodeNames()}: }{retrieve node/term names (ie, row names of nodeInfo) in the object}
#' \item{\code{term_id()}: }{retrieve term id (ie, column 'term_id' of nodeInfo) in the object, if any}
#' \item{\code{term_name()}: }{retrieve term id (ie, column 'term_name' of nodeInfo) in the object, if any}
#' \item{\code{term_namespace()}: }{retrieve term id (ie, column 'term_namespace' of nodeInfo) in the object, if any}
#' \item{\code{term_distance()}: }{retrieve term id (ie, column 'term_distance' of nodeInfo) in the object, if any}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' \item{\code{as(matrix, "Onto")}: }{convert a matrix to an object of class Onto}
#' \item{\code{as(dgCMatrix, "Onto")}: }{convert a dgCMatrix (a sparse matrix) to an object of class Onto}
#' \item{\code{[i]}: }{get the subset of the same class}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("Onto")}: }{show the class definition}
#' \item{\code{showMethods(classes="Onto")}: }{show the method definition upon this class}
#' \item{\code{getSlots("Onto")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("Onto")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="Onto")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name Onto-class
#' @seealso \code{\link{Onto-method}}
#' @examples
#' # create an object of class Onto, only given a matrix
#' adjM <- matrix(runif(25),nrow=5,ncol=5)
#' as(adjM, "Onto")
#'
#' # create an object of class Onto, given a matrix plus information on nodes
#' # 1) create nodeI: an object of class InfoDataFrame
#' data <- data.frame(term_id=paste("Term", 1:5, sep="_"), term_name=I(LETTERS[1:5]), term_namespace=rep("Namespace",5), term_distance=1:5, row.names=paste("Term", 1:5, sep="_"))
#' nodeI <- new("InfoDataFrame", data=data)
#' nodeI
#' # 2) create an object of class Onto
#' # VERY IMPORTANT: make sure having consistent names between nodeInfo and adjMatrix
#' adjM <- matrix(runif(25),nrow=5,ncol=5)
#' colnames(adjM) <- rownames(adjM) <- rowNames(nodeI)
#' x <- new("Onto", adjMatrix=adjM, nodeInfo=nodeI)
#' x
#' # 3) look at various methods defined on class Onto
#' dim(x)
#' adjMatrix(x)
#' nodeInfo(x)
#' nInfo(x)
#' nodeNames(x)
#' term_id(x)
#' term_namespace(x)
#' term_distance(x)
#' # 4) get the subset
#' x[1:2]

#' @rdname Onto-class
#' @aliases Onto
#' @exportClass Onto
setClass(
    Class="Onto",
    representation(
        nodeInfo = "InfoDataFrame",
        adjMatrix = "AdjData"
    ),
    prototype = prototype(
        nodeInfo = new("InfoDataFrame",dimLabels=c("termNames", "termColumns")),
        adjMatrix = matrix()
    ),
    validity = function(object){
        msg <- NULL
        # dimension for adjMatrix
        adim <- dim(object)
        if(adim[1]!=adim[2]){
            msg <- append(msg, "dimensions differ for adjacent matrix")
        }
        if( dim(nodeInfo(object))[1] != 0 ){
            if (adim[1] != dim(nodeInfo(object))[1]){
                msg <- append(msg, "term numbers differ between nodeInfo and adjMatrix")
            }
            if (!identical(nodeNames(object), rownames(adjMatrix(object)))){
                msg <- append(msg, "term names differ between nodeInfo and adjMatrix")
            }
        }
        if (is.null(msg)) TRUE else msg
    }
)

########################################
#' @title Methods defined for S4 class Onto
#' @description Methods defined for class \code{Onto}.
#' @param x an object of class \code{Onto}
#' @param object an object of class \code{Onto}
#' @param i an index
#' @param j an index
#' @param drop a logic for matrices and arrays. If TRUE the result is coerced to the lowest possible dimension. This only works for extracting elements, not for the replacement
#' @param ... additional parameters
#' @docType methods
#' @keywords S4 methods
#' @name Onto-method
#' @rdname Onto-method
#' @seealso \code{\link{Onto-class}}

#' @rdname Onto-method
#' @aliases dim,Onto-method
#' @export
setMethod("dim", "Onto", function(x) dim(x@adjMatrix))

setGeneric("adjMatrix", function(x) standardGeneric("adjMatrix"))
#' @rdname Onto-method
#' @aliases adjMatrix
#' @export
setMethod("adjMatrix", "Onto", function(x) x@adjMatrix)

setGeneric("nodeInfo", function(x) standardGeneric("nodeInfo"))
#' @rdname Onto-method
#' @aliases nodeInfo
#' @export
setMethod("nodeInfo", "Onto", function(x) x@nodeInfo)

setGeneric("nInfo", function(object) standardGeneric("nInfo"))
#' @rdname Onto-method
#' @aliases nInfo
#' @export
setMethod("nInfo", signature(object="Onto"), function(object){
    data <- Data(nodeInfo(object))
    if(sum(dim(data))==0){
        cat("No data is available\n", sep="")
    }else{
        data
    }
})

setGeneric("nodeNames", function(object) standardGeneric("nodeNames"))
#' @rdname Onto-method
#' @aliases nodeNames
#' @export
setMethod("nodeNames", signature(object="Onto"), function(object) rowNames(nodeInfo(object)))

setGeneric("term_id", function(object) standardGeneric("term_id"))
#' @rdname Onto-method
#' @aliases term_id
#' @export
setMethod("term_id", signature(object="Onto"),
    function(object){
        if(is.null(nInfo(object)$term_id)){
            stop(sprintf("This method '%s' cannot be used for this object '%s' of S4 class '%s'", "term_id()", deparse(substitute(object)), class(object)))
        }else{
            as.vector(nInfo(object)$term_id)
        }
    }
)

setGeneric("term_name", function(object) standardGeneric("term_name"))
#' @rdname Onto-method
#' @aliases term_name
#' @export
setMethod("term_name", signature(object="Onto"),
    function(object){
        if(is.null(nInfo(object)$term_name)){
            stop(sprintf("This method '%s' cannot be used for this object '%s' of S4 class '%s'", "term_name()", deparse(substitute(object)), class(object)))
        }else{
            as.vector(nInfo(object)$term_name)
        }
    }
)

setGeneric("term_namespace", function(object) standardGeneric("term_namespace"))
#' @rdname Onto-method
#' @aliases term_namespace
#' @export
setMethod("term_namespace", signature(object="Onto"),
    function(object){
        if(is.null(nInfo(object)$term_namespace)){
            stop(sprintf("This method '%s' cannot be used for this object '%s' of S4 class '%s'", "term_namespace()", deparse(substitute(object)), class(object)))
        }else{
            as.vector(nInfo(object)$term_namespace)
        }
    }
)

setGeneric("term_distance", function(object) standardGeneric("term_distance"))
#' @rdname Onto-method
#' @aliases term_distance
#' @export
setMethod("term_distance", signature(object="Onto"),
    function(object){
        if(is.null(nInfo(object)$term_distance)){
            stop(sprintf("This method '%s' cannot be used for this object '%s' of S4 class '%s'", "term_distance()", deparse(substitute(object)), class(object)))
        }else{
            as.vector(nInfo(object)$term_distance)
        }
    }
)

#' @rdname Onto-method
#' @name matrix2Onto
setAs("matrix", "Onto", function(from) {
    ## for nodeInfo    
    rn <- rownames(from)
    if(is.null(rn)){
        rn <- 1:nrow(from)
        rownames(from) <- rn
    }
    nodeI <- new("InfoDataFrame", data=data.frame(term_id=rn, row.names=rn))
    ## for Onto
    new("Onto", adjMatrix=from, nodeInfo=nodeI)
})


#' @rdname Onto-method
#' @name dgCMatrix2Onto
setAs("dgCMatrix", "Onto", function(from) {
    ## for nodeInfo    
    rn <- rownames(from)
    if(is.null(rn)){
        rn <- 1:nrow(from)
        rownames(from) <- rn
    }
    nodeI <- new("InfoDataFrame", data=data.frame(term_id=rn, row.names=rn))
    ## for Onto
    new("Onto", adjMatrix=from, nodeInfo=nodeI)
})


#' @rdname Onto-method
#' @export
setMethod("show", 
    signature=signature(object="Onto"),
    function(object) {
        cat("An object of S4 class '", class(object), "'\n", sep="")
        adim <- dim(object)
        if (length(adim)>1){
            cat("@adjMatrix:", if (length(adim)>1) paste("a direct matrix of",adim[[1]], "terms (parents/from) X",adim[[2]], "terms (children/to)") else NULL, "\n")
        }
        ## nodeInfo
        if( dim(nodeInfo(object))[1] != 0 ){
            cat("@nodeInfo (", class(nodeInfo(object)), ")\n", sep="")
            .showInfoDataFrame(
                nodeInfo(object),
                labels=list(
                    object="nodeInfo",
                    termNames="nodeNames",
                    varLabels="nodeAttr"
                )
            )
        }else{
            cat("nodeInfo (NULL)\n", sep="")
        }
    }
)

#' @rdname Onto-method
#' @aliases [,Onto-method
#' @export
setMethod("[", signature(x="Onto"), 
    function(x, i, j, ..., drop = FALSE) {
        if (missing(drop)){
            drop <- FALSE
        }
        if (missing(i) && missing(j)) {
            if (length(list(...))!=0){
                stop("specify terms to subset")
            }
            return(x)
        }
        
        if (!missing(i)) {
            nD <- nodeInfo(x)[i,,..., drop=drop]
        }else{
            nD <- nodeInfo(x)
        }
        
        if (!missing(i)){
            aD <- adjMatrix(x)[i,i]
        }else{
            aD <- adjMatrix(x)
        }
        
        x <- new("Onto", adjMatrix=aD, nodeInfo=nD)
    }
)


################################################################################
################################################################################
#' @title Definition for S4 class Dnetwork
#' @description \code{Dnetwork} is an S4 class to store a domain network, such as the one from semantic similairty between pairs of domains by \code{\link{dcDAGdomainSim}}. It has 2 slots: nodeInfo and adjMatrix
#' @return Class Dnetwork
#' @slot nodeInfo An object of S4 class \code{\link{InfoDataFrame}}, describing information on nodes/domains.
#' @slot adjMatrix An object of S4 class \code{\link{AdjData}}, containing symmetric adjacency data matrix for an indirect domain network
#' @section Creation:
#' An object of this class can be created via: \code{new("Dnetwork", nodeInfo, adjMatrix)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{dim()}: }{retrieve the dimension in the object}
#' \item{\code{adjMatrix()}: }{retrieve the slot 'adjMatrix' in the object}
#' \item{\code{nodeInfo()}: }{retrieve the slot 'nodeInfo' (as class InfoDataFrame) in the object}
#' \item{\code{nInfo()}: }{retrieve nodeInfo (as data.frame) in the object}
#' \item{\code{nodeNames()}: }{retrieve node/term names (ie, row names of nodeInfo) in the object}
#' \item{\code{id()}: }{retrieve domain id (ie, column 'id' of nodeInfo) in the object, if any}
#' \item{\code{level()}: }{retrieve domain level (ie, column 'level' of nodeInfo) in the object, if any}
#' \item{\code{description()}: }{retrieve domain description (ie, column 'description' of nodeInfo) in the object, if any}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' \item{\code{as(matrix, "Dnetwork")}: }{convert a matrix to an object of class Dnetwork}
#' \item{\code{as(dgCMatrix, "Dnetwork")}: }{convert a dgCMatrix (a sparse matrix) to an object of class Dnetwork}
#' \item{\code{[i]}: }{get the subset of the same class}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("Dnetwork")}: }{show the class definition}
#' \item{\code{showMethods(classes="Dnetwork")}: }{show the method definition upon this class}
#' \item{\code{getSlots("Dnetwork")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("Dnetwork")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="Dnetwork")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name Dnetwork-class
#' @seealso \code{\link{Dnetwork-method}}
#' @examples
#' # create an object of class Dnetwork, only given a matrix
#' adjM <- matrix(runif(25),nrow=5,ncol=5)
#' as(adjM, "Dnetwork")
#'
#' # create an object of class Dnetwork, given a matrix plus information on nodes
#' # 1) create nodeI: an object of class InfoDataFrame
#' data <- data.frame(id=paste("Domain", 1:5, sep="_"), level=rep("SCOP",5), description=I(LETTERS[1:5]), row.names=paste("Domain", 1:5, sep="_"))
#' nodeI <- new("InfoDataFrame", data=data)
#' nodeI
#' # 2) create an object of class Dnetwork
#' # VERY IMPORTANT: make sure having consistent names between nodeInfo and adjMatrix
#' adjM <- matrix(runif(25),nrow=5,ncol=5)
#' colnames(adjM) <- rownames(adjM) <- rowNames(nodeI)
#' x <- new("Dnetwork", adjMatrix=adjM, nodeInfo=nodeI)
#' x
#' # 3) look at various methods defined on class Dnetwork
#' dim(x)
#' adjMatrix(x)
#' nodeInfo(x)
#' nInfo(x)
#' nodeNames(x)
#' id(x)
#' level(x)
#' description(x)
#' # 4) get the subset
#' x[1:2]

#' @rdname Dnetwork-class
#' @aliases Dnetwork
#' @exportClass Dnetwork
setClass(
    Class="Dnetwork",
    representation(
        nodeInfo = "InfoDataFrame",
        adjMatrix = "AdjData"
    ),
    prototype = prototype(
        nodeInfo = new("InfoDataFrame",dimLabels=c("domainNames", "domainColumns")),
        adjMatrix = matrix()
    ),
    validity = function(object){
        msg <- NULL
        # dimension for adjMatrix
        adim <- dim(object)
        if(adim[1]!=adim[2]){
            msg <- append(msg, "dimensions differ for adjacent matrix")
        }
        if( dim(nodeInfo(object))[1] != 0 ){
            if (adim[1] != dim(nodeInfo(object))[1]){
                msg <- append(msg, "domain numbers differ between nodeInfo and adjMatrix")
            }
            if (!identical(nodeNames(object), rownames(adjMatrix(object)))){
                msg <- append(msg, "domain names differ between nodeInfo and adjMatrix")
            }
        }
        if (is.null(msg)) TRUE else msg
    }
)

########################################
#' @title Methods defined for S4 class Dnetwork
#' @description Methods defined for class \code{Dnetwork}.
#' @param x an object of class \code{Dnetwork}
#' @param object an object of class \code{Dnetwork}
#' @param i an index
#' @param j an index
#' @param drop a logic for matrices and arrays. If TRUE the result is coerced to the lowest possible dimension. This only works for extracting elements, not for the replacement
#' @param ... additional parameters
#' @docType methods
#' @keywords S4 methods
#' @name Dnetwork-method
#' @rdname Dnetwork-method
#' @seealso \code{\link{Dnetwork-class}}

#' @rdname Dnetwork-method
#' @aliases dim,Dnetwork-method
#' @export
setMethod("dim", "Dnetwork", function(x) dim(x@adjMatrix))

#' @rdname Dnetwork-method
#' @aliases adjMatrix,Dnetwork-method
#' @export
setMethod("adjMatrix", "Dnetwork", function(x) x@adjMatrix)

#' @rdname Dnetwork-method
#' @aliases nodeInfo,Dnetwork-method
#' @export
setMethod("nodeInfo", "Dnetwork", function(x) x@nodeInfo)

#' @rdname Dnetwork-method
#' @aliases nInfo,Dnetwork-method
#' @export
setMethod("nInfo", signature(object="Dnetwork"), function(object){
    data <- Data(nodeInfo(object))
    if(sum(dim(data))==0){
        cat("No data is available\n", sep="")
    }else{
        data
    }
})

#' @rdname Dnetwork-method
#' @aliases nodeNames,Dnetwork-method
#' @export
setMethod("nodeNames", signature(object="Dnetwork"), function(object) rowNames(nodeInfo(object)))

setGeneric("id", function(object) standardGeneric("id"))
#' @rdname Dnetwork-method
#' @aliases id
#' @export
setMethod("id", signature(object="Dnetwork"),
    function(object){
        if(is.null(nInfo(object)$id)){
            stop(sprintf("This method '%s' cannot be used for this object '%s' of S4 class '%s'", "id()", deparse(substitute(object)), class(object)))
        }else{
            as.vector(nInfo(object)$id)
        }
    }
)

setGeneric("level", function(object) standardGeneric("level"))
#' @rdname Dnetwork-method
#' @aliases level
#' @export
setMethod("level", signature(object="Dnetwork"),
    function(object){
        if(is.null(nInfo(object)$level)){
            stop(sprintf("This method '%s' cannot be used for this object '%s' of S4 class '%s'", "level()", deparse(substitute(object)), class(object)))
        }else{
            as.vector(nInfo(object)$level)
        }
    }
)

setGeneric("description", function(object) standardGeneric("description"))
#' @rdname Dnetwork-method
#' @aliases description
#' @export
setMethod("description", signature(object="Dnetwork"),
    function(object){
        if(is.null(nInfo(object)$description)){
            stop(sprintf("This method '%s' cannot be used for this object '%s' of S4 class '%s'", "description()", deparse(substitute(object)), class(object)))
        }else{
            as.vector(nInfo(object)$description)
        }
    }
)

#' @rdname Dnetwork-method
#' @name matrix2Dnetwork
setAs("matrix", "Dnetwork", function(from) {
    ## for nodeInfo
    rn <- rownames(from)
    if(is.null(rn)){
        rn <- 1:nrow(from)
        rownames(from) <- rn
    }
    nodeI <- new("InfoDataFrame", data=data.frame(id=rn, row.names=rn))
    ## for Dnetwork
    new("Dnetwork", adjMatrix=from, nodeInfo=nodeI)
})

#' @rdname Dnetwork-method
#' @name dgCMatrix2Dnetwork
setAs("dgCMatrix", "Dnetwork", function(from) {
    ## for nodeInfo    
    rn <- rownames(from)
    if(is.null(rn)){
        rn <- 1:nrow(from)
        rownames(from) <- rn
    }
    nodeI <- new("InfoDataFrame", data=data.frame(id=rn, row.names=rn))
    ## for Dnetwork
    new("Dnetwork", adjMatrix=from, nodeInfo=nodeI)
})


#' @rdname Dnetwork-method
#' @export
setMethod("show", 
    signature=signature(object="Dnetwork"),
    function(object) {
        cat("An object of S4 class '", class(object), "'\n", sep="")
        adim <- dim(object)
        if (length(adim)>1){
            cat("@adjMatrix:", if (length(adim)>1) paste("a weighted symmetric matrix of", adim[[1]], "domains X",adim[[2]], "domains") else NULL, "\n")
        }
        ## nodeInfo
        if( dim(nodeInfo(object))[1] != 0 ){
            cat("@nodeInfo (", class(nodeInfo(object)), ")\n", sep="")
            .showInfoDataFrame(
                nodeInfo(object),
                labels=list(
                    object="nodeInfo",
                    termNames="nodeNames",
                    varLabels="nodeAttr"
                )
            )
        }else{
            cat("nodeInfo (NULL)\n", sep="")
        }
    }
)

#' @rdname Dnetwork-method
#' @aliases [,Dnetwork-method
#' @export
setMethod("[", signature(x="Dnetwork"), 
    function(x, i, j, ..., drop = FALSE) {
        if (missing(drop)){
            drop <- FALSE
        }
        if (missing(i) && missing(j)) {
            if (length(list(...))!=0){
                stop("specify domains to subset")
            }
            return(x)
        }
        
        if (!missing(i)) {
            nD <- nodeInfo(x)[i,,..., drop=drop]
        }else{
            nD <- nodeInfo(x)
        }
        
        if (!missing(i)){
            aD <- adjMatrix(x)[i,i]
        }else{
            aD <- adjMatrix(x)
        }
        
        x <- new("Dnetwork", adjMatrix=aD, nodeInfo=nD)
    }
)


################################################################################
################################################################################
#' @title Definition for S4 class Cnetwork
#' @description \code{Cnetwork} is an S4 class to store a contact network, such as the one from RWR-based contact between samples/terms by \code{\link{dcRWRpipeline}}. It has 2 slots: nodeInfo and adjMatrix
#' @return Class Cnetwork
#' @slot nodeInfo An object of S4 class \code{\link{InfoDataFrame}}, describing information on nodes/domains.
#' @slot adjMatrix An object of S4 class \code{\link{AdjData}}, containing symmetric adjacency data matrix for an indirect domain network
#' @section Creation:
#' An object of this class can be created via: \code{new("Cnetwork", nodeInfo, adjMatrix)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{dim()}: }{retrieve the dimension in the object}
#' \item{\code{adjMatrix()}: }{retrieve the slot 'adjMatrix' in the object}
#' \item{\code{nodeInfo()}: }{retrieve the slot 'nodeInfo' (as class InfoDataFrame) in the object}
#' \item{\code{nInfo()}: }{retrieve nodeInfo (as data.frame) in the object}
#' \item{\code{nodeNames()}: }{retrieve node/term names (ie, row names of nodeInfo) in the object}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' \item{\code{as(matrix, "Cnetwork")}: }{convert a matrix to an object of class Cnetwork}
#' \item{\code{as(dgCMatrix, "Cnetwork")}: }{convert a dgCMatrix (a sparse matrix) to an object of class Cnetwork}
#' \item{\code{[i]}: }{get the subset of the same class}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("Cnetwork")}: }{show the class definition}
#' \item{\code{showMethods(classes="Cnetwork")}: }{show the method definition upon this class}
#' \item{\code{getSlots("Cnetwork")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("Cnetwork")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="Cnetwork")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name Cnetwork-class
#' @seealso \code{\link{Cnetwork-method}}
#' @examples
#' # create an object of class Cnetwork, only given a matrix
#' adjM <- matrix(runif(25),nrow=5,ncol=5)
#' as(adjM, "Cnetwork")
#'
#' # create an object of class Cnetwork, given a matrix plus information on nodes
#' # 1) create nodeI: an object of class InfoDataFrame
#' data <- data.frame(id=paste("Domain", 1:5, sep="_"), level=rep("SCOP",5), description=I(LETTERS[1:5]), row.names=paste("Domain", 1:5, sep="_"))
#' nodeI <- new("InfoDataFrame", data=data)
#' nodeI
#' # 2) create an object of class Cnetwork
#' # VERY IMPORTANT: make sure having consistent names between nodeInfo and adjMatrix
#' adjM <- matrix(runif(25),nrow=5,ncol=5)
#' colnames(adjM) <- rownames(adjM) <- rowNames(nodeI)
#' x <- new("Cnetwork", adjMatrix=adjM, nodeInfo=nodeI)
#' x
#' # 3) look at various methods defined on class Cnetwork
#' dim(x)
#' adjMatrix(x)
#' nodeInfo(x)
#' nInfo(x)
#' nodeNames(x)
#' # 4) get the subset
#' x[1:2]

#' @rdname Cnetwork-class
#' @aliases Cnetwork
#' @exportClass Cnetwork
setClass(
    Class="Cnetwork",
    representation(
        nodeInfo = "InfoDataFrame",
        adjMatrix = "AdjData"
    ),
    prototype = prototype(
        nodeInfo = new("InfoDataFrame",dimLabels=c("sampleNames", "sampleColumns")),
        adjMatrix = matrix()
    ),
    validity = function(object){
        msg <- NULL
        # dimension for adjMatrix
        adim <- dim(object)
        if(adim[1]!=adim[2]){
            msg <- append(msg, "dimensions differ for adjacent matrix")
        }
        if( dim(nodeInfo(object))[1] != 0 ){
            if (adim[1] != dim(nodeInfo(object))[1]){
                msg <- append(msg, "sample numbers differ between nodeInfo and adjMatrix")
            }
            if (!identical(nodeNames(object), rownames(adjMatrix(object)))){
                msg <- append(msg, "sample names differ between nodeInfo and adjMatrix")
            }
        }
        if (is.null(msg)) TRUE else msg
    }
)

########################################
#' @title Methods defined for S4 class Cnetwork
#' @description Methods defined for class \code{Cnetwork}.
#' @param x an object of class \code{Cnetwork}
#' @param object an object of class \code{Cnetwork}
#' @param i an index
#' @param j an index
#' @param drop a logic for matrices and arrays. If TRUE the result is coerced to the lowest possible dimension. This only works for extracting elements, not for the replacement
#' @param ... additional parameters
#' @docType methods
#' @keywords S4 methods
#' @name Cnetwork-method
#' @rdname Cnetwork-method
#' @seealso \code{\link{Cnetwork-class}}

#' @rdname Cnetwork-method
#' @aliases dim,Cnetwork-method
#' @export
setMethod("dim", "Cnetwork", function(x) dim(x@adjMatrix))

#' @rdname Cnetwork-method
#' @aliases adjMatrix,Cnetwork-method
#' @export
setMethod("adjMatrix", "Cnetwork", function(x) x@adjMatrix)

#' @rdname Cnetwork-method
#' @aliases nodeInfo,Cnetwork-method
#' @export
setMethod("nodeInfo", "Cnetwork", function(x) x@nodeInfo)

#' @rdname Cnetwork-method
#' @aliases nInfo,Cnetwork-method
#' @export
setMethod("nInfo", signature(object="Cnetwork"), function(object){
    data <- Data(nodeInfo(object))
    if(sum(dim(data))==0){
        cat("No data is available\n", sep="")
    }else{
        data
    }
})

#' @rdname Cnetwork-method
#' @aliases nodeNames,Cnetwork-method
#' @export
setMethod("nodeNames", signature(object="Cnetwork"), function(object) rowNames(nodeInfo(object)))

#' @rdname Cnetwork-method
#' @name matrix2Cnetwork
setAs("matrix", "Cnetwork", function(from) {
    ## for nodeInfo
    rn <- rownames(from)
    if(is.null(rn)){
        rn <- 1:nrow(from)
        rownames(from) <- rn
    }
    nodeI <- new("InfoDataFrame", data=data.frame(id=rn, row.names=rn))
    ## for Cnetwork
    new("Cnetwork", adjMatrix=from, nodeInfo=nodeI)
})

#' @rdname Cnetwork-method
#' @name dgCMatrix2Cnetwork
setAs("dgCMatrix", "Cnetwork", function(from) {
    ## for nodeInfo    
    rn <- rownames(from)
    if(is.null(rn)){
        rn <- 1:nrow(from)
        rownames(from) <- rn
    }
    nodeI <- new("InfoDataFrame", data=data.frame(name=rn, row.names=rn))
    ## for Cnetwork
    new("Cnetwork", adjMatrix=from, nodeInfo=nodeI)
})


#' @rdname Cnetwork-method
#' @export
setMethod("show", 
    signature=signature(object="Cnetwork"),
    function(object) {
        cat("An object of S4 class '", class(object), "'\n", sep="")
        adim <- dim(object)
        if (length(adim)>1){
            cat("@adjMatrix:", if (length(adim)>1) paste("a weighted symmetric matrix of", adim[[1]], "samples/terms X",adim[[2]], "samples/terms") else NULL, "\n")
        }
        ## nodeInfo
        if( dim(nodeInfo(object))[1] != 0 ){
            cat("@nodeInfo (", class(nodeInfo(object)), ")\n", sep="")
            .showInfoDataFrame(
                nodeInfo(object),
                labels=list(
                    object="nodeInfo",
                    termNames="nodeNames",
                    varLabels="nodeAttr"
                )
            )
        }else{
            cat("nodeInfo (NULL)\n", sep="")
        }
    }
)

#' @rdname Cnetwork-method
#' @aliases [,Cnetwork-method
#' @export
setMethod("[", signature(x="Cnetwork"), 
    function(x, i, j, ..., drop = FALSE) {
        if (missing(drop)){
            drop <- FALSE
        }
        if (missing(i) && missing(j)) {
            if (length(list(...))!=0){
                stop("specify samples/terms to subset")
            }
            return(x)
        }
        
        if (!missing(i)) {
            nD <- nodeInfo(x)[i,,..., drop=drop]
        }else{
            nD <- nodeInfo(x)
        }
        
        if (!missing(i)){
            aD <- adjMatrix(x)[i,i]
        }else{
            aD <- adjMatrix(x)
        }
        
        x <- new("Cnetwork", adjMatrix=aD, nodeInfo=nD)
    }
)


################################################################################
################################################################################
#' @title Definition for S4 class Coutput
#' @description \code{Coutput} is an S4 class to store output by \code{\link{dcRWRpipeline}}.
#' @return Class Coutput
#' @slot ratio A symmetrix matrix, containing ratio
#' @slot zscore A symmetrix matrix, containing z-scores
#' @slot pvalue A symmetrix matrix, containing p-values
#' @slot adjp A symmetrix matrix, containing adjusted p-values
#' @slot cnetwork An object of S4 class \code{\link{Cnetwork}}, storing contact network.
#' @section Creation:
#' An object of this class can be created via: \code{new("Coutput", ratio, zscore, pvalue, adjp, cnetwork)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{ratio()}: }{retrieve the slot 'ratio' in the object}
#' \item{\code{zscore()}: }{retrieve the slot 'zscore' in the object}
#' \item{\code{pvalue()}: }{retrieve the slot 'pvalue' in the object}
#' \item{\code{adjp()}: }{retrieve the slot 'adjp' in the object}
#' \item{\code{cnetwork()}: }{retrieve the slot 'cnetwork' in the object}
#' \item{\code{write()}: }{write the object into a local file}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("Coutput")}: }{show the class definition}
#' \item{\code{showMethods(classes="Coutput")}: }{show the method definition upon this class}
#' \item{\code{getSlots("Coutput")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("Coutput")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="Coutput")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name Coutput-class
#' @rdname Coutput-class
#' @seealso \code{\link{Coutput-method}}
#' @examples
#' \dontrun{
#' # 1) load onto.GOMF (as 'Onto' object)
#' g <- dcRDataLoader('onto.GOMF')
#'
#' # 2) load SCOP superfamilies annotated by GOMF (as 'Anno' object)
#' Anno <- dcRDataLoader('SCOP.sf2GOMF')
#'
#' # 3) prepare for ontology appended with annotation information
#' dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths", verbose=TRUE)
#'
#' # 4) calculate pair-wise semantic similarity between 10 randomly chosen domains 
#' alldomains <- unique(unlist(nInfo(dag)$annotations))
#' domains <- sample(alldomains,10)
#' dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
#' dnetwork
#' 
#' # 5) estimate RWR dating based sample/term relationships
#' # define sets of seeds as data
#' # each seed with equal weight (i.e. all non-zero entries are '1')
#' data <- data.frame(aSeeds=c(1,0,1,0,1), bSeeds=c(0,0,1,0,1))
#' rownames(data) <- id(dnetwork)[1:5]
#' # calcualte their two contact graph
#' coutput <- dcRWRpipeline(data=data, g=dnetwork, parallel=FALSE)
#' coutput
#'
#' # 6) write into the file 'Coutput.txt' in your local directory
#' write(coutput, file='Coutput.txt', saveBy="adjp")
#'
#' # 7) retrieve several slots directly
#' ratio(coutput)
#' zscore(coutput)
#' pvalue(coutput)
#' adjp(coutput)
#' cnetwork(coutput)
#' }

#' @rdname Coutput-class
#' @aliases Coutput
#' @exportClass Coutput
setClass(
    Class="Coutput",
    representation(
        ratio    = "matrix",
        zscore   = "matrix",
        pvalue   = "matrix",
        adjp     = "matrix",
        cnetwork = "Cnetwork"
    ),
    prototype = prototype(
        ratio    = matrix(),
        zscore   = matrix(),
        pvalue   = matrix(),
        adjp     = matrix(),
        cnetwork = new("Cnetwork")
    ),
    validity = function(object){
        if(dim(object@ratio)[1]!=dim(object@zscore)[1]){
            return("Dimension is not matched")
        }else{
            return(TRUE)
        }
    }
)

########################################
#' @title Methods defined for S4 class Coutput
#' @description Methods defined for S4 class \code{Coutput}.
#' @param object an object of S4 class \code{Coutput}. Usually this is an output from \code{\link{dcRWRpipeline}}
#' @param x an object of S4 class \code{Coutput}. Usually this is part of the output from \code{\link{dcRWRpipeline}}
#' @param saveBy which statistics will be saved. It can be "pvalue" for p value, "adjp" for adjusted p value, "zscore" for z-score, "ratio" for ratio
#' @param file a character specifying a file name written into. By default, it is 'Coutput.txt'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' write(x) also returns a symmetrix matrix storing the specific statistics
#' @docType methods
#' @keywords S4 methods
#' @name Coutput-method
#' @rdname Coutput-method
#' @seealso \code{\link{Coutput-class}}

#' @rdname Coutput-method
#' @aliases show,Coutput-method
#' @export
setMethod("show", "Coutput",
    function(object) {
        cat(sprintf("An object of S4 class '%s', containing following slots:", class(object)), "\n", sep="")
        cat(sprintf("  @ratio: a matrix of %d X %d, containing ratio", dim(object@ratio)[1], dim(object@ratio)[2]), "\n", sep="")
        cat(sprintf("  @zscore: a matrix of %d X %d, containing z-scores", dim(object@zscore)[1], dim(object@zscore)[2]), "\n", sep="")
        cat(sprintf("  @pvalue: a matrix of %d X %d, containing p-values", dim(object@pvalue)[1], dim(object@pvalue)[2]), "\n", sep="")
        cat(sprintf("  @adjp: a matrix of %d X %d, containing adjusted p-values", dim(object@adjp)[1], dim(object@adjp)[2]), "\n", sep="")
        cat(sprintf("  @cnetwork: an object of S4 class 'Cnetwork', containing %d interacting nodes", dim(object@cnetwork)[1]), "\n", sep="")
    }
)

setGeneric("ratio", function(x) standardGeneric("ratio"))
#' @rdname Coutput-method
#' @aliases ratio
#' @export
setMethod("ratio", "Coutput", function(x) x@ratio)

#' @rdname Coutput-method
#' @aliases zscore,Coutput-method
#' @export
setMethod("zscore", "Coutput", function(x) x@zscore)

#' @rdname Coutput-method
#' @aliases pvalue,Coutput-method
#' @export
setMethod("pvalue", "Coutput", function(x) x@pvalue)

#' @rdname Coutput-method
#' @aliases adjp,Coutput-method
#' @export
setMethod("adjp", "Coutput", function(x) x@adjp)

setGeneric("cnetwork", function(x) standardGeneric("cnetwork"))
#' @rdname Coutput-method
#' @aliases cnetwork
#' @export
setMethod("cnetwork", "Coutput", function(x) x@cnetwork)

#' @rdname Coutput-method
#' @aliases write,Coutput-method
#' @export
setMethod("write", "Coutput", 
    function(x, file="Coutput.txt", saveBy=c("adjp","pvalue","zscore","ratio"), verbose=T){
        if(file=='' || is.na(file) || is.null(file)){
            file <- "Coutput.txt"
        }
        
        saveBy <- match.arg(saveBy)
        switch(saveBy, 
            adjp={res <- x@adjp},
            pvalue={res <- x@pvalue},
            zscore={res <- x@zscore},
            ratio={res <- x@ratio}
        )
        
        res <- cbind(rownames(res), res)
        utils::write.table(res, file=file, col.names=T, row.names=F, sep="\t")
        
        if(verbose){
            message(sprintf("A file ('%s') has been written into your local directory ('%s')", file, getwd()), appendLF=T)
        }
        
        invisible(res)
    }
)

