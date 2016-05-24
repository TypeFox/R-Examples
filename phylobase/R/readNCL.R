### This file contains the source code for the functions:
###  - readNCL (generic function)
###  - readNexus (wrapper for readNCL importing Nexus files)
###  - readNewick (wrapper for readNCL importing Newick files)

##' Create a \code{phylo4}, \code{phylo4d} or \code{data.frame} object
##' from a NEXUS or a Newick file
##'
##' \code{readNexus} reads a NEXUS file and outputs a \code{phylo4},
##' \code{phylo4d} or \code{data.frame} object.
##'
##' \code{readNewick} reads a Newick file and outputs a \code{phylo4}
##' or \code{phylo4d} object.
##'
##' \code{readNexus} is used internally by both \code{readNexus} and
##' \code{readNewick} to extract data held in a tree files,
##' specifically in NEXUS files from DATA, CHARACTER or TREES
##' blocks.
##'
##' The \code{type} argument specifies which of these is returned:
##'
##' \describe{
##'
##'   \item{data}{will only return a \code{data.frame} of the contents
##' of all DATA and CHARACTER blocks.}
##'
##'   \item{tree}{will only return a \code{phylo4} object of the
##' contents of the TREES block.}
##'
##'   \item{all}{if only data or a tree are present in the file, this
##' option will act as the options above, returning either a
##' \code{data.frame} or a \code{phylo4} object respectively. If both
##' are present then a \code{phylo4d} object is returned containing
##' both.}
##'
##' }
##'
##' The function returns \code{NULL} if the \code{type} of
##' data requested is not present in the file, or if neither data nor
##' tree blocks are present.
##'
##' Depending on the context \code{readNexus} will call either the
##' \code{phylo4} or \code{phylo4d} constructor. The \code{phylo4d}
##' constructor will be used with \code{type="all"}, or if the option
##' \code{check.node.labels="asdata"} is invoked.
##'
##' \code{readNewick} imports Newick formatted tree files and will
##' return a \code{phylo4} or a \code{phylo4d} object if the option
##' \code{check.node.labels="asdata"} is invoked.
##'
##' For both \code{readNexus} and \code{readNewick}, the options for
##' \code{check.node.labels} can take the values:
##'
##' \describe{
##'
##' \item{keep}{the node labels of the trees will be passed as node
##' labels in the \code{phylo4} object}
##'
##' \item{drop}{the node labels of the trees will be ignored in the
##' \code{phylo4} object}
##'
##' \item{asdata}{the node labels will be passed as data and a
##' \code{phylo4d} object will be returned.}
##'
##' }
##'
##' If you use the option \code{asdata} on a file with no node labels,
##' a warning message is issued, and is thus equivalent to the value
##' \code{drop}.
##'
##' For both \code{readNexus} and \code{readNewick}, additional
##' arguments can be passed to the constructors such as \code{annote},
##' \code{missing.data} or \code{extra.data}. See the \sQuote{Details}
##' section of \code{\link{phylo4d-methods}} for the complete list of
##' options.
##'
##' @name Import Nexus and Newick files
##' @docType methods
##' @param file a NEXUS file for \code{readNexus} or a file that
##'     contains Newick formatted trees for \code{readNewick}.
##' @param simplify If TRUE, if there are multiple trees in the file,
##' only the first one is returned; otherwise a list of
##' \code{phylo4(d)} objects is returned if the file contains multiple
##' trees.
##' @param type Determines which type of objects to return, if present
##'     in the file (see Details).
##' @param spacesAsUnderscores In the NEXUS file format white spaces
##'     are not allowed in taxa labels and are represented by
##'     underscores. Therefore, NCL converts underscores found in taxa
##'     labels in the NEXUS file into white spaces
##'     (e.g. \code{species_1} will become \code{"species 1"}. If you
##'     want to preserve the underscores, set as TRUE, the default).
##' @param char.all If \code{TRUE}, returns all characters, even those
##'     excluded in the NEXUS file
##' @param polymorphic.convert If \code{TRUE}, converts polymorphic
##'     characters to missing data
##' @param levels.uniform If \code{TRUE}, uses the same levels for all
##'     characters
##' @param quiet If \code{FALSE} the output of the NCL interface is
##'     printed. This is mainly for debugging purposes. This option
##'     can considerably slow down the process if the tree is big or
##'     there are many trees in the file.
##' @param check.node.labels Determines how the node labels in the
##'     NEXUS or Newick files should be treated in the phylo4 object,
##'     see Details for more information.
##' @param return.labels Determines whether state names (if
##'     \code{TRUE}) or state codes should be returned.
##' @param file.format character indicating the format of the
##'     specified file (either \dQuote{\code{newick}} or
##'     \dQuote{\code{nexus}}). It's more convenient to just use
##'     \code{readNexus}  or \code{readNewick}.
##' @param check.names logical. If \sQuote{TRUE} then the names of the
##'     characters from the NEXUS file are checked to ensure that they
##'     are syntactically valid variable names and are not duplicated.
##'     If necessary they are adjusted using \sQuote{make.names}.
##' @param convert.edge.length logical. If \code{TRUE} negative edge
##'     lengths are replaced with 0. At this time \code{phylobase}
##'     does not accept objects with negative branch lengths, this
##'     workaround allows to import trees with negative branch
##'     lengths.
##' @param \dots Additional arguments to be passed to phylo4 or
##'     phylo4d constructor (see Details)
##' @return Depending on the value of \code{type} and the contents of
##'     the file, one of: a \code{data.frame}, a \linkS4class{phylo4}
##'     object, a \linkS4class{phylo4d} object or \code{NULL}.  If
##'     several trees are included in the NEXUS file and the option
##'     \code{simplify=FALSE} a list of \linkS4class{phylo4} or
##'     \linkS4class{phylo4d} objects is returned.
##' @note Underscores in state labels (i.e. trait or taxon names) will
##' be translated to spaces. Unless \code{check.names=FALSE}, trait
##' names will be converted to valid R names (see
##' \code{\link{make.names}}) on input to R, so spaces will be
##' translated to periods.
##' @author Brian O'Meara, Francois Michonneau, Derrick Zwickl
##' @seealso the \linkS4class{phylo4d} class, the \linkS4class{phylo4}
##'     class
##' @export
##' @rdname readNexus
##' @aliases readNCL
##' @keywords misc

readNCL <- function(file, simplify=FALSE, type=c("all", "tree","data"),
                    spacesAsUnderscores = TRUE, char.all=FALSE,
                    polymorphic.convert=TRUE, levels.uniform=FALSE, quiet=TRUE,
                    check.node.labels=c("keep", "drop", "asdata"), return.labels=TRUE,
                    file.format=c("nexus", "newick"), check.names=TRUE,
                    convert.edge.length=FALSE, ...) {


    type <- match.arg(type)
    file.format <- match.arg(file.format)

    check.node.labels <- match.arg(check.node.labels)

    if (type == "all" || type == "data") {
        returnData <- TRUE
    }
    else {
        returnData <- FALSE
    }
    if (type == "all" || type == "tree") {
        returnTrees <- TRUE
    }
    else {
        returnTrees <- FALSE
    }

    ## GetNCL returns a list containing:
    ##  $taxaNames: names of the taxa (from taxa block, implied or declared)
    ##  $treeNames: the names of the trees
    ##  $trees: a vector of (untranslated) Newick strings
    ##  $dataTypes: data type for each character block of the nexus file (length = number of chr blocks)
    ##  $nbCharacters: number of characters in each block (length = number of chr blocks)
    ##  $charLabels: the labels for the characters, i.e. the headers of the data frame to be returned
    ##    (length = number of chr blocks * sum of number of characters in each block)
    ##  $nbStates: the number of states of each character (equals 0 for non-standard types, length = number
    ##    of characters)
    ##  $stateLabels: the labels for the states of the characters, i.e. the levels of the factors to be returned
    ##  $dataChr: string that contains the data to be returned

    ncl <- rncl::rncl(file = file, file.format = file.format, spacesAsUnderscores = spacesAsUnderscores,
                      char.all = char.all, polymorphic.convert = polymorphic.convert,
                      levels.uniform = levels.uniform)


    ## Return Error message
    if (length(ncl) == 1 && names(ncl) == "ErrorMsg") {
        stop(ncl$ErrorMsg)
    }

    if (!quiet) message(ncl)

    ## Disclaimer
    if (!length(grep("\\{", ncl$dataChr)) && return.labels && !polymorphic.convert) {
        stop("At this stage, it's not possible to use the combination: ",
             "return.labels=TRUE and polymorphic.convert=FALSE for datasets ",
             "that contain polymorphic characters.")
    }

    if (returnData && length(ncl$dataChr)) {
        tipData <- vector("list", length(ncl$dataChr))
        for (iBlock in 1:length(ncl$dataTypes)) {
            chrCounter <- ifelse(iBlock == 1, 0, sum(ncl$nbCharacters[1:(iBlock-1)]))
            if (ncl$dataTypes[iBlock] == "Continuous") {
                for (iChar in 1:ncl$nbCharacters[iBlock]) {
                    i <- chrCounter + iChar
                    tipData[[i]] <- eval(parse(text=ncl$dataChr[i]))
                    names(tipData)[i] <- ncl$charLabels[i]
                }
            }
            else {
                if (ncl$dataTypes[iBlock] == "Standard") {
                    iForBlock <- integer(0)
                    for (iChar in 1:ncl$nbCharacters[iBlock]) {
                        i <- chrCounter + iChar
                        iForBlock <- c(iForBlock, i)
                        lblCounterMin <- ifelse(i == 1, 1, sum(ncl$nbStates[1:(i-1)]) + 1)
                        lblCounter <- seq(lblCounterMin, length.out=ncl$nbStates[i])
                        tipData[[i]] <- eval(parse(text=ncl$dataChr[i]))
                        names(tipData)[i] <- ncl$charLabels[i]
                        tipData[[i]] <- as.factor(tipData[[i]])

                        lbl <- ncl$stateLabels[lblCounter]
                        if (return.labels) {
                            if (any(nchar(gsub("\\s|_", "", lbl)) == 0)) {
                                warning("state labels are missing for \'", ncl$charLabels[i],
                                        "\', the option return.labels is thus ignored.")
                            }
                            else {
                                levels(tipData[[i]]) <- lbl
                            }
                        }
                    }
                    if (levels.uniform) {
                        allLevels <- character(0)
                        for (j in iForBlock) {
                            allLevels <- union(allLevels, levels(tipData[[j]]))
                        }
                        for (j in iForBlock) {
                            levels(tipData[[j]]) <- allLevels
                        }
                    }
                }
                else {
                    warning("This datatype is not currently supported by phylobase")
                    next
                    ## FIXME: different datatypes in a same file isn't going to work
                }
            }
        }
        tipData <- data.frame(tipData, check.names=check.names)
        if (length(ncl$taxaNames) == nrow(tipData)) {
            rownames(tipData) <- ncl$taxaNames
        }
        else stop("phylobase doesn't deal with multiple taxa block at this time.")
    }
    else {
        tipData <- NULL
    }

    if (returnTrees && length(ncl$trees) > 0) {
        listTrees <- vector("list", length(ncl$trees))

        for (i in 1:length(ncl$trees)) {

            isRooted <- is_rooted(ncl$parentVector[[i]])

            edgeMat <- get_edge_matrix(ncl$parentVector[[i]], isRooted)

            edgeLgth <- get_edge_length(ncl$branchLengthVector[[i]],
                                        ncl$parentVector[[i]],
                                        isRooted)

            tipLbl <- ncl$taxonLabelVector[[i]]

            if (convert.edge.length) {
                edgeLgth[edgeLgth < 0] <- 0
            }

            if (check.node.labels == "asdata" &&
                !has_node_labels(ncl$nodeLabelsVector[[i]])) {
                warning("Could not use value \"asdata\" for ",
                        "check.node.labels because there are no ",
                        "labels associated with the tree")
                check.node.labels <- "drop"
            }


            if (has_node_labels(ncl$nodeLabelsVector[[i]]) &&
                !identical(check.node.labels, "drop")) {
                nodeLbl <- ncl$nodeLabelsVector[[i]]
                rootNd <- attr(edgeMat, "root")
                nodeLbl[rootNd] <- nodeLbl[1]
                node_pos <- (length(tipLbl)+1):length(nodeLbl)
                nodeLbl <- nodeLbl[node_pos]

                if (identical(check.node.labels, "asdata")) {
                    tr <- phylo4(x = edgeMat,
                                 edge.length = edgeLgth,
                                 tip.label = tipLbl)
                    nodeDt <- label_to_data(nodeLbl, row.names = node_pos)
                    tr <- phylo4d(tr, node.data = nodeDt)
                } else {

                    tr <- phylo4(x = edgeMat,
                                 edge.length = edgeLgth,
                                 tip.label = tipLbl,
                                 node.label = nodeLbl)
                }
            } else {
                tr <- phylo4(x = edgeMat,
                             edge.length = edgeLgth,
                             tip.label = tipLbl)


            }

            listTrees[[i]] <- tr
            if (simplify) break
        }

        if (length(listTrees) == 1 || simplify)
            listTrees <- listTrees[[1]]

    } else {
        listTrees <- NULL
    }

###
    switch(type,
           "data" = {
               if (is.null(tipData)) {
                   toRet <- NULL
               }
               else {
                   toRet <- tipData
               }
           },
           "tree" = {
               if (is.null(listTrees)) {
                   toRet <- NULL
               }
               else {
                   toRet <- listTrees
               }
           },
           "all" = {
               if (is.null(tipData) && is.null(listTrees)) {
                   toRet <- NULL
               }
               else if (is.null(tipData)) {
                   toRet <- listTrees
               }
               else if (is.null(listTrees)) {
                   toRet <- tipData
               }
               else {
                   if (length(listTrees) > 1) {
                       toRet <- lapply(listTrees, function(tr)
                                       addData(tr, tip.data=tipData, ...))
                   }
                   else toRet <- addData(listTrees, tip.data=tipData, ...)
               }
           })
    toRet
}


## check if the implicit root is dichotomous
is_rooted <- function(parentVector) {
    tab_edg <- table(parentVector)
    if (tabulate(parentVector)[which(parentVector == 0)] > 2)
        FALSE
    else TRUE
}


## Returns the edge matrix from the parentVector (the i^th element is
## the descendant element of node i)
get_edge_matrix <- function(parentVector, isRooted) {
    edgeMat <- cbind(ancestor = parentVector,
                     descendant = 1:length(parentVector))
    rootNd <- edgeMat[which(edgeMat[, 1] == 0), 2]
    if (!isRooted) {
        edgeMat <- edgeMat[-which(edgeMat[, 1] == 0), ]
    }
    attr(edgeMat, "root") <- rootNd
    edgeMat
}

## Returns the edge lengths (missing are represented by -999)
get_edge_length <- function(branchLengthVector, parentVector, isRooted) {
    edgeLgth <- branchLengthVector
    if (isRooted) {
        edgeLgth[which(parentVector == 0)] <- NA
    } else {
        edgeLgth <- edgeLgth[which(parentVector != 0)]
    }
    edgeLgth[edgeLgth == -999] <- NA
    edgeLgth
}

## Tests whether there are node labels
has_node_labels <- function(nodeLabelsVector) {
    any(nzchar(nodeLabelsVector))
}

##' @rdname readNexus
##' @aliases readNexus
##' @export
readNexus <- function (file, simplify=FALSE, type=c("all", "tree", "data"),
                       char.all=FALSE, polymorphic.convert=TRUE,
                       levels.uniform=FALSE, quiet=TRUE,
                       check.node.labels=c("keep", "drop", "asdata"),
                       return.labels=TRUE, check.names=TRUE, convert.edge.length=FALSE,
                       ...) {

    return(readNCL(file=file, simplify=simplify, type=type, char.all=char.all,
                   polymorphic.convert=polymorphic.convert, levels.uniform=levels.uniform,
                   quiet=quiet, check.node.labels=check.node.labels,
                   return.labels=return.labels, file.format="nexus",
                   check.names=check.names, convert.edge.length=convert.edge.length, ...))
}

##' @rdname readNexus
##' @aliases readNewick
##' @export
readNewick <- function(file, simplify=FALSE, quiet=TRUE,
                       check.node.labels=c("keep", "drop", "asdata"),
                       convert.edge.length=FALSE, ...) {

    return(readNCL(file=file, simplify=simplify, quiet=quiet,
                   check.node.labels=check.node.labels, file.format="newick",
                   convert.edge.length=convert.edge.length, ...))
}
