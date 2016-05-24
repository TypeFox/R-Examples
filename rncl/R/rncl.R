##' Returns a list of the elements contained in a NEXUS file used to
##' build phylogenetic objects in R
##'
##' NEXUS is a common file format used in phylogenetics to represent
##' phylogenetic trees, and other types of phylogenetic data. This
##' function uses NCL (the NEXUS Class Library) to parse NEXUS, Newick
##' or other common phylogenetic file formats, and returns the
##' relevant elements as a list. \code{phylo} (from the ape package)
##' or \code{phylo4} (from the phylobase package) can be constructed
##' from the elements contained in this list.
##'
##' @title Get the elements from a NEXUS (or Newick) file
##' @param file path to a NEXUS or Newick file
##' @param file.format a character string indicating the type of file
##' to be parsed.
##' @param spacesAsUnderscores In the NEXUS file format white spaces
##' are not allowed and are represented by underscores. Therefore, NCL
##' converts underscores found in taxon labels in the NEXUS file into
##' white spaces (e.g. \code{species_1} will become \code{"species
##' 1"}). If you want to preserve the underscores, set as \code{TRUE}
##' (default). This option affects taxon labels, character labels and
##' state labels.
##' @param char.all If \code{TRUE} (default), returns all characters,
##' even those excluded in the NEXUS file (only when NEXUS file
##' contains DATA block).
##' @param polymorphic.convert If TRUE (default), converts polymorphic
##' characters to missing data (only when NEXUS file contains DATA
##' block).
##' @param levels.uniform If TRUE (default), uses the same levels for
##' all characters (only when NEXUS file contains DATA block).
##' @param ... additional parameters (currently not in use).
##' @references Maddison DR, Swofford DL, Maddison WP (1997). "NEXUS:
##' An extensible file format for systematic information". Systematic
##' Biology 46(4) : 590-621.
##' doi:\href{http://dx.doi.org/10.1093/sysbio/46.4.590}{10.1093/sysbio/46.4.590}
##'
##' Lewis, P. O. 2003. NCL: a C++ class library for interpreting data
##' files in NEXUS format. Bioinformatics 19 (17) : 2330-2331.
##' @author Francois Michonneau
##' @seealso For examples on how to use the elements of the list
##' returned by this function to build tree objects, inspect the
##' source code of this package, in particular how
##' \code{read_newick_phylo} and \code{read_nexus_phylo} work. For a
##' more complex example that also use the data contained in NEXUS
##' files, inspect the source code of the \code{readNCL} function in
##' the phylobase package.
##' @return A list that contains the elements extracted from a NEXUS
##' or a Newick file.
##'
##' \itemize{
##'
##'   \item {\code{taxaNames}} {A vector of the taxa names listed in
##' the TAXA block of the NEXUS file or inferred from the tree strings
##' (if block missing or Newick file).}
##'
##'   \item {\code{treeNames}} {A vector listing the names of the trees}
##'
##'   \item {\code{taxonLabelVector}} {A list containing as many
##' elements as there are trees in the file. Each element is a
##' character vector that lists the taxon names encountered in the
##' tree string *in the order they appear*, and therefore may not
##' match the order they are listed in the translation table.}
##'
##'   \item {\code{parentVector}} { A list containing as many elements
##' as there are trees in the file. Each element is a numeric vector
##' listing the parent node for the node given by its position in the
##' vector. If the beginning of the vector is 5 5 6, the parent node
##' of node 1 is 5, the parent of node 2 is 5 and the parent of node 3
##' is 6. The implicit root of the tree is identified with 0 (node
##' without a parent).}
##'
##'   \item{\code{branchLengthVector}} { A list containing as many
##' elements as there are trees in the file. Each element is a numeric
##' vector listing the edge/branch lengths for the edges in the same
##' order as nodes are listed in the corresponding \code{parentVector}
##' element. Values of -999 indicate that the value is missing for this
##' particular edge. The implicit root as a length of 0.}
##'
##'   \item{\code{nodeLabelsVector}} { A list containing as many
##' elements as there are trees in the file. Each element is a
##' character vector listing the node labels in the same order as the
##' nodes are specified in the same order as nodes are listed in the
##' corresponding \code{parentVector} element.}
##'
##'   \item{\code{trees}} { A character vector listing the tree
##' strings where tip labels have been replaced by their indices in
##' the \code{taxaNames} vector. They do not correspond to the numbers
##' listed in the translation table that might be associated with the
##' tree.}
##'
##'   \item{\code{dataTypes}} { A character vector indicating the type
##' of data associated with the tree (e.g., \dQuote{standard}). }
##'
##'   \item{\code{nbCharacters}} { A numeric vector indicating how
##' many characters/traits are available. }
##'
##'   \item{\code{charLabels}} { A character vector listing the names
##' of the characters/traits that are available. }
##'
##'   \item {\code{nbStates}} { A numeric vector listing the number of
##' possible states for each character/trait.}
##'
##'   \item {\code{stateLabels}} { A character vector listing in
##' order, all possible states for each character/trait.}
##'
##'   \item {\code{dataChr}} { A character vector with as many
##' elements as there are characters/traits in the dataset. Each
##' element is string that can be parsed by R to create a factor
##' vector representing the data found in the file.}
##'
##'   \item {\code{isRooted}} { A list with as many elements as there
##' are trees in the file. Each element is a logical indicating
##' whether the tree is rooted. NCL definition of a rooted tree
##' differs from the one APE uses in some cases. }
##'
##'   \item {\code{hasPolytomies}} { A list with as many elements as
##' there are trees in the file. Each element is a logical indicating
##' whether the tree contains polytomies.}
##'
##'   \item {\code{hasSingletons}} { A list with as many elements as
##' there are trees in the file. Each element is a logical indicating
##' whether the tree contains singleton nodes, in other words nodes
##' with a single descendant (also known as knuckles).}
##'
##' }
##'
##'
##' @export
rncl <- function(file, file.format = c("nexus", "newick"),
                 spacesAsUnderscores = TRUE, char.all=TRUE,
                 polymorphic.convert=TRUE, levels.uniform=TRUE, ...) {

    file <- path.expand(file)
    if (!file.exists(file)) {
        stop(file, " doesn't exist.")
    }

    file.format <- match.arg(file.format)
    if (file.format == "newick") file.format <- "relaxedphyliptree"

    fileName <- list(fileName=file, fileFormat=file.format)

    ## Order of the logical parameters for GetNCL R (and C++) arguments
    ## - char.all (charall)
    ## - polymorphic.convert (polyconvert)
    ## - levels.uniform (levelsUnif)
    ## - (returnTrees)
    ## - (returnData)
    parameters <- c(char.all, polymorphic.convert, levels.uniform, TRUE, TRUE)

    ncl <- RNCL(fileName, parameters)

    ## Return Error message
    if (exists("ErrorMsg", where=ncl)) {
        stop(ncl$ErrorMsg)
    }

    if (spacesAsUnderscores) {
        ncl$taxonLabelVector <- lapply(ncl$taxonLabelVector, function(x) {
                                           gsub("\\s", "_", x)
                                       })
        ncl$taxaNames <- gsub("\\s", "_", ncl$taxaNames)

        ncl$charLabels <- gsub("\\s", "_", ncl$charLabels)

        ncl$stateLabels <- gsub("\\s", "_", ncl$stateLabels)

    }

    ncl
}

## Returns the edge matrix from the parentVector (the i^th element is
## the descendant element of node i)
get_edge_matrix <- function(parentVector) {
    edgeMat <- cbind(parentVector, 1:length(parentVector))
    rootNd <- edgeMat[which(edgeMat[, 1] == 0), 2]
    edgeMat <- edgeMat[-which(edgeMat[, 1] == 0), ]
    attr(edgeMat, "root") <- rootNd
    edgeMat
}

## Returns the edge lengths (missing are represented by -999)
get_edge_length <- function(branchLengthVector, parentVector) {
    edgeLgth <- branchLengthVector[which(parentVector != 0)]
    edgeLgth[edgeLgth == -999] <- NA
    edgeLgth
}

## Tests whether there are node labels
has_node_labels <- function(nodeLabelsVector) {
    any(nzchar(nodeLabelsVector))
}


## Pieces together the elements needed to build a phylo object, but
## they are not converted as such to allow for singletons (and
## possibly other kinds of trees that phylo doesn't support)
build_raw_phylo <- function(ncl, missing_edge_length) {
    if (length(ncl$trees) > 0) {
        listTrees <- vector("list", length(ncl$trees))
        names(listTrees) <- ncl$treeNames

        for (i in 1:length(ncl$trees)) {
            edgeMat <- get_edge_matrix(ncl$parentVector[[i]])
            rootNd <- attr(edgeMat, "root")
            attr(edgeMat, "root") <- NULL
            attr(edgeMat, "dimnames") <- NULL

            edgeLgth <- get_edge_length(ncl$branchLength[[i]], ncl$parentVector[[i]])

            tipLbl <- ncl$taxonLabelVector[[i]]

            nNodes <- length(ncl$parentVector[[i]]) - length(tipLbl)

            tr <- list(edge=edgeMat, tip.label=tipLbl, Nnode=nNodes)

            if (!all(is.na(edgeLgth))) {
                if (any(is.na(edgeLgth))) {
                    if(!(identical(length(missing_edge_length), 1L))) {
                        stop("A single numerical value should be provided for the missing edge length.")
                    }
                    if (is.na(missing_edge_length)) {
                        warning("missing edge lengths are not allowed in phylo class. All removed.")
                    } else {
                        if(!identical(mode(missing_edge_length), "numeric")) {
                            stop("A single numerical value should be provided for the missing edge lengths.")
                        }
                        warning("missing edge lengths replaced by ", sQuote(missing_edge_length), ".")
                        edgeLgth[is.na(edgeLgth)] <- missing_edge_length
                        tr <- c(tr,  list(edge.length = edgeLgth))
                    }
                } else {
                    tr <- c(tr, list(edge.length=edgeLgth))
                }
            }

            if (has_node_labels(ncl$nodeLabelsVector[[i]])) {
                ndLbl <- ncl$nodeLabelsVector[[i]]
                ndLbl[rootNd] <- ndLbl[1]
                ndLbl <- ndLbl[min(tr$edge[, 1]):length(ndLbl)]
                tr <- c(tr, list(node.label=ndLbl))
            }

            listTrees[[i]] <- tr
        }

    } else {
        return(NULL)
    }
    listTrees
}

## polishes things up
build_phylo <- function(ncl, simplify=FALSE, missing_edge_length) {
    trees <- build_raw_phylo(ncl, missing_edge_length)
    if (!is.null(trees)) {
        trees <- lapply(trees, function(tr) {
                            if (any(tabulate(tr$edge[, 1]) == 1)) {
                                tr <- collapse_singles(tr)
                            }
                            class(tr) <- "phylo"
                            tr
                        })
        if (length(trees) == 1 || simplify) {
            trees <- trees[[1]]
        } else {
            class(trees) <- "multiPhylo"
        }
    }
    trees
}

##' Create phylo objects from NEXUS or Newick files
##'
##' These functions read NEXUS or Newick files and return an object of
##' class phylo/multiPhylo.
##' @title Read phylogenetic trees from files
##' @param file Path of NEXUS or Newick file
##' @param simplify If the file includes more than one tree, returns
##'     only the first tree; otherwise, returns a multiPhylo object
##' @param missing_edge_length If the tree contains missing edge
##'     lengths, the value to be attributed to these edge lengths. By
##'     default, (\code{missing_edge_length = NA}) if at least edge
##'     length is missing, they are all removed. Otherwise, the value
##'     must be a single numeric value. In any case, a warning will
##'     be generated if the tree contains missing edge lengths.
##' @param ... additional parameters to be passed to the rncl function
##' @return A phylo or a multiPhylo object
##' @author Francois Michonneau
##' @seealso rncl-package
##' @rdname read_nexus_phylo
##' @note \code{make_phylo} will soon be deprecated, use
##' \code{read_nexus_phylo} or \code{read_newick_phylo} instead.
##' @export

make_phylo <- function(file, simplify=FALSE, missing_edge_length = NA, ...) {
    .Deprecated(msg = paste0("Use ", sQuote("read_nexus_phylo"),
                " or ", sQuote("read_newick_phylo"), " instead"))
    internal_make_phylo(file = file, simplify=simplify, missing_edge_length = missing_edge_length, ...)
}

internal_make_phylo <- function(file, simplify=FALSE, missing_edge_length = NA, ...) {
    ncl <- rncl(file=file, ...)
    build_phylo(ncl, simplify=simplify, missing_edge_length = missing_edge_length)
}

##' @rdname read_nexus_phylo
##' @export
read_nexus_phylo <- function(file, simplify=FALSE, missing_edge_length = NA, ...) {
    internal_make_phylo(file=file, simplify=simplify, file.format="nexus",
               missing_edge_length = missing_edge_length, ...)
}

##' @rdname read_nexus_phylo
##' @export
read_newick_phylo <- function(file, simplify=FALSE, missing_edge_length = NA, ...) {
    internal_make_phylo(file=file, simplify=simplify, file.format="newick",
               missing_edge_length = missing_edge_length, ...)
}
