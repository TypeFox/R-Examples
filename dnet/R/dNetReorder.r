#' Function to reorder the multiple graph colorings within a sheet-shape rectangle grid
#'
#' \code{dNetReorder} is reorder the multiple graph colorings within a sheet-shape rectangle grid
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param data an input data matrix used to color-code vertices/nodes. One column corresponds to one graph node coloring. The input matrix must have row names, and these names should include all node names of input graph, i.e. V(g)$name, since there is a mapping operation. After mapping, the length of the patern vector should be the same as the number of nodes of input graph. The way of how to color-code is to map values in the pattern onto the whole colormap (see the next arguments: colormap, ncolors, zlim and colorbar)
#' @param feature the type of the features used. It can be one of either 'edge' for the edge feature or 'node' for the node feature. See 'Note' for explanations. 
#' @param node.normalise the normalisation of the nodes. It can be one of either 'none' for no normalisation or 'degree' for a node being penalised by its degree.
#' @param xdim an integer specifying x-dimension of the grid
#' @param ydim an integer specifying y-dimension of the grid
#' @param amplifier an integer specifying the amplifier (3 by default) of the number of component planes. The product of the component number and the amplifier constitutes the number of rectangles in the sheet grid
#' @param metric distance metric used to define the similarity between component planes. It can be "none", which means directly using column-wise vectors of codebook/data matrix. Otherwise, first calculate the covariance matrix from the codebook/data matrix. The distance metric used for calculating the covariance matrix between component planes can be: "pearson" for pearson correlation, "spearman" for spearman rank correlation, "kendall" for kendall tau rank correlation, "euclidean" for euclidean distance, "manhattan" for cityblock distance, "cos" for cosine similarity, "mi" for mutual information.
#' @param init an initialisation method. It can be one of "uniform", "sample" and "linear" initialisation methods
#' @param algorithm the training algorithm. Currently, only "sequential" algorithm has been implemented
#' @param alphaType the alpha type. It can be one of "invert", "linear" and "power" alpha types
#' @param neighKernel the training neighbor kernel. It can be one of "gaussian", "bubble", "cutgaussian", "ep" and "gamma" kernels
#' @return 
#' an object of class "sReorder", a list with following components:
#' \itemize{
#'  \item{\code{nHex}: the total number of rectanges in the grid}
#'  \item{\code{xdim}: x-dimension of the grid}
#'  \item{\code{ydim}: y-dimension of the grid}
#'  \item{\code{uOrder}: the unique order/placement for each component plane that is reordered to the "sheet"-shape grid with rectangular lattice}
#'  \item{\code{coord}: a matrix of nHex x 2, with each row corresponding to the coordinates of each "uOrder" rectangle in the 2D map grid}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note 
#' According to which features are used and whether nodes should be penalised by degrees, the feature data are constructed differently from the input data and input graph:
#' \itemize{
#' \item{When the node features are used, the feature data is the input data (or penalised data) with the same dimension.}
#' \item{When the edge featrues are used, each entry (i.e. given an edge and a sample) in the feature data is the absolute difference between its two-end nodes (or after being penalised).}
#' \item{After that, the constructed feature are subject to sample correlation analysis by supraHex. That is, a map grid (with sheet shape consisting of a rectangular lattice) is used to train either column-wise vectors of the feature data matrix or the covariance matrix thereof.}
#' \item{As a result, similar samples are placed closer to each other within this map grid. More precisely, to ensure the unique placement, each sample mapped to the "sheet"-shape grid with rectangular lattice is determinied iteratively in an order from the best matched to the next compromised one. If multiple samples are hit in the same rectangular lattice, the worse one is always sacrificed by moving to the next best one till all samples are placed somewhere exclusively on their own.}
#' }
#' The size of "sheet"-shape rectangle grid depends on the input arguments: 
#' \itemize{
#' \item{How the input parameters are used to determine nHex is taken priority in the following order: "xdim & ydim" > "nHex" > "data".}
#' \item{If both of xdim and ydim are given, \eqn{nHex=xdim*ydim}.}
#' \item{If only data is input, \eqn{nHex=5*sqrt(dlen)}, where dlen is the number of rows of the input data.}
#' \item{After nHex is determined, xy-dimensions of rectangle grid are then determined according to the square root of the two biggest eigenvalues of the input data.}
#' }
#' @export
#' @seealso \code{\link{visNetReorder}}
#' @include dNetReorder.r
#' @examples
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#'
#' # 2) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#'
#' # 3) reorder the module with vertices being color-coded by input data
#' nnodes <- vcount(subg)
#' nsamples <- 10
#' data <- matrix(runif(nnodes*nsamples), nrow=nnodes, ncol=nsamples)
#' rownames(data) <- V(subg)$name
#' sReorder <- dNetReorder(g=subg, data, feature="node", node.normalise="none")

dNetReorder <- function (g, data, feature=c("node", "edge"), node.normalise=c("none", "degree"), xdim=NULL, ydim=NULL, amplifier=NULL, metric=c("none","pearson","spearman","kendall","euclidean","manhattan","cos","mi"), init=c("linear","uniform","sample"), algorithm=c("sequential","batch"), alphaType=c("invert","linear","power"), neighKernel=c("gaussian","bubble","cutgaussian","ep","gamma"))
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    feature <- match.arg(feature)
    node.normalise <- match.arg(node.normalise)
    
    metric <- match.arg(metric)
    init <- match.arg(init)
    algorithm <- match.arg(algorithm)
    alphaType <- match.arg(alphaType)
    neighKernel <- match.arg(neighKernel)

    ## check input graph
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    if(is.null(V(ig)$name)){
        V(ig)$name <- as.character(V(ig))
    }

    ## check input data
    if(is.matrix(data) | is.data.frame(data) | is.vector(data)){
        data <- as.matrix(data)
    }else if(is.null(data)){
        stop("The input data must be not NULL.\n")
    }

    if(is.null(rownames(data))) {
        stop("The function must require the row names of the input data.\n")
    }else if(any(is.na(rownames(data)))){
        warning("Data with NA as row names will be removed")
        data <- data[!is.na(rownames(data)),]
    }
    cnames <- colnames(data)
    if(is.null(cnames)){
        cnames <- seq(1,ncol(data))
    }
    
    ## check mapping between input data and graph
    ind <- match(rownames(data), V(ig)$name)
    nodes_mapped <- V(ig)$name[ind[!is.na(ind)]]
    if(length(nodes_mapped)!=vcount(ig)){
        stop("The function must require that the row names of input data could all be mapped onto the input graph.\n")
    }
    data <- as.matrix(data[nodes_mapped,])
    
    ######################################################################################
    
    if(node.normalise=="degree"){
        node.degree <- igraph::degree(ig)
    }else if(node.normalise=="none"){
        node.degree <- sapply(igraph::degree(ig), function(x) ifelse(x>0,1,0))
    }
    
    if(feature=="node"){
    	# force those isolated genes to have degree=1
    	node.degree[node.degree==0] <- 1
        fdata <- 1/node.degree * data
    }else if(feature=="edge"){
        tmp_n1 <- get.edgelist(ig,names=T)[,1]
        tmp_n2 <- get.edgelist(ig,names=T)[,2]
        fdata <- abs(1/node.degree[tmp_n1] * data[tmp_n1,] - 1/node.degree[tmp_n2] * data[tmp_n2,])
        rownames(fdata) <- E(ig)
    }

    ## define the topology of a map grid (with "sheet" shape consisting of "rect" lattice)
    if(is.null(amplifier)){
        amplifier <- 3
    }else if (amplifier <= 2){
        amplifier <- 2
    }
    nHex <- ceiling(amplifier*ncol(data))
    if(is.null(xdim) | is.null(ydim)){
        xdim <- ceiling(sqrt(nHex))
        ydim <- ceiling(sqrt(nHex))
    }else if(xdim*ydim < ncol(data)){
        xdim <- ceiling(sqrt(nHex))
        ydim <- ceiling(sqrt(nHex))
    }

    sReorder <- sCompReorder(fdata, xdim=xdim, ydim=ydim, amplifier=amplifier, metric=metric, init=init, algorithm=algorithm, alphaType=alphaType, neighKernel=neighKernel)
    
    invisible(sReorder)
}