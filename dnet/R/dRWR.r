#' Function to implement Random Walk with Restart (RWR) on the input graph
#'
#' \code{dRWR} is supposed to implement Random Walk with Restart (RWR) on the input graph. If the seeds (i.e. a set of starting nodes) are given, it intends to calculate the affinity score of all nodes in the graph to the seeds. If the seeds are not give, it will pre-compute affinity matrix for nodes in the input graph with respect to each starting node (as a seed) by looping over every node in the graph. Parallel computing is also supported for Linux or Mac operating systems.
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param setSeeds an input matrix used to define sets of starting seeds. One column corresponds to one set of seeds that a walker starts with. The input matrix must have row names, coming from node names of input graph, i.e. V(g)$name, since there is a mapping operation. The non-zero entries mean that the corresonding rows (i.e. the gene/row names) are used as the seeds, and non-zero values can be viewed as how to weight the relative importance of seeds. By default, this option sets to "NULL", suggesting each node in the graph will be used as a set of the seed to pre-compute affinity matrix for the input graph. This default does not scale for large input graphs since it will loop over every node in the graph; however, the pre-computed affinity matrix can be extensively reused for obtaining affinity scores between any combinations of nodes/seeds, allows for some flexibility in the downstream use, in particular when sampling a large number of random node combinations for statistical testing
#' @param restart the restart probability used for RWR. The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a sparse matrix, called 'PTmatrix':
#' \itemize{
#'  \item{When the seeds are NOT given: a pre-computated affinity matrix with the dimension of n X n, where n is the number of nodes in the input graph. Columns stand for starting nodes walking from, and rows for ending nodes walking to. Therefore, a column for a starting node represents a steady-state affinity vector that the starting node will visit all the ending nodes in the graph}
#'  \item{When the seeds are given: an affinity matrix with the dimension of n X nset, where n is the number of nodes in the input graph, and nset for the number of the sets of seeds (i.e. the number of columns in setSeeds). Each column stands for the steady probability vector, storing the affinity score of all nodes in the graph to the starting nodes/seeds. This steady probability vector can be viewed as the "influential impact" over the graph imposed by the starting nodes/seeds.}
#' }
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @import Matrix
#' @seealso \code{\link{dRWRcontact}}, \code{\link{dRWRpipeline}}, \code{\link{dCheckParallel}}
#' @include dRWR.r
#' @examples
#' \dontrun{
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#'
#' # 2) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#' V(subg)$name <- 1:vcount(subg)
#'
#' # 3) obtain the pre-computated affinity matrix
#' PTmatrix <- dRWR(g=subg, normalise="laplacian", restart=0.75, parallel=FALSE)
#' # visualise affinity matrix
#' visHeatmapAdv(PTmatrix, Rowv=FALSE, Colv=FALSE, colormap="wyr", KeyValueName="Affinity")
#' 
#' # 4) obtain affinity matrix given sets of seeds
#' # define sets of seeds
#' # each seed with equal weight (i.e. all non-zero entries are '1')
#' aSeeds <- c(1,0,1,0,1)
#' bSeeds <- c(0,0,1,0,1)
#' setSeeds <- data.frame(aSeeds,bSeeds)
#' rownames(setSeeds) <- 1:5
#' # calcualte affinity matrix
#' PTmatrix <- dRWR(g=subg, normalise="laplacian", setSeeds=setSeeds, restart=0.75, parallel=FALSE)
#' PTmatrix
#' }

dRWR <- function(g, normalise=c("laplacian","row","column","none"), setSeeds=NULL, restart=0.75, normalise.affinity.matrix=c("none","quantile"), parallel=TRUE, multicores=NULL, verbose=T)
{
    
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }

    if(verbose){
        now <- Sys.time()
        message(sprintf("First, get the adjacency matrix of the input graph (%s) ...", as.character(now)), appendLF=T)
    }
    
    if ("weight" %in% list.edge.attributes(ig)){
        adjM <- get.adjacency(ig, type="both", attr="weight", edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        if(verbose){
            message(sprintf("\tNotes: using weighted graph!"), appendLF=T)
        }
    }else{
        adjM <- get.adjacency(ig, type="both", attr=NULL, edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        if(verbose){
            message(sprintf("\tNotes: using unweighted graph!"), appendLF=T)
        }
    }
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Then, normalise the adjacency matrix using %s normalisation (%s) ...", normalise, as.character(now)), appendLF=T)
    }
    
    A <- adjM!=0
    ## library(Matrix) 
    ## vignette("Intro2Matrix")
    if(normalise == "row"){
        #D <- diag(apply(A,1,sum)^(-1))
        D <- Matrix::Diagonal(x=(Matrix::rowSums(A))^(-1))
        nadjM <- adjM %*% D
    }else if(normalise == "column"){
        #D <- diag(apply(A,1,sum)^(-1))
        D <- Matrix::Diagonal(x=(Matrix::colSums(A))^(-1))
        nadjM <- D %*% adjM
    }else if(normalise == "laplacian"){
        #D <- diag(apply(A,1,sum)^(-0.5))
        D <- Matrix::Diagonal(x=(Matrix::colSums(A))^(-0.5))
        nadjM <- D %*% adjM %*% D
    }else{
        nadjM <- adjM
    }
    #nadjM <- as.matrix(nadjM)
    
    #nig <- graph.adjacency(nadjM, mode=c("undirected"), weighted=T, diag=F, add.colnames=NULL, add.rownames=NA)
    ## update the vertex attributes
    nattr <- list.vertex.attributes(ig)
    for(attr in nattr){
        #nig <- set.vertex.attribute(nig, name=attr, index=V(nig), get.vertex.attribute(ig,attr))
    }
    ## update the edge attributes    
    eattr <- list.edge.attributes(ig)
    for(attr in eattr){
        if (!("weight" %in% attr)){
            #nig <- set.edge.attribute(nig, name=attr, index=E(nig), get.edge.attribute(ig,attr))
        }
    }
    
    
    ####################################################
    # A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d seed sets (%s)", i, B, as.character(Sys.time())), appendLF=T)
            }
        }
    }
    
    ## A function to make sure the sum of elements in each steady probability vector is one
    sum2one <- function(PTmatrix){
        col_sum <- apply(PTmatrix, 2, sum)
        #col_sum <- Matrix::colSums(PTmatrix, sparseResult=F)
        col_sum_matrix <- matrix(rep(col_sum, nrow(PTmatrix)), ncol=ncol(PTmatrix), nrow=nrow(PTmatrix), byrow =T)
        res <- as.matrix(PTmatrix)/col_sum_matrix
        res[is.na(res)] <- 0
        return(res)
    }
    
    ## a function to normalize columns of a matrix to have the same Quantiles
    normalizeQuantiles <- function (A, ties=TRUE) {
        n <- dim(A)
        if(is.null(n)) return(A)
        if(n[2] == 1) return(A)
        O <- S <- array(, n)
        nobs <- rep(n[1], n[2])
        i <- (0:(n[1] - 1))/(n[1] - 1)
        for(j in 1:n[2]){
            Si <- sort(A[, j], method = "quick", index.return = TRUE)
            nobsj <- length(Si$x)
            if (nobsj < n[1]){
                nobs[j] <- nobsj
                isna <- is.na(A[, j])
                S[, j] <- stats::approx((0:(nobsj - 1))/(nobsj - 1), Si$x, 
                    i, ties = "ordered")$y
                O[!isna, j] <- ((1:n[1])[!isna])[Si$ix]
            }else{
                S[, j] <- Si$x
                O[, j] <- Si$ix
            }
        }
        m <- rowMeans(S)
        for (j in 1:n[2]){
            if(ties) r<-rank(A[, j])
            if(nobs[j] < n[1]){
                isna <- is.na(A[, j])
                if(ties){ 
                    A[!isna, j] <- stats::approx(i, m, (r[!isna] - 1)/(nobs[j] - 1), ties = "ordered")$y
                }else{ 
                    A[O[!isna, j], j] <- stats::approx(i, m, (0:(nobs[j] - 1))/(nobs[j] - 1), ties = "ordered")$y
                }
            }else{
                if(ties){
                    A[, j] <- stats::approx(i, m, (r - 1)/(n[1] - 1), ties = "ordered")$y
                }else{
                    A[O[, j], j] <- m
                }
            }
        }
    
        return(A)
    }
    
    ####################################################
    
    ################## RWR
    ## restarting prob
    if(is.null(restart) || is.na(restart) || restart<0 || restart>100){
        r <- 0.75
    }else if(restart>1 && restart<100){
        r <- restart/100
    }else{
        r <- restart
    }
    ## stopping critera
    stop_delta <- 1e-6   # L1 norm of successive estimates 'PT' below the threshold 'stop_delta'
    stop_step <- 50      # maximum steps of iterations
    
    if(is.null(setSeeds)){
        P0matrix <- Matrix::Matrix(diag(vcount(ig)), sparse=T)
        rownames(P0matrix) <- V(ig)$name
        colnames(P0matrix) <- V(ig)$name
        
    }else{
        ## check input data
        if(is.matrix(setSeeds) | is.data.frame(setSeeds)){
            data <- as.matrix(setSeeds)
        }else if(is.vector(setSeeds)){
            data <- as.matrix(setSeeds, ncol=1)
        }

        if(is.null(rownames(data))) {
            stop("The function must require the row names of the input setSeeds.\n")
        }else if(any(is.na(rownames(data)))){
            warning("setSeeds with NA as row names will be removed")
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
            warning("The row names of input setSeeds do not contain all those in the input graph.\n")
        }
        P0matrix <- matrix(0,nrow=nrow(nadjM),ncol=ncol(data))
        P0matrix[ind[!is.na(ind)],] <- as.matrix(data[!is.na(ind),])
        
        ## make sure the sum of elements in each steady probability vector is one
        P0matrix <- sum2one(P0matrix)
        
        if(0){
        ## make sure the sum of elements in each starting probability vector is one
        P0matrix <- sapply(1:ncol(P0matrix), function(i){
            if(sum(P0matrix[,i]!=0)){
                P0matrix[,i]/sum(P0matrix[,i])
            }else{
                P0matrix[,i]
            }
        })
        rownames(P0matrix) <- V(ig)$name
        colnames(P0matrix) <- cnames
        }
        
        ## convert to sparse matrix
        P0matrix <- Matrix::Matrix(P0matrix, sparse=T)
    }
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Third, RWR of %d sets of seeds using %1.1e restart probability (%s) ...", ncol(P0matrix), restart, as.character(now)), appendLF=T)
    }
    
    if(restart==1){
        ## just seeds themselves
        PTmatrix <- P0matrix
    }else{
        
        ###### parallel computing
        flag_parallel <- F
        if(parallel==TRUE){

            flag_parallel <- dCheckParallel(multicores=multicores, verbose=verbose)
            if(flag_parallel){
                j <- 1
                PTmatrix <- foreach::`%dopar%` (foreach::foreach(j=1:ncol(P0matrix), .inorder=T, .combine='cbind'), {
                    progress_indicate(j, ncol(P0matrix), 10, flag=T)
                    P0 <- P0matrix[,j]
                    ## Initializing variables
                    step <- 0
                    delta <- 1
                    PT <- P0
                    ## Iterative update till convergence (delta<=1e-10)
                    while (delta>stop_delta && step<=stop_step){
                        PX <- (1-r) * nadjM %*% PT + r * P0
                        # p-norm of v: sum((abs(v).p)^(1/p))
                        delta <- sum(abs(PX-PT))
                        PT <- PX
                        step <- step+1
                    }
                    as.matrix(PT)
                })
                    
                PTmatrix[PTmatrix<1e-6] <- 0
                #PTmatrix <- Matrix::Matrix(PTmatrix, sparse=T)
            }
        }
        
        ###### non-parallel computing
        if(flag_parallel==F){
            PTmatrix <- Matrix::Matrix(0, nrow=nrow(P0matrix), ncol=ncol(P0matrix), sparse=T)
            for(j in 1:ncol(P0matrix)){
                #P0 <- as.matrix(P0matrix[,j],ncol=1)
                P0 <- P0matrix[,j]
        
                ## Initializing variables
                step <- 0
                delta <- 1
        
                PT <- P0
                ## Iterative update till convergence (delta<=1e-10)
                while (delta>stop_delta && step<=stop_step){
                    PX <- (1-r) * nadjM %*% PT + r * P0
    
                    # p-norm of v: sum((abs(v).p)^(1/p))
                    delta <- sum(abs(PX-PT))
    
                    PT <- PX
                    step <- step+1
                }
                #PTmatrix[,j] <- as.matrix(PT, ncol=1)
                PT[PT<1e-6] <- 0
                #PTmatrix[,j] <- Matrix::Matrix(PT, sparse=T)
                PTmatrix[,j] <- PT
            
                progress_indicate(j, ncol(P0matrix), 10, flag=T)
        
            }
        }
    }
    
    ## make sure the sum of elements in each steady probability vector is one
    if(verbose){
        now <- Sys.time()
        message(sprintf("Fourth, rescale steady probability vector (%s) ...", as.character(now)), appendLF=T)
    }
    PTmatrix <- sum2one(PTmatrix) # input/output: full matrix
    
    if(0){
        ## make sure the sum of elements in each steady probability vector is one
        PTmatrix <- sapply(1:ncol(PTmatrix), function(i){
            if(sum(PTmatrix[,i])!=0){
                PTmatrix[,i]/sum(PTmatrix[,i])
            }else{
                PTmatrix[,i]
            }
        })
        rownames(PTmatrix) <- rownames(P0matrix)
        colnames(PTmatrix) <- colnames(P0matrix)
    }
    
    ####################################################################################
    if(ncol(PTmatrix) == 1){
        normalise.affinity.matrix <- "none"
    }
    if(normalise.affinity.matrix=="quantile"){
        PTmatrix <- normalizeQuantiles(PTmatrix)
    }
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Finally, output %d by %d affinity matrix normalised by %s (%s) ...", nrow(PTmatrix), ncol(PTmatrix), normalise.affinity.matrix, as.character(now)), appendLF=T)
    }
    
    PTmatrix[PTmatrix<1e-6] <- 0
    #PTmatrix <- signif(PTmatrix, digits=7)
    PTmatrix <- Matrix::Matrix(PTmatrix, sparse=T)
    rownames(PTmatrix) <- rownames(P0matrix)
    colnames(PTmatrix) <- colnames(P0matrix)
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)

    invisible(PTmatrix)
}


