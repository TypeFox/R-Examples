#' Function to setup a pipeine to estimate RWR-based contact strength between samples from an input gene-sample data matrix and an input graph
#'
#' \code{dRWRpipeline} is supposed to estimate sample relationships (ie. contact strength between samples) from an input gene-sample matrix and an input graph. The pipeline includes: 1) random walk restart (RWR) of the input graph using the input matrix as seeds; 2) calculation of contact strength (inner products of RWR-smoothed columns of input matrix); 3) estimation of the contact signficance by a randomalisation procedure. It supports two methods how to use RWR: 'direct' for directly applying RWR in the given seeds; 'indirectly' for first pre-computing affinity matrix of the input graph, and then deriving the affinity score. Parallel computing is also supported for Linux or Mac operating systems.
#'
#' @param data an input gene-sample data matrix used for seeds. Each value in input gene-sample matrix does not necessarily have to be binary (non-zeros will be used as a weight, but should be non-negative for easy interpretation). 
#' @param g an object of class "igraph" or "graphNEL"
#' @param method the method used to calculate RWR. It can be 'direct' for directly applying RWR, 'indirect' for indirectly applying RWR (first pre-compute affinity matrix and then derive the affinity score)
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for RWR. The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param permutation how to do permutation. It can be 'degree' for degree-preserving permutation, 'random' for permutation in random
#' @param num.permutation the number of permutations used to for generating the distribution of contact strength under randomalisation
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param adjp.cutoff the cutoff of adjusted pvalue to construct the contact graph
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' an object of class "dContact", a list with following components:
#' \itemize{
#'  \item{\code{ratio}: a symmetric matrix storing ratio (the observed against the expected) between pairwise samples}
#'  \item{\code{zscore}: a symmetric matrix storing zscore between pairwise samples}
#'  \item{\code{pval}: a symmetric matrix storing pvalue between pairwise samples}
#'  \item{\code{adjpval}: a symmetric matrix storing adjusted pvalue between pairwise samples}
#'  \item{\code{cgraph}: the constructed contact graph (as a 'igraph' object) under the cutoff of adjusted value}
#'  \item{\code{Amatrix}: a pre-computated affinity matrix when using 'inderect' method; NULL otherwise}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The choice of which method to use RWR depends on the number of seed sets and the number of permutations for statistical test. If the total product of both numbers are huge, it is better to use 'indrect' method (for a single run). However, if the user wants to re-use pre-computed affinity matrix (ie. re-use the input graph a lot), then it is highly recommended to sequentially use \code{\link{dRWR}} and \code{\link{dRWRcontact}} instead. 
#' @export
#' @import Matrix
#' @seealso \code{\link{dRWR}}, \code{\link{dRWRcontact}}, \code{\link{dCheckParallel}}
#' @include dRWRpipeline.r
#' @examples
#' \dontrun{
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#'
#' # 2) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#' V(subg)$name <- 1:vcount(subg)
#' 
#' # 3) estimate RWR dating based sample relationships
#' # define sets of seeds as data
#' # each seed with equal weight (i.e. all non-zero entries are '1')
#' aSeeds <- c(1,0,1,0,1)
#' bSeeds <- c(0,0,1,0,1)
#' data <- data.frame(aSeeds,bSeeds)
#' rownames(data) <- 1:5
#' # calcualte their two contact graph
#' dContact <- dRWRpipeline(data=data, g=subg, parallel=FALSE)
#' dContact
#' }

dRWRpipeline <- function(data, g, method=c("direct","indirect"), normalise=c("laplacian","row","column","none"), restart=0.75, normalise.affinity.matrix=c("none","quantile"), permutation=c("random","degree"), num.permutation=10, p.adjust.method=c("BH","BY","bonferroni","holm","hochberg","hommel"), adjp.cutoff=0.05, parallel=TRUE, multicores=NULL, verbose=T)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    method <- match.arg(method)    
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    permutation <- match.arg(permutation)
    p.adjust.method <- match.arg(p.adjust.method)
    
    ## check input data
    if(is.matrix(data) | is.data.frame(data)){
        data <- as.matrix(data)
        if(ncol(data)<2){
            stop("The input data must be matrix with at least two columns.\n")
        }
    }else if(is.null(data)){
        stop("The input data must be matrix.\n")
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
    ## check input graph
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    
    ####################################################
    # A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d (%s)", i, B, as.character(Sys.time())), appendLF=T)
            }
        }
    }
    
    ## A function to degree-preserving randomisation
    dp_randomisation <- function(ig, data){
        dg <- igraph::degree(ig)
        at <- unique(stats::quantile(dg, seq(from=0,to=1,by=0.1)))
        groups <- sapply(dg, function(x){
            if(length(which(x>at))==0){
                at[1]
            }else{
                at[max(which(x>at))+1]
            }
        })
        dg_random <- 1:length(dg)
        for(k in 1:length(at)){
            ind <- which(groups==at[k])
            dg_random[ind] <- ind[sample(1:length(ind))]
        }
        res <- data[dg_random,]
        rownames(res) <- rownames(data)
        return(res)
    }
    
    ## A function to make sure the sum of elements in each steady probability vector is one
    sum2one <- function(PTmatrix){
        col_sum <- apply(PTmatrix, 2, sum)
        col_sum_matrix <- matrix(rep(col_sum, nrow(PTmatrix)), ncol=ncol(PTmatrix), nrow=nrow(PTmatrix), byrow =T)
        res <- PTmatrix/col_sum_matrix
        res[is.na(res)] <- 0
        return(res)
    }
    
    ####################################################
    if(verbose){
        message(sprintf("First, RWR on input graph (%d nodes and %d edges) using input matrix (%d rows and %d columns) as seeds (%s)...", vcount(ig), ecount(ig), nrow(data), ncol(data), as.character(Sys.time())), appendLF=T)
    }
    
    ## check mapping between input data and graph
    ind <- match(rownames(data), V(ig)$name)
    nodes_mapped <- V(ig)$name[ind[!is.na(ind)]]
    P0matrix <- matrix(0, nrow=vcount(ig),ncol=ncol(data))
    P0matrix[ind[!is.na(ind)],] <- as.matrix(data[!is.na(ind),])    
    rownames(P0matrix) <- V(ig)$name
    colnames(P0matrix) <- cnames

    if(method=='indirect'){
    
        if(verbose){
            message(sprintf("\tusing %s method to do RWR (%s)...", method, as.character(Sys.time())), appendLF=T)
        }
        sAmatrix <- suppressWarnings(suppressMessages(dRWR(g=ig, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores)))
        
        PTmatrix <- sAmatrix %*% Matrix::Matrix(P0matrix, sparse=T)
        
        ## make sure the sum of elements in each steady probability vector is one
        PTmatrix <- sum2one(as.matrix(PTmatrix))
        
    }else if(method=='direct'){
    
        sAmatrix <- NULL
        
        ## RWR of a Matrix
        PTmatrix <- suppressWarnings(suppressMessages(dRWR(g=ig, normalise=normalise, setSeeds=P0matrix, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores)))
        
        PTmatrix <- Matrix::Matrix(PTmatrix, sparse=T)
    }
    
    ####################################################
    if(verbose){
        message(sprintf("Second, calculate contact strength (%s)...", as.character(Sys.time())), appendLF=T)
    }
    
    obs <- as.matrix(t(as.matrix(PTmatrix)) %*% PTmatrix)
    #diag(obs)<-0
    #obs <- 2*obs/sum(obs)
    
    B <- num.permutation
    if(verbose){
        message(sprintf("Third, generate the distribution of contact strength based on %d permutations on nodes respecting %s (%s)...", B, permutation, as.character(Sys.time())), appendLF=T)
    }
    
    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
    
        flag_parallel <- dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            b <- 1
            exp_b <- foreach::`%dopar%` (foreach::foreach(b=1:B, .inorder=T), {
                progress_indicate(b, B, 10, flag=T)
                if(permutation=="degree"){
                    seeds_random <- dp_randomisation(ig, P0matrix)
                }else if(permutation=="random"){
                    seeds_random <- P0matrix[sample(1:nrow(P0matrix)),]
                    rownames(seeds_random) <- rownames(P0matrix)
                }
                if(method=='indirect'){
                    PT_random <- sAmatrix %*% Matrix::Matrix(seeds_random, sparse=T)
                    ## make sure the sum of elements in each steady probability vector is one
                    PT_random <- sum2one(as.matrix(PT_random))
                }else if(method=='direct'){
                    PT_random <- suppressWarnings(suppressMessages(dRWR(g=ig, normalise=normalise, setSeeds=seeds_random, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores)))
                    PTmatrix <- Matrix::Matrix(PTmatrix, sparse=T)
                }
                as.matrix(t(as.matrix(PT_random)) %*% PT_random)
            })
        }
    }
    
    ###### non-parallel computing
    if(flag_parallel==F){        
        #exp_b <- list()
        exp_b <- lapply(1:B, function(b){
        
            progress_indicate(b, B, 10, flag=T)
            if(permutation=="degree"){
                seeds_random <- dp_randomisation(ig, P0matrix)
            }else if(permutation=="random"){
                seeds_random <- P0matrix[sample(1:nrow(P0matrix)),]
                rownames(seeds_random) <- rownames(P0matrix)
            }
        
            if(method=='indirect'){
                PT_random <- sAmatrix %*% Matrix::Matrix(seeds_random, sparse=T)
                ## make sure the sum of elements in each steady probability vector is one
                PT_random <- sum2one(as.matrix(PT_random))
            }else if(method=='direct'){
                PT_random <- suppressWarnings(suppressMessages(dRWR(g=ig, normalise=normalise, setSeeds=seeds_random, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix)))
                PTmatrix <- Matrix::Matrix(PTmatrix, sparse=T)
            }
            #exp_b[[b]] <- as.matrix(t(as.matrix(PT_random)) %*% PT_random)
            as.matrix(t(as.matrix(PT_random)) %*% PT_random)
        })
    }

    if(verbose){
        message(sprintf("Last, estimate the significance of contact strength: zscore, pvalue, and %s adjusted-pvalue (%s)...", p.adjust.method, as.character(Sys.time())), appendLF=T)
    }
    
    n <- ncol(obs)
    ## for zscore
    exp_mean <- matrix(0, ncol=n, nrow=n)
    exp_square <- matrix(0, ncol=n, nrow=n)
    for(b in 1:B){
        
        exp <- exp_b[[b]]
        #diag(exp)<-0
        #exp <- 2*exp/sum(exp)
    
        exp_mean <- exp_mean + exp
        exp_square <- exp_square + exp^2
    }
    exp_mean <- exp_mean/B
    exp_square <- exp_square/B
    exp_std <- sqrt(exp_square-exp_mean^2)
    zscore <- (obs-exp_mean)/exp_std
    zscore[is.na(zscore)] <- 0
    zscore[is.infinite(zscore)] <- max(zscore[!is.infinite(zscore)])
    ratio <- obs/exp_mean
    
    ## for pvalue
    num <- matrix(0, ncol=n, nrow=n)
    for(b in 1:B){
        num <- num + (obs < exp_b[[b]])
    }
    pval <- num/B
    colnames(pval) <- colnames(obs)
    rownames(pval) <- rownames(obs)
    
    ## for adjusted pvalue
    adjpval <- pval
    ## lower part
    flag_lower <- lower.tri(pval, diag=F)
    adjpval[flag_lower] <- stats::p.adjust(pval[flag_lower], method=p.adjust.method)
    ## upper part
    flag_upper <- upper.tri(pval, diag=F)
    adjpval[flag_upper] <- stats::p.adjust(pval[flag_upper], method=p.adjust.method)

    if(verbose){
        message(sprintf("Also, construct the contact graph under the cutoff %1.1e of adjusted-pvalue (%s)...", adjp.cutoff, as.character(Sys.time())), appendLF=T)
    }
    flag <- adjpval < adjp.cutoff
    adjmatrix <- flag
    adjmatrix[flag] <- zscore[flag]
    adjmatrix <- as.matrix(adjmatrix)
    cgraph <- igraph::graph.adjacency(adjmatrix, mode="undirected", weighted=T, diag=F, add.colnames=NULL, add.rownames=NA)

    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    dContact <- list(ratio      = ratio, 
                     zscore     = zscore, 
                     pval       = pval,
                     adjpval    = adjpval, 
                     cgraph     = cgraph, 
                     Amatrix    = sAmatrix,
                     call       = match.call(), 
                     method     = "dnet")
    class(dContact) <- "dContact"
    
    invisible(dContact)
}
