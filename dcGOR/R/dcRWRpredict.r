#' Function to perform RWR-based ontology term predictions from input known annotations and an input graph
#'
#' \code{dcRWRpredict} is supposed to perform ontology term predictions based on Random Walk with Restart (RWR) from input known annotations and an input graph.
#'
#' @param data an input gene-term data matrix containing known annotations used for seeds. Each value in input matrix does not necessarily have to be binary (non-zeros will be used as a weight, but should be non-negative for easy interpretation). Also, data can be a list, each containing the known annotated genes
#' @param g an object of class "igraph" or \code{\link{Dnetwork}}
#' @param output.file an output file containing predicted results. If not NULL, a tab-delimited text file will be also written out; otherwise, there is no output file (by default)
#' @param ontology the ontology identity. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPON" for Human Phenotype ONset and clinical course, "MP" for Mammalian Phenotype, "EC" for Enzyme Commission, "KW" for UniProtKB KeyWords, "UP" for UniProtKB UniPathway. For details on the eligibility for pairs of input domain and ontology, please refer to the online Documentations at \url{http://supfam.org/dcGOR/docs.html}. If NA, then the user has to input a customised RData-formatted file (see \code{RData.ontology.customised} below)
#' @param method the method used to calculate RWR. It can be 'direct' for directly applying RWR, 'indirect' for indirectly applying RWR (first pre-compute affinity matrix and then derive the affinity score)
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for RWR. The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param leave.one.out logical to indicate whether the leave-one-out test is used for predictions. By default, it sets to true for doing leave-one-out test (that is, known seeds are removed)
#' @param propagation how to propagate the score. It can be "max" for retaining the maximum score (by default), "sum" for additively accumulating the score
#' @param scale.method the method used to scale the predictive scores. It can be: "none" for no scaling, "linear" for being linearily scaled into the range between 0 and 1, "log" for the same as "linear" but being first log-transformed before being scaled. The scaling between 0 and 1 is done via: \eqn{\frac{S - S_{min}}{S_{max} - S_{min}}}, where \eqn{S_{min}} and \eqn{S_{max}} are the minimum and maximum values for \eqn{S}
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.ontology.customised a file name for RData-formatted file containing an object of S4 class 'Onto' (i.g. ontology). By default, it is NULL. It is only needed when the user wants to perform customised analysis using their own ontology. See \code{\link{dcBuildOnto}} for how to creat this object
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{dcRDataLoader}} for details
#' @return 
#' a data frame containing three columns: 1st column the same as the input file (e.g. 'SeqID'), 2nd for 'Term' (predicted ontology terms), 3rd for 'Score' (along with predicted scores)
#' @note
#' When 'output.file' is specified, a tab-delimited text file is written out, with the column names: 1st column the same as the input file (e.g. 'SeqID'), 2nd for 'Term' (predicted ontology terms), 3rd for 'Score' (along with predicted scores).
#' The choice of which method to use RWR depends on the number of seed sets and whether using leave-one-out test. If the total product of both numbers are huge, it is better to use 'indrect' method (for a single run). Also, when using leave-one-out test, it has to be use 'indrect' method.
#' @export
#' @importFrom dnet dRWR dCheckParallel
#' @import Matrix
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcAlgoPropagate}}, \code{\link{dcList2Matrix}}
#' @include dcRWRpredict.r
#' @examples
#' \dontrun{
#' # 1) define an input network
#' ## 1a) an igraph object that contains a functional protein association network in human.
#' ### The network is extracted from the STRING database (version 9.1). 
#' ### Only those associations with medium confidence (score>=400) are retained
#' org.Hs.string <- dnet::dRDataLoader(RData='org.Hs.string')
#' ## 1b) restrict to those edges with confidence score>=999
#' ### keep the largest connected component
#' network <- igraph::subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=999])
#' g <- dnet::dNetInduce(g=network, nodes_query=V(network)$name, largest.comp=TRUE)
#' ## Notably, in reality, 1b) can be replaced by: 
#' #g <- igraph::subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=400])
#' ## 1c) make sure there is a 'weight' edge attribute
#' E(g)$weight <- E(g)$combined_score
#' ### use EntrezGene ID as default 'name' node attribute
#' V(g)$name <- V(g)$geneid
#' g
#'
#' # 2) define the known annotations as seeds
#' anno.file <- "http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt"
#' data <- dcSparseMatrix(anno.file)
#'
#' # 3) perform RWR-based ontology term predictions
#' res <- dcRWRpredict(data=data, g=g, ontology="HPPA", parallel=FALSE)
#' res[1:10,]

#' # 4) calculate Precision and Recall
#' GSP.file <- anno.file
#' prediction.file <- res
#' res_PR <- dcAlgoPredictPR(GSP.file=GSP.file, prediction.file=prediction.file, ontology="HPPA")
#' res_PR
#' 
#' # 5) Plot PR-curve
#' plot(res_PR[,2], res_PR[,1], xlim=c(0,1), ylim=c(0,1), type="b", xlab="Recall", ylab="Precision")
#' }

dcRWRpredict <- function(data, g, output.file=NULL, ontology=c(NA,"GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPON","MP","EC","KW","UP"), method=c("indirect","direct"), normalise=c("laplacian","row","column","none"), restart=0.75, normalise.affinity.matrix=c("none","quantile"), leave.one.out=T, propagation=c("max","sum"), scale.method=c("log","linear","none"), parallel=TRUE, multicores=NULL, verbose=T, RData.ontology.customised=NULL, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/dcGOR")
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    ontology <- match.arg(ontology)
    method <- match.arg(method)    
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    scale.method <- match.arg(scale.method)
    
    ## check input data
    if(is.matrix(data) | is.data.frame(data)){
        data <- as.matrix(data)
    }else if(is.list(data)){
        if(is.null(names(data))){
            names(data) <- seq(1, length(data))
        }
        data_names <- names(data)
        res_list <- lapply(1:length(data), function(i){
            x <- data[[i]]
            if(length(x) > 0){
                return( cbind(x, rep(data_names[i],length(x))) )
            }else{
                return(NULL)
            }
        })
        res_mat <- base::do.call(base::rbind, res_list)
        data <- suppressWarnings(dcSparseMatrix(res_mat))
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
    ig <- g
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' object.\n")
    }
    
    ####################################################
    if(verbose){
        message(sprintf("First, RWR on input graph (%d nodes and %d edges) using input matrix (%d rows and %d columns) as seeds (%s)...", vcount(ig), ecount(ig), nrow(data), ncol(data), as.character(Sys.time())), appendLF=T)
    }
    
    ## check mapping between input data and graph
    P0matrix <- matrix(0, nrow=vcount(ig),ncol=ncol(data))
    ind <- match(rownames(data), V(ig)$name)
    P0matrix[ind[!is.na(ind)],] <- as.matrix(data[!is.na(ind),])    
    rownames(P0matrix) <- V(ig)$name
    colnames(P0matrix) <- colnames(data)
    P0matrix <- Matrix::Matrix(P0matrix, sparse=T)
    
    if(leave.one.out){
        method <- 'indirect'
    }
    
    if(verbose){
        if(leave.one.out){
            message(sprintf("\tusing '%s' method to do RWR with leave-one-out test (%s)...", method, as.character(Sys.time())), appendLF=T)
        }else{
            message(sprintf("\tusing '%s' method to do RWR without leave-one-out test (%s)...", method, as.character(Sys.time())), appendLF=T)
        }
        
        message(paste(c("\n##############################"), collapse=""), appendLF=T)
        message(paste(c("'dnet::dRWR' is being called..."), collapse=""), appendLF=T)
        message(paste(c("##############################\n"), collapse=""), appendLF=T)
        
    }
    
    if(method=='indirect'){
        sAmatrix <- suppressWarnings(dnet::dRWR(g=ig, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores))
        
        sAmatrix <- as.matrix(sAmatrix)
        if(leave.one.out){
            # remove self-self affinity
            base::diag(sAmatrix) <- 0
            
            ## calculate summed affinity matrix (ASmatrix: seqID X termID)
            ASmatrix <- sAmatrix %*% P0matrix
        
            ## obtain averaged affinity matrix (AAmatrix: seqID X termID)
            col_sum <- Matrix::colSums(P0matrix)
            col_sum_matrix <- matrix(rep(col_sum,nrow(ASmatrix)), ncol=ncol(ASmatrix), nrow=nrow(ASmatrix), byrow=T)
            ##################
            ### for seeds, the average should be 1-less
            col_sum_adj_matrix <- col_sum_matrix - as.matrix(P0matrix)
            AAmatrix <- ASmatrix / col_sum_adj_matrix
            ##################
            
            AAmatrix[is.na(AAmatrix)] <- 0
            AAmatrix <- Matrix::Matrix(AAmatrix, sparse=T)
            
        }else{
            ## calculate summed affinity matrix (ASmatrix: seqID X termID)
            ASmatrix <- sAmatrix %*% P0matrix
        
            ## obtain averaged affinity matrix (AAmatrix: seqID X termID)
            col_sum <- Matrix::colSums(P0matrix)
            col_sum_matrix <- matrix(rep(col_sum,nrow(ASmatrix)), ncol=ncol(ASmatrix), nrow=nrow(ASmatrix), byrow=T)
            AAmatrix <- ASmatrix / col_sum_matrix
            
            AAmatrix[is.na(AAmatrix)] <- 0
            AAmatrix <- Matrix::Matrix(AAmatrix, sparse=T)
        }
        
    }else if(method=='direct'){
        ## RWR of a Matrix
        AAmatrix <- suppressWarnings(dnet::dRWR(g=ig, normalise=normalise, setSeeds=P0matrix, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores))
        
        colnames(AAmatrix) <- colnames(P0matrix)
        rownames(AAmatrix) <- rownames(P0matrix)
    }
    
    ####################################################
    if(verbose){
        message(paste(c("##############################"), collapse=""), appendLF=T)
        message(paste(c("'dnet::dRWR' has been completed!"), collapse=""), appendLF=T)
        message(paste(c("##############################\n"), collapse=""), appendLF=T)
    
        message(sprintf("Second, propagate '%s' ontology annotations via '%s' operation (%s)...", ontology, propagation, as.character(Sys.time())), appendLF=T)
    }
    
    seq_names <- rownames(AAmatrix)
    term_names <- colnames(AAmatrix)
    idx <- which(AAmatrix>1e-6, arr.ind=TRUE, useNames=FALSE)
    output <- cbind(SeqID=seq_names[idx[,1]], Term=term_names[idx[,2]], Score=signif(AAmatrix[idx], digits=4))
    
    if(verbose){
        message(paste(c("\n##############################"), collapse=""), appendLF=T)
        message(paste(c("'dcAlgoPropagate' is being called..."), collapse=""), appendLF=T)
        message(paste(c("##############################\n"), collapse=""), appendLF=T)
    }
    
    res_RData <- dcAlgoPropagate(input.file=output, ontology=ontology, propagation=propagation, output.file=NA, verbose=verbose, RData.ontology.customised=RData.ontology.customised, RData.location=RData.location)
    score <- res_RData$hscore
    
    if(verbose){
        message(paste(c("##############################"), collapse=""), appendLF=T)
        message(paste(c("'dcAlgoPropagate' has been completed!"), collapse=""), appendLF=T)
        message(paste(c("##############################\n"), collapse=""), appendLF=T)
    
        message(sprintf("Third, rescale predictive score using '%s' method (%s)...", scale.method, as.character(Sys.time())), appendLF=T)
    }
    
    if(scale.method=='none'){
        score_scaled <- score
    }else if(scale.method=='linear'){
        score_scaled <- lapply(score, function(x){
            y <- (x-min(x))/diff(range(x))
            y <- y[!is.na(y)]
            signif(y, digits=4)
        })
    }else if(scale.method=='log'){
        score_scaled <- lapply(score, function(x){
            x <- x[x>0] ## make sure that all scores are positive
            x <- log(x)
            y <- (x-min(x))/diff(range(x))
            y <- y[!is.na(y)]
            signif(y, digits=4)
        })
    }
    
    ## output the prediction
    output <- dcList2Matrix(score_scaled)
    colnames(output) <- c("SeqID", "Term", "Score")
    if(!is.null(output)){
        if(!is.null(output.file)){
            utils::write.table(output, file=output.file, quote=F, row.names=F, sep="\t")
            if(file.exists(output.file)){
                message(sprintf("The predictions have been saved into '%s'.", file.path(getwd(),output.file)), appendLF=T)
            }
        }
    }else{
        return(NULL)
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)

    invisible(output)
}
