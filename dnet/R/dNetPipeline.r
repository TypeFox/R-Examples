#' Function to setup the pipeline for finding maximum-scoring subgraph from an input graph and the signficance imposed on its nodes
#'
#' \code{dNetPipeline} is supposed to finish ab inito maximum-scoring subgraph identification for the input graph with the node information on the significance (p-value or fdr). It returns an object of class "igraph" or "graphNEL". 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param pval a vector containing input p-values (or fdr). For each element, it must have the name that could be mapped onto the input graph. Also, the names in input "pval" should contain all those in the input graph "g", but the reverse is not necessary
#' @param method the method used for the transformation. It can be either "pdf" for the method based on the probability density function of the fitted model, or "cdf" for the method based on the cumulative distribution function of the fitted model
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint. If given, those p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively
#' @param nsize the desired number of nodes constrained to the resulting subgraph. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subgraph. Notably, the given significance threshold will be overwritten by this option. 
#' @param plot logical to indicate whether the histogram plot, contour plot and scatter plot should be drawn. By default, it sets to false for no plotting
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a subgraph with a maximum score, an object of class "igraph" or "graphNEL"
#' @note The pipeline sequentially consists of: 
#' \itemize{
#' \item{ia) if the method is either "pdf" or "cdf", \code{\link{dBUMfit}} used to fit the p-value distribution under beta-uniform mixture model, and \code{\link{dBUMscore}} used to calculate the scores according to the fitted BUM and the significance threshold.}
#' \item{ib) if the method is either "customised", then the user input list of fdr (or p-values) and the significance threshold will be directly used for score transformation by \code{\link{dFDRscore}}.}
#' \item{ii) if there is the desired number of nodes constrained to the resulting subgraph, a wide range of significance thresholds (including rough stage with large intervals, and finetune stage with smaller intervals) will be scanned to find the significance threshold to meet the desired number of nodes.}
#' \item{iii) \code{\link{dNetFind}} used to find maximum-scoring subgraph from the input graph and scores imposed on its nodes.}
#' }
#' @export
#' @seealso \code{\link{dBUMfit}}, \code{\link{dBUMscore}}, \code{\link{dFDRscore}}, \code{\link{dNetFind}}
#' @include dBUMfit.r
#' @examples
#' \dontrun{
#' # 1) generate an vector consisting of random values from beta distribution
#' x <- rbeta(1000, shape1=0.5, shape2=1)
#' names(x) <- as.character(1:length(x))
#'
#' # 2) generate a random graph according to the ER model
#' g <- erdos.renyi.game(1000, 1/100)
#'
#' # 3) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#'
#' # 4) find maximum-scoring subgraph based on the given significance threshold
#' # 4a) assume the input is a list of p-values (controlling fdr=0.1)
#' subgraph <- dNetPipeline(g=subg, pval=x, significance.threshold=0.1)
#' # 4b) assume the input is a list of customised significance (eg FDR directly)
#' subgraph <- dNetPipeline(g=subg, pval=x, method="customised", significance.threshold=0.1)
#' 
#' # 5) find maximum-scoring subgraph with the desired node number nsize=20
#' subgraph <- dNetPipeline(g=subg, pval=x, nsize=20)
#' }

dNetPipeline <- function(g, pval, method=c("pdf","cdf","customised"), significance.threshold=NULL, nsize=NULL, plot=F, verbose=T)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    method <- match.arg(method)
    
    # force those zeros to be miminum of non-zeros
    tmp <- as.numeric(format(.Machine)['double.xmin'])
    pval[pval < tmp] <- tmp
    
    ####################
    
    if(method!="customised"){
        if(verbose){
            message(sprintf("First, fit the input p-value distribution under beta-uniform mixture model..."), appendLF=T)
        }
        ## fit a p-value distribution under beta-uniform mixture model
        if(plot){
            fit <- dBUMfit(pval, ntry=1, hist.bum=T, contour.bum=T, verbose=verbose)
        }else{
            fit <- dBUMfit(pval, ntry=1, hist.bum=F, contour.bum=F, verbose=verbose)
        }
    }else if(method=="customised"){
        if(verbose){
            message(sprintf("First, consider the input fdr (or p-value) distribution"), appendLF=T)
        }
    }
    
    if(verbose){
        message(sprintf("Second, determine the significance threshold..."), appendLF=T)
    }
    ## Determine the final fdr threshold
    if(is.null(nsize)){
        fdr_final <- significance.threshold
    }else{
        fdr_final <- NULL
        ####################################
        if(verbose){
            message(sprintf("\tVia constraint on the size of subnetwork to be identified (%d nodes)", nsize), appendLF=T)
        }
        ## Constraint on the network size
        
        if(verbose){
            message(sprintf("\tScanning significance threshold at rough stage..."), appendLF=T)
        }
        ## at rough phase
        fdr_rough <- NULL
        nsize_rough <- 0
        
        st <- ceiling(log10(min(pval[pval!=0])))
        if(st < -300){
        	all_i <- c(st, seq(from=-300,to=-250,by=50), seq(from=-200,to=-120,by=20), seq(from=-100,to=-20,by=10), seq(from=-10,to=0,by=1))
        }else if(st < -200){
        	all_i <- c(st, seq(from=-200,to=-120,by=20), seq(from=-100,to=-20,by=10), seq(from=-10,to=0,by=1))
        }else if(st < -100){
        	all_i <- c(st, seq(from=-100,to=-20,by=10), seq(from=-10,to=0,by=1))
        }else{
        	all_i <- seq(from=st,to=0)
        }
        
        for(i in 1:length(all_i)){
            fdr_test <- 10^all_i[i]
            
            if(method!="customised"){
                scores_test <- dBUMscore(fit=fit, method=method, fdr=fdr_test, scatter.bum=F)
            }else if(method=="customised"){
                scores_test <- dFDRscore(pval, fdr.threshold=fdr_test, scatter=F)
                if(is.null(scores_test)){
                	break
                }
            }
    
            module_test <- suppressWarnings(dNetFind(g, scores_test))
            nsize_test <- vcount(module_test)
            
            if(verbose){
                message(sprintf("\t\tsignificance threshold: %1.2e, corresponding to the network size (%d nodes)", fdr_test, nsize_test), appendLF=T)
            }
            
            fdr_rough <- fdr_test
            nsize_rough <- nsize_test
            if(nsize_test >= nsize){    
                break
            }
        }
        
        if(nsize_rough==nsize){
            fdr_final <- fdr_rough
        }else{
            if(verbose){
                message(sprintf("\tScanning significance threshold at finetuning stage..."), appendLF=T)
            }
            ## at finetune phase
            fdr_final <- NULL
            for(fdr_test in seq(from=fdr_rough/10+fdr_rough/20,to=fdr_rough-fdr_rough/20,by=fdr_rough/20)){
            
                if(method!="customised"){
                    scores_test <- dBUMscore(fit=fit, method=method, fdr=fdr_test, scatter.bum=F)
                }else if(method=="customised"){
                    scores_test <- dFDRscore(pval, fdr.threshold=fdr_test, scatter=F)
                }
            
                module_test <- suppressWarnings(dNetFind(g, scores_test))
                nsize_test <- vcount(module_test)
            
                if(verbose){
                    message(sprintf("\t\tsignificance threshold: %1.2e, corresponding to the network size (%d nodes)", fdr_test, nsize_test), appendLF=T)
                }
            	
            	fdr_final <- fdr_test
                if(nsize_test >= nsize){
                    break
                }
            }
        }
        
    }
    
    if(is.null(fdr_final)){
        if(verbose){
            message(sprintf("\tNo significance threshold required"), appendLF=T)
        }
    }else{
        if(verbose){
            message(sprintf("\tsignificance threshold: %1.2e", fdr_final), appendLF=T)
        }
    }
    
    if(method!="customised"){
        if(verbose){
            message(sprintf("Third, calculate the scores according to the fitted BUM and FDR threshold (if any)..."), appendLF=T)
        }
        ## calculate the scores according to the fitted BUM and fdr threshold (fdr_final) 
        if(plot){
            scores <- dBUMscore(fit=fit, method=method, fdr=fdr_final, scatter.bum=T)
        }else{
            scores <- dBUMscore(fit=fit, method=method, fdr=fdr_final, scatter.bum=F)
        }
    }else if(method=="customised"){
        if(verbose){
            message(sprintf("Third, calculate the scores according to the input fdr (or p-value) and the threshold (if any)..."), appendLF=T)
        }
        ## calculate the scores according to the input fdr (or p-value) and the threshold (fdr_final) 
        if(plot){
            scores <- dFDRscore(pval, fdr.threshold=fdr_final, scatter=T)
        }else{
            scores <- dFDRscore(pval, fdr.threshold=fdr_final, scatter=F)
        }
    }
    ####################
    
    if(verbose){
        message(sprintf("\tAmongst %d scores, there are %d positives.", length(scores), sum(scores>0)), appendLF=T)
        message(sprintf("Finally, find the subgraph from the input graph with %d nodes and %d edges...", vcount(g), ecount(g)), appendLF=T)
    }
    ## find the subgraph with the maximum score
    subgraph <- dNetFind(g, scores)
    if(verbose){
        message(sprintf("\tSize of the subgraph: %d nodes and %d edges", vcount(subgraph), ecount(subgraph)), appendLF=T)
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    return(subgraph)

}