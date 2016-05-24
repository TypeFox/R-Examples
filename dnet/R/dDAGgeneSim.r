#' Function to calculate pair-wise semantic similarity between genes based on a direct acyclic graph (DAG) with annotated data
#'
#' \code{dDAGgeneSim} is supposed to calculate pair-wise semantic similarity between genes based on a direct acyclic graph (DAG) with annotated data. It first calculates semantic similarity between terms and then derives semantic similarity between genes from terms-term semantic similarity. Parallel computing is also supported for Linux or Mac operating systems.
#'
#' @param g an object of class "igraph" or "graphNEL". It must contain a vertex attribute called 'annotations' for storing annotation data (see example for howto)
#' @param genes the genes between which pair-wise semantic similarity is calculated. If NULL, all genes annotatable in the input dag will be used for calculation, which is very prohibitively expensive!
#' @param method.gene the method used for how to derive semantic similarity between genes from semantic similarity between terms. It can be "average" for average similarity between any two terms (one from gene 1, the other from gene 2), "max" for the maximum similarity between any two terms, "BM.average" for best-matching (BM) based average similarity (i.e. for each term of either gene, first calculate maximum similarity to any term in the other gene, then take average of maximum similarity; the final BM-based average similiary is the pre-calculated average between two genes in pair), "BM.max" for BM based maximum similarity (i.e. the same as "BM.average", but the final BM-based maximum similiary is the maximum of the pre-calculated average between two genes in pair), "BM.complete" for BM-based complete-linkage similarity (inspired by complete-linkage concept: the least of any maximum similarity between a term of one gene and a term of the other gene). When comparing BM-based similarity between genes, "BM.average" and "BM.max" are sensitive to the number of terms invovled; instead, "BM.complete" is much robust in this aspect. By default, it uses "BM.average".
#' @param method.term the method used to measure semantic similarity between terms. It can be "Resnik" for information content (IC) of most informative common ancestor (MICA) (see \url{http://arxiv.org/pdf/cmp-lg/9511007.pdf}), "Lin" for 2*IC at MICA divided by the sum of IC at pairs of terms (see \url{https://www.cse.iitb.ac.in/~cs626-449/Papers/WordSimilarity/3.pdf}), "Schlicker" for weighted version of 'Lin' by the 1-prob(MICA) (see \url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}), "Jiang" for 1 - difference between the sum of IC at pairs of terms and 2*IC at MICA (see \url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}), "Pesquita" for graph information content similarity related to Tanimoto-Jacard index (ie. summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 (see \url{http://www.ncbi.nlm.nih.gov/pubmed/18460186}))
#' @param force logical to indicate whether the only most specific terms (for each gene) will be used. By default, it sets to true. It is always advisable to use this since it is computationally fast but without compromising accuracy (considering the fact that true-path-rule has been applied when running \code{\link{dDAGannotate}})
#' @param fast logical to indicate whether a vectorised fast computation is used. By default, it sets to true. It is always advisable to use this vectorised fast computation; since the conventional computation is just used for understanding scripts
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a sparse matrix containing pair-wise semantic similarity between input genes. This sparse matrix can be converted to the full matrix via the function \code{as.matrix}
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @import Matrix
#' @seealso \code{\link{dDAGtermSim}}, \code{\link{dDAGinduce}}, \code{\link{dDAGtip}}, \code{\link{dCheckParallel}}
#' @include dDAGgeneSim.r
#' @examples
#' \dontrun{
#' # 1) load HPPA as igraph object
#' ig.HPPA <-dRDataLoader(RData='ig.HPPA')
#' g <- ig.HPPA
#'
#' # 2) load human genes annotated by HPPA
#' org.Hs.egHPPA <- dRDataLoader(RData='org.Hs.egHPPA')
#'
#' # 3) prepare for ontology and its annotation information 
#' dag <- dDAGannotate(g, annotations=org.Hs.egHPPA, path.mode="all_paths", verbose=TRUE)
#'
#' # 4) calculate pair-wise semantic similarity between 5 randomly chosen genes 
#' allgenes <- unique(unlist(V(dag)$annotations))
#' genes <- sample(allgenes,5)
#' sim <- dDAGgeneSim(g=dag, genes=genes, method.gene="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
#' sim
#' }

dDAGgeneSim <- function (g, genes=NULL, method.gene=c("BM.average","BM.max","BM.complete","average","max"), method.term=c("Resnik","Lin","Schlicker","Jiang","Pesquita"), force=TRUE, fast=TRUE, parallel=TRUE, multicores=NULL, verbose=TRUE)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################

    method.gene <- match.arg(method.gene)
    method.term <- match.arg(method.term)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    if(is.null(V(ig)$annotations) | is.null(V(ig)$IC)){
        stop("The function requires that input graph has already contained annotation data. Please first run 'dDAGannotate'.\n")
    }

    ####################################################
    ## A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d (%s)", i, B, as.character(Sys.time())), appendLF=T)
            }
        }
    }
    
    ####################################################

    if(verbose){
        message(sprintf("First, extract all annotatable genes (%s)...", as.character(Sys.time())), appendLF=T)
    }

    anno <- V(ig)$annotations
    allgenes <- sort(as.numeric(unique(unlist(anno))))
    
    ## checking input genes
    genes <- genes[!is.na(genes)]
    if(is.null(genes) || is.na(genes)){
        genes <- allgenes
    }else{
        flag <- genes %in% allgenes
        if(sum(flag)!=0){
            genes <- genes[flag]
        }else{
            genes <- allgenes
        }
    }
    
    ## pre-compute a sparse matrix of input genes x terms
    allterms <- 1:length(anno)
    sGT <- Matrix::Matrix(0, nrow=length(genes), ncol=length(allterms), sparse=T)
    for(j in 1:length(allterms)){
        ind <- match(anno[[j]], genes)
        flag <- ind[!is.na(ind)]
        if(length(flag)!=0){
            sGT[flag,j] <- 1
        }
    }
    colnames(sGT) <- V(ig)$name
    rownames(sGT) <- genes
    
    if(verbose){
        message(sprintf("\tthere are %d input genes amongst %d annotatable genes", length(genes), length(allgenes)), appendLF=T)
    }
    
    ## a list of genes, each containing terms annotated by
    genes2terms <- sapply(1:length(genes), function(x){
        res <- names(which(sGT[x,]==1))
        if(force){
            subg <- dDAGinduce(ig, nodes_query=res, path.mode="all_paths")
            res <- dDAGtip(subg)
        }
        return(res)
    })
    names(genes2terms) <- genes
    terms <- unique(unlist(genes2terms))
    
    ## also instore index for terms (in genes2terms)
    genes2terms_index <- sapply(genes2terms, function(x){
        match(x, terms)
    })
    
    if(verbose){
        if(force){
            message(sprintf("Second, pre-compute semantic similarity between %d terms (forced to be the most specific for each gene) using %s method (%s)...", length(terms), method.term, as.character(Sys.time())), appendLF=T)
        }else{
            message(sprintf("Second, pre-compute semantic similarity between %d terms using %s method (%s)...", length(terms), method.term, as.character(Sys.time())), appendLF=T)
        }
    }
    ## pre-compute semantic similarity between terms in subject
    sim.term <- suppressMessages(dDAGtermSim(ig, terms=terms, method=method.term, parallel=parallel, multicores=multicores, verbose=T))
    
    if(verbose){
        message(sprintf("Last, calculate pair-wise semantic similarity between %d genes using %s method (%s)...", length(genes), method.gene, as.character(Sys.time())), appendLF=T)
    }
    num_genes <- length(genes2terms)
    
    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
        flag_parallel <- dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            if(method.gene=='average'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=T, .combine=rbind), {
                    ind1 <- genes2terms_index[[i]]
                    progress_indicate(i, num_genes, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_genes
                        ind_js <- genes2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            mean(sim12[,which(new_ind_js==k)])
                        })
                        x <- rep(0, num_genes)
                        x[js] <- res
                        x
                    }
                })
            }else if(method.gene=='max'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=T, .combine=rbind), {
                    ind1 <- genes2terms_index[[i]]
                    progress_indicate(i, num_genes, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_genes
                        ind_js <- genes2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            max(sim12[,which(new_ind_js==k)])
                        })
                        x <- rep(0, num_genes)
                        x[js] <- res
                        x
                    }
                })
            }else if(method.gene=='BM.average'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=T, .combine=rbind), {
                    ind1 <- genes2terms_index[[i]]
                    progress_indicate(i, num_genes, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_genes
                        ind_js <- genes2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            x <- as.matrix(sim12[,which(new_ind_js==k)])
                            0.5*(mean(apply(x,1,max)) + mean(apply(x,2,max)))
                        })
                        x <- rep(0, num_genes)
                        x[js] <- res
                        x
                    }
                })
            }else if(method.gene=='BM.max'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=T, .combine=rbind), {
                    ind1 <- genes2terms_index[[i]]
                    progress_indicate(i, num_genes, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_genes
                        ind_js <- genes2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            x <- as.matrix(sim12[,which(new_ind_js==k)])
                            max(mean(apply(x,1,max)), mean(apply(x,2,max)))
                        })
                        x <- rep(0, num_genes)
                        x[js] <- res
                        x
                    }
                })
            }else if(method.gene=='BM.complete'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=T, .combine=rbind), {
                    ind1 <- genes2terms_index[[i]]
                    progress_indicate(i, num_genes, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_genes
                        ind_js <- genes2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            x <- as.matrix(sim12[,which(new_ind_js==k)])
                            min(c(apply(x,1,max),apply(x,2,max)))
                        })
                        x <- rep(0, num_genes)
                        x[js] <- res
                        x
                    }
                })
            }

            ## add the last row
            sim <- rbind(sim, rep(0, num_genes))

            sim <- sim + Matrix::t(sim)
            sim <- Matrix::Matrix(sim, sparse=T)
        }
    }
    
    ###### non-parallel computing
    if(flag_parallel==F){
        ## calculate pair-wise semantic similarity between input genes
        sim <- Matrix::Matrix(0, nrow=length(genes), ncol=length(genes), sparse=T)
    
        ## print with possibly greater accuracy:
        ##op <- options(digits.secs = 6)
        ##options(op)
     
        if(method.gene=='average'){
            for(i in 1:(num_genes-1)){
                ind1 <- genes2terms_index[[i]]
                progress_indicate(i, num_genes, 10, flag=T)
                if(fast){
                    js <- (i+1):num_genes
                    ind_js <- genes2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        mean(sim12[,which(new_ind_js==k)])
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_genes){
                        ind2 <- genes2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- mean(sim12)
                    }
                }
            }
        }else if(method.gene=='max'){
            for(i in 1:(num_genes-1)){
                ind1 <- genes2terms_index[[i]]
                progress_indicate(i, num_genes, 10, flag=T)
                if(fast){
                    js <- (i+1):num_genes
                    ind_js <- genes2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        max(sim12[,which(new_ind_js==k)])
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_genes){
                        ind2 <- genes2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- max(sim12)
                    }
                }
            }
        }else if(method.gene=='BM.average'){
            for(i in 1:(num_genes-1)){
                ind1 <- genes2terms_index[[i]]
                progress_indicate(i, num_genes, 10, flag=T)
                if(fast){
                    js <- (i+1):num_genes
                    ind_js <- genes2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        x <- as.matrix(sim12[,which(new_ind_js==k)])
                        0.5*(mean(apply(x,1,max)) + mean(apply(x,2,max)))
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_genes){
                        ind2 <- genes2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- 0.5*(mean(apply(sim12,1,max)) + mean(apply(sim12,2,max)))
                    }
                }
            }
        
        }else if(method.gene=='BM.max'){
            for(i in 1:(num_genes-1)){
                ind1 <- genes2terms_index[[i]]
                progress_indicate(i, num_genes, 10, flag=T)
                if(fast){
                    js <- (i+1):num_genes
                    ind_js <- genes2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        x <- as.matrix(sim12[,which(new_ind_js==k)])
                        max(mean(apply(x,1,max)), mean(apply(x,2,max)))
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_genes){
                        ind2 <- genes2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- max(mean(apply(sim12,1,max)), mean(apply(sim12,2,max)))
                    }
                }
            }
        }else if(method.gene=='BM.complete'){
            for(i in 1:(num_genes-1)){
                ind1 <- genes2terms_index[[i]]
                progress_indicate(i, num_genes, 10, flag=T)
                if(fast){
                    js <- (i+1):num_genes
                    ind_js <- genes2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        x <- as.matrix(sim12[,which(new_ind_js==k)])
                        min(c(apply(x,1,max),apply(x,2,max)))
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_genes){
                        ind2 <- genes2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- min(c(apply(sim12,1,max),apply(sim12,2,max)))
                    }
                }
            }
        }
        sim <- sim + Matrix::t(sim)
    
    }
    
    rownames(sim) <- colnames(sim) <- genes
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(sim)
}
