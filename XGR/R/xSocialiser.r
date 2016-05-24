#' Function to calculate pair-wise semantic similarity given the input data and the ontology and its annotation
#'
#' \code{xSocialiser} is supposed to calculate pair-wise semantic similarity given the input data and the ontology direct acyclic graph (DAG) and its annotation. It returns an object of class "igraph", a network representation of socialized genes/SNPs. It first calculates semantic similarity between terms and then derives semantic similarity from term-term semantic similarity. Parallel computing is also supported for Linux or Mac operating systems.
#'
#' @param data an input vector containing a list of genes or SNPs of interest between which pair-wise semantic similarity is calculated/socialized
#' @param annotation the vertices/nodes for which annotation data are provided. It can be a sparse Matrix of class "dgCMatrix" (with variants/genes as rows and terms as columns), or a list of nodes/terms each containing annotation data, or an object of class 'GS' (basically a list for each node/term with annotation data)
#' @param g an object of class "igraph" to represent DAG. It must have node/vertice attributes: "name" (i.e. "Term ID"), "term_id" (i.e. "Term ID"), "term_name" (i.e "Term Name") and "term_distance" (i.e. Term Distance: the distance to the root; always 0 for the root itself)
#' @param measure the measure used to derive semantic similarity between genes/SNPs from semantic similarity between terms. Take the semantic similartity between SNPs as an example. It can be "average" for average similarity between any two terms (one from SNP 1, the other from SNP 2), "max" for the maximum similarity between any two terms, "BM.average" for best-matching (BM) based average similarity (i.e. for each term of either SNP, first calculate maximum similarity to any term in the other SNP, then take average of maximum similarity; the final BM-based average similiary is the pre-calculated average between two SNPs in pair), "BM.max" for BM based maximum similarity (i.e. the same as "BM.average", but the final BM-based maximum similiary is the maximum of the pre-calculated average between two SNPs in pair), "BM.complete" for BM-based complete-linkage similarity (inspired by complete-linkage concept: the least of any maximum similarity between a term of one SNP and a term of the other SNP). When comparing BM-based similarity between SNPs, "BM.average" and "BM.max" are sensitive to the number of terms invovled; instead, "BM.complete" is much robust in this aspect. By default, it uses "BM.average"
#' @param method.term the method used to measure semantic similarity between terms. It can be "Resnik" for information content (IC) of most informative common ancestor (MICA) (see \url{http://dl.acm.org/citation.cfm?id=1625914}), "Lin" for 2*IC at MICA divided by the sum of IC at pairs of terms (see \url{https://www.cse.iitb.ac.in/~cs626-449/Papers/WordSimilarity/3.pdf}), "Schlicker" for weighted version of 'Lin' by the 1-prob(MICA) (see \url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}), "Jiang" for 1 - difference between the sum of IC at pairs of terms and 2*IC at MICA (see \url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}), "Pesquita" for graph information content similarity related to Tanimoto-Jacard index (ie. summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 (see \url{http://www.ncbi.nlm.nih.gov/pubmed/18460186}))
#' @param rescale logical to indicate whether the resulting values are rescaled to the range [0,1]. By default, it sets to true
#' @param force logical to indicate whether the only most specific terms (for each SNP) will be used. By default, it sets to true. It is always advisable to use this since it is computationally fast but without compromising accuracy (considering the fact that true-path-rule has been applied when running \code{\link{xDAGanno}})
#' @param fast logical to indicate whether a vectorised fast computation is used. By default, it sets to true. It is always advisable to use this vectorised fast computation; since the conventional computation is just used for understanding scripts
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to true
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' It returns an object of class "igraph", with nodes for input genes/SNPs and edges for pair-wise semantic similarity between them. If no similarity is calculuated, it returns NULL.
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @import Matrix
#' @seealso \code{\link{xDAGsim}}, \code{\link{xSocialiserGenes}}, \code{\link{xSocialiserSNPs}}
#' @include xSocialiser.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' 
#' # 1) SNP-based enrichment analysis using GWAS Catalog traits (mapped to EF)
#' # 1a) ig.EF (an object of class "igraph" storing as a directed graph)
#' g <- xRDataLoader('ig.EF')
#' g
#'
#' # 1b) load GWAS SNPs annotated by EF (an object of class "dgCMatrix" storing a spare matrix)
#' anno <- xRDataLoader(RData='GWAS2EF')
#' 
#' # 1c) prepare the input SNPs of interest (eg 8 randomly chosen SNPs)
#' allSNPs <- rownames(anno)
#' data <- sample(allSNPs,8)
#' 
#' # 1d) perform calculate pair-wise semantic similarity between 8 randomly chosen SNPs
#' sim <- xSocialiser(data=data, annotation=anno, g=g, parallel=FALSE, verbose=TRUE)
#' sim
#'
#' # 1e) save similarity results to the file called 'EF_similarity.txt'
#' output <- igraph::get.data.frame(sim, what="edges")
#' utils::write.table(output, file="EF_similarity.txt", sep="\t", row.names=FALSE)
#'
#' # 1f) visualise the SNP network
#' ## extract edge weight (with 2-digit precision)
#' x <- signif(as.numeric(E(sim)$weight), digits=2)
#' ## rescale into an interval [1,4] as edge width
#' edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
#' ## do visualisation
#' xVisNet(g=sim, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
#' }

xSocialiser <- function(data, annotation, g, measure=c("BM.average","BM.max","BM.complete","average","max"), method.term=c("Resnik","Lin","Schlicker","Jiang","Pesquita"), rescale=TRUE, force=TRUE, fast=TRUE, parallel=TRUE, multicores=NULL, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=TRUE, verbose=T)
{

    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    measure <- match.arg(measure)
    method.term <- match.arg(method.term)
    path.mode <- match.arg(path.mode)
    
    
    if(class(annotation)=="GS"){
        originAnnos <- annotation$gs
    }else if(class(annotation)=="list"){
        originAnnos <- annotation
    }else if(class(annotation)=="dgCMatrix"){
		D <- annotation
		originAnnos <- sapply(1:ncol(D), function(j){
			names(which(D[,j]!=0))
		})
		names(originAnnos) <- colnames(annotation)
    }else{
    	stop("The input annotation must be either 'GS' or 'list' or 'dgCMatrix' object.\n")
    }
    annotation <- originAnnos

    ig <- g
    if (class(ig) != "igraph"){
        stop("The function must apply to the 'igraph' object.\n")
    }else{
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("First, generate a subgraph induced (via '%s' mode) by the annotation data (%s) ...", path.mode, as.character(now)), appendLF=T)
		}
    	
    	# obtain the induced subgraph according to the input annotation data based on shortest paths (i.e. the most concise subgraph induced)
		ig <- xDAGanno(g=ig, annotation=annotation, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose)
	
	}
	
    anno <- V(ig)$anno
    allSNPs <- sort(unique(unlist(anno)))
    
    ## checking input SNPs
    SNPs <- data[!is.na(data)]
    if(is.null(SNPs) || is.na(SNPs)){
        #SNPs <- allSNPs
    }else{
        flag <- SNPs %in% allSNPs
        if(sum(flag)!=0){
            SNPs <- SNPs[flag]
        }else{
        	SNPs <- SNPs[flag]
            #SNPs <- allSNPs
        }
    }
    
    ## less than 2 SNPs/Genes are annotable
    if(length(SNPs)<=1 || is.na(SNPs)){
		res <- NULL
    	warning("The function returns NULL as no similarity is found at all.\n")
    	return(res)
    }
    
    ## pre-compute a sparse matrix of input SNPs x terms
    allterms <- 1:length(anno)
    sGT <- Matrix::Matrix(0, nrow=length(SNPs), ncol=length(allterms), sparse=T)
    for(j in 1:length(allterms)){
        ind <- match(anno[[j]], SNPs)
        flag <- ind[!is.na(ind)]
        if(length(flag)!=0){
            sGT[flag,j] <- 1
        }
    }
    colnames(sGT) <- V(ig)$name
    rownames(sGT) <- SNPs
    
    if(verbose){
        message(sprintf("\tthere are %d inputs amongst %d annotatable", length(SNPs), length(allSNPs)), appendLF=T)
    }
    
    ## a list of SNPs, each containing terms annotated by
    SNPs2terms <- sapply(1:length(SNPs), function(x){
        res <- names(which(sGT[x,]==1))
        if(force){
            subg <- dnet::dDAGinduce(ig, nodes_query=res, path.mode="all_paths")
            res <- dnet::dDAGtip(subg)
        }
        return(res)
    })
    names(SNPs2terms) <- SNPs
    terms <- unique(unlist(SNPs2terms))
    
    ## also instore index for terms (in SNPs2terms)
    SNPs2terms_index <- sapply(SNPs2terms, function(x){
        match(x, terms)
    })
    
    ##############################################################################################
    ## A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d (%s)", i, B, as.character(Sys.time())), appendLF=T)
            }
        }
    }
    ##############################################################################################

    if(verbose){
    	now <- Sys.time()
        if(force){
            message(sprintf("Next, pre-compute semantic similarity between %d terms (forced to be the most specific for each gene/SNP) using '%s' method (%s)...", length(terms), method.term, as.character(now)), appendLF=T)
        }else{
            message(sprintf("Next, pre-compute semantic similarity between %d terms using '%s' method (%s)...", length(terms), method.term, as.character(now)), appendLF=T)
        }
    }

    ## pre-compute semantic similarity between terms in subject
    sim.term <- suppressMessages(xDAGsim(ig, terms=terms, method.term=method.term, fast=fast, parallel=parallel, multicores=multicores, verbose=T))
    
    if (class(sim.term) == "igraph"){
    	sim.term <- xConverter(sim.term, from='igraph', to='dgCMatrix')
    }
    
    ##############################################################################################
    
    if(verbose){
    	now <- Sys.time()
        message(sprintf("Last, calculate pair-wise semantic similarity between %d genes/SNPs using '%s' measure (%s)...", length(SNPs), measure, as.character(now)), appendLF=T)
    }
    num_SNPs <- length(SNPs2terms)
	
    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
    
        flag_parallel <- dnet::dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            if(measure=='average'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_SNPs-1), .inorder=T, .combine=rbind), {
                    ind1 <- SNPs2terms_index[[i]]
                    progress_indicate(i, num_SNPs, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_SNPs
                        ind_js <- SNPs2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            mean(sim12[,which(new_ind_js==k)])
                        })
                        x <- rep(0, num_SNPs)
                        x[js] <- res
                        x
                    }
                })
            }else if(measure=='max'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_SNPs-1), .inorder=T, .combine=rbind), {
                    ind1 <- SNPs2terms_index[[i]]
                    progress_indicate(i, num_SNPs, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_SNPs
                        ind_js <- SNPs2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            max(sim12[,which(new_ind_js==k)])
                        })
                        x <- rep(0, num_SNPs)
                        x[js] <- res
                        x
                    }
                })
            }else if(measure=='BM.average'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_SNPs-1), .inorder=T, .combine=rbind), {
                    ind1 <- SNPs2terms_index[[i]]
                    progress_indicate(i, num_SNPs, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_SNPs
                        ind_js <- SNPs2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            x <- as.matrix(sim12[,which(new_ind_js==k)])
                            0.5*(mean(apply(x,1,max)) + mean(apply(x,2,max)))
                        })
                        x <- rep(0, num_SNPs)
                        x[js] <- res
                        x
                    }
                })
            }else if(measure=='BM.max'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_SNPs-1), .inorder=T, .combine=rbind), {
                    ind1 <- SNPs2terms_index[[i]]
                    progress_indicate(i, num_SNPs, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_SNPs
                        ind_js <- SNPs2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            x <- as.matrix(sim12[,which(new_ind_js==k)])
                            max(mean(apply(x,1,max)), mean(apply(x,2,max)))
                        })
                        x <- rep(0, num_SNPs)
                        x[js] <- res
                        x
                    }
                })
            }else if(measure=='BM.complete'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_SNPs-1), .inorder=T, .combine=rbind), {
                    ind1 <- SNPs2terms_index[[i]]
                    progress_indicate(i, num_SNPs, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_SNPs
                        ind_js <- SNPs2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            x <- as.matrix(sim12[,which(new_ind_js==k)])
                            min(c(apply(x,1,max),apply(x,2,max)))
                        })
                        x <- rep(0, num_SNPs)
                        x[js] <- res
                        x
                    }
                })
            }

            ## add the last row
            sim <- rbind(sim, rep(0, num_SNPs))

            sim <- sim + Matrix::t(sim)
            sim <- Matrix::Matrix(sim, sparse=T)
        }
    }
    
    ###### non-parallel computing
    if(flag_parallel==F){
        ## calculate pair-wise semantic similarity between input SNPs
        sim <- Matrix::Matrix(0, nrow=length(SNPs), ncol=length(SNPs), sparse=T)
    
        ## print with possibly greater accuracy:
        ##op <- options(digits.secs = 6)
        ##options(op)
     
        if(measure=='average'){
            for(i in 1:(num_SNPs-1)){
                ind1 <- SNPs2terms_index[[i]]
                progress_indicate(i, num_SNPs, 10, flag=T)
                if(fast){
                    js <- (i+1):num_SNPs
                    ind_js <- SNPs2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        mean(sim12[,which(new_ind_js==k)])
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_SNPs){
                        ind2 <- SNPs2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- mean(sim12)
                    }
                }
            }
        }else if(measure=='max'){
            for(i in 1:(num_SNPs-1)){
                ind1 <- SNPs2terms_index[[i]]
                progress_indicate(i, num_SNPs, 10, flag=T)
                if(fast){
                    js <- (i+1):num_SNPs
                    ind_js <- SNPs2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        max(sim12[,which(new_ind_js==k)])
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_SNPs){
                        ind2 <- SNPs2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- max(sim12)
                    }
                }
            }
        }else if(measure=='BM.average'){
            for(i in 1:(num_SNPs-1)){
                ind1 <- SNPs2terms_index[[i]]
                progress_indicate(i, num_SNPs, 10, flag=T)
                if(fast){
                    js <- (i+1):num_SNPs
                    ind_js <- SNPs2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        x <- as.matrix(sim12[,which(new_ind_js==k)])
                        0.5*(mean(apply(x,1,max)) + mean(apply(x,2,max)))
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_SNPs){
                        ind2 <- SNPs2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- 0.5*(mean(apply(sim12,1,max)) + mean(apply(sim12,2,max)))
                    }
                }
            }
        
        }else if(measure=='BM.max'){
            for(i in 1:(num_SNPs-1)){
                ind1 <- SNPs2terms_index[[i]]
                progress_indicate(i, num_SNPs, 10, flag=T)
                if(fast){
                    js <- (i+1):num_SNPs
                    ind_js <- SNPs2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        x <- as.matrix(sim12[,which(new_ind_js==k)])
                        max(mean(apply(x,1,max)), mean(apply(x,2,max)))
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_SNPs){
                        ind2 <- SNPs2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- max(mean(apply(sim12,1,max)), mean(apply(sim12,2,max)))
                    }
                }
            }
        }else if(measure=='BM.complete'){
            for(i in 1:(num_SNPs-1)){
                ind1 <- SNPs2terms_index[[i]]
                progress_indicate(i, num_SNPs, 10, flag=T)
                if(fast){
                    js <- (i+1):num_SNPs
                    ind_js <- SNPs2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        x <- as.matrix(sim12[,which(new_ind_js==k)])
                        min(c(apply(x,1,max),apply(x,2,max)))
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_SNPs){
                        ind2 <- SNPs2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- min(c(apply(sim12,1,max),apply(sim12,2,max)))
                    }
                }
            }
        }
        sim <- sim + Matrix::t(sim)
    
    }
    rownames(sim) <- colnames(sim) <- SNPs
    sim[as.matrix(is.na(sim))] <- 0
    #sim <- signif(sim, digits=2)
	
	if(rescale){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Also rescale similarity into the [0,1] range (%s)", nrow(sim), as.character(now)), appendLF=T)
		}
	
		# rescale to [0 1]
		d_min <- min(sim)
		d_max <- max(sim)
		sim <- apply(sim, 2, function(x){
			(x - d_min)/(d_max - d_min)
		})
		sim <- Matrix::Matrix(sim, sparse=T)
	}
	
    ####################################################################################
    
    if (class(sim) == "dgCMatrix" | class(sim) == "dsCMatrix"){
    	res <- xConverter(sim, from="dgCMatrix", to="igraph", verbose=F)
    }
    
    ## no edges
    if(ecount(res)==0){
    	res <- NULL
    	warning("The function returns NULL as no similarity is found at all.\n")
    }
    
    invisible(res)
}
