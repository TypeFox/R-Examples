#' Function to calculate pair-wise semantic similarity between domains based on a direct acyclic graph (DAG) with annotated data
#'
#' \code{dcDAGdomainSim} is supposed to calculate pair-wise semantic similarity between domains based on a direct acyclic graph (DAG) with annotated data. It first calculates semantic similarity between terms and then derives semantic similarity between domains from terms-term semantic similarity. Parallel computing is also supported for Linux or Mac operating systems.
#'
#' @param g an object of class "igraph" or \code{\link{Onto}}. It must contain a node attribute called 'annotations' for storing annotation data (see example for howto)
#' @param domains the domains between which pair-wise semantic similarity is calculated. If NULL, all domains annotatable in the input dag will be used for calcluation, which is very prohibitively expensive!
#' @param method.domain the method used for how to derive semantic similarity between domains from semantic similarity between terms. It can be "average" for average similarity between any two terms (one from domain 1, the other from domain 2), "max" for the maximum similarity between any two terms, "BM.average" for best-matching (BM) based average similarity (i.e. for each term of either domain, first calculate maximum similarity to any term in the other domain, then take average of maximum similarity; the final BM-based average similiary is the pre-calculated average between two domains in pair), "BM.max" for BM based maximum similarity (i.e. the same as "BM.average", but the final BM-based maximum similiary is the maximum of the pre-calculated average between two domains in pair), "BM.complete" for BM-based complete-linkage similarity (inspired by complete-linkage concept: the least of any maximum similarity between a term of one domain and a term of the other domain). When comparing BM-based similarity between domains, "BM.average" and "BM.max" are sensitive to the number of terms invovled; instead, "BM.complete" is much robust in this aspect. By default, it uses "BM.average".
#' @param method.term the method used to measure semantic similarity between terms. It can be "Resnik" for information content (IC) of most informative common ancestor (MICA) (see \url{http://arxiv.org/pdf/cmp-lg/9511007.pdf}), "Lin" for 2*IC at MICA divided by the sum of IC at pairs of terms (see \url{http://webdocs.cs.ualberta.ca/~lindek/papers/sim.pdf}), "Schlicker" for weighted version of 'Lin' by the 1-prob(MICA) (see \url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}), "Jiang" for 1 - difference between the sum of IC at pairs of terms and 2*IC at MICA (see \url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}), "Pesquita" for graph information content similarity related to Tanimoto-Jacard index (ie. summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 (see \url{http://www.ncbi.nlm.nih.gov/pubmed/18460186}))
#' @param force logical to indicate whether the only most specific terms (for each domain) will be used. By default, it sets to true. It is always advisable to use this since it is computationally fast but without compromising accuracy (considering the fact that true-path-rule has been applied when running \code{\link{dcDAGannotate}})
#' @param fast logical to indicate whether a vectorised fast computation is used. By default, it sets to true. It is always advisable to use this vectorised fast computation; since the conventional computation is just used for understanding scripts
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' an object of S4 class \code{\link{Dnetwork}}. It is a weighted and undirect graph, with following slots:
#' \itemize{
#'  \item{\code{nodeInfo}: an object of S4 class, describing information on nodes/domains}
#'  \item{\code{adjMatrix}: an object of S4 class \code{\link{AdjData}}, containing symmetric adjacency data matrix for pair-wise semantic similarity between domains}
#' }
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @importFrom dnet dDAGinduce dDAGtip dDAGtermSim dCheckParallel visNet
#' @import Matrix
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcDAGannotate}}, \code{\link{dcConverter}}, \code{\link{Dnetwork-class}}
#' @include dcDAGdomainSim.r
#' @examples
#' \dontrun{
#' # 1) Semantic similarity between SCOP domain superfamilies (sf)
#' ## 1a) load onto.GOMF (as 'Onto' object)
#' g <- dcRDataLoader('onto.GOMF')
#' ## 1b) load SCOP superfamilies annotated by GOMF (as 'Anno' object)
#' Anno <- dcRDataLoader('SCOP.sf2GOMF')
#' ## 1c) prepare for ontology appended with annotation information
#' dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths", verbose=FALSE)
#' ## 1d) calculate pair-wise semantic similarity between 8 randomly chosen domains
#' alldomains <- unique(unlist(nInfo(dag)$annotations))
#' domains <- sample(alldomains,8)
#' dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
#' dnetwork
#' ## 1e) convert it to an object of class 'igraph'
#' ig <- dcConverter(dnetwork, from='Dnetwork', to='igraph')
#' ig
#' ## 1f) visualise the domain network
#' ### extract edge weight (with 2-digit precision)
#' x <- signif(E(ig)$weight, digits=2)
#' ### rescale into an interval [1,4] as edge width
#' edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
#' ### do visualisation
#' dnet::visNet(g=ig, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
#'
#' ###########################################################
#' # 2) Semantic similarity between Pfam domains (Pfam)
#' ## 2a) load onto.GOMF (as 'Onto' object)
#' g <- dcRDataLoader('onto.GOMF')
#' ## 2b) load Pfam domains annotated by GOMF (as 'Anno' object)
#' Anno <- dcRDataLoader('Pfam2GOMF')
#' ## 2c) prepare for ontology appended with annotation information
#' dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths", verbose=FALSE)
#' ## 2d) calculate pair-wise semantic similarity between 8 randomly chosen domains
#' alldomains <- unique(unlist(nInfo(dag)$annotations))
#' domains <- sample(alldomains,8)
#' dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
#' dnetwork
#' ## 2e) convert it to an object of class 'igraph'
#' ig <- dcConverter(dnetwork, from='Dnetwork', to='igraph')
#' ig
#' ## 2f) visualise the domain network
#' ### extract edge weight (with 2-digit precision)
#' x <- signif(E(ig)$weight, digits=2)
#' ### rescale into an interval [1,4] as edge width
#' edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
#' ### do visualisation
#' dnet::visNet(g=ig, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
#'
#' ###########################################################
#' # 3) Semantic similarity between InterPro domains (InterPro)
#' ## 3a) load onto.GOMF (as 'Onto' object)
#' g <- dcRDataLoader('onto.GOMF')
#' ## 3b) load InterPro domains annotated by GOMF (as 'Anno' object)
#' Anno <- dcRDataLoader('InterPro2GOMF')
#' ## 3c) prepare for ontology appended with annotation information
#' dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths", verbose=FALSE)
#' ## 3d) calculate pair-wise semantic similarity between 8 randomly chosen domains
#' alldomains <- unique(unlist(nInfo(dag)$annotations))
#' domains <- sample(alldomains,8)
#' dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
#' dnetwork
#' ## 3e) convert it to an object of class 'igraph'
#' ig <- dcConverter(dnetwork, from='Dnetwork', to='igraph')
#' ig
#' ## 3f) visualise the domain network
#' ### extract edge weight (with 2-digit precision)
#' x <- signif(E(ig)$weight, digits=2)
#' ### rescale into an interval [1,4] as edge width
#' edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
#' ### do visualisation
#' dnet::visNet(g=ig, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
#' 
#' ###########################################################
#' # 4) Semantic similarity between Rfam RNA families (Rfam)
#' ## 4a) load onto.GOBP (as 'Onto' object)
#' g <- dcRDataLoader('onto.GOBP')
#' ## 4b) load Rfam families annotated by GOBP (as 'Anno' object)
#' Anno <- dcRDataLoader('Rfam2GOBP')
#' ## 4c) prepare for ontology appended with annotation information
#' dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths", verbose=FALSE)
#' ## 4d) calculate pair-wise semantic similarity between 8 randomly chosen RNAs
#' alldomains <- unique(unlist(nInfo(dag)$annotations))
#' domains <- sample(alldomains,8)
#' dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
#' dnetwork
#' ## 4e) convert it to an object of class 'igraph'
#' ig <- dcConverter(dnetwork, from='Dnetwork', to='igraph')
#' ig
#' ## 4f) visualise the domain network
#' ### extract edge weight (with 2-digit precision)
#' x <- signif(E(ig)$weight, digits=2)
#' ### rescale into an interval [1,4] as edge width
#' edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
#' ### do visualisation
#' dnet::visNet(g=ig, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
#'
#' ###########################################################
#' # 5) Advanced usage: customised data for ontology and annotations
#' # 5a) customise ontology
#' g <- dcBuildOnto(relations.file="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_edges.txt", nodes.file="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_nodes.txt", output.file="ontology.RData")
#' # 5b) customise Anno
#' Anno <- dcBuildAnno(domain_info.file="http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt", term_info.file="http://dcgor.r-forge.r-project.org/data/InterPro/GO.txt", association.file="http://dcgor.r-forge.r-project.org/data/InterPro/Domain2GOMF.txt", output.file="annotations.RData")
#' ## 5c) prepare for ontology appended with annotation information
#' dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths", verbose=FALSE)
#' ## 5d) calculate pair-wise semantic similarity between 8 randomly chosen domains
#' alldomains <- unique(unlist(nInfo(dag)$annotations))
#' domains <- sample(alldomains,8)
#' dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
#' dnetwork
#' ## 5e) convert it to an object of class 'igraph'
#' ig <- dcConverter(dnetwork, from='Dnetwork', to='igraph')
#' ig
#' ## 5f) visualise the domain network
#' ### extract edge weight (with 2-digit precision)
#' x <- signif(E(ig)$weight, digits=2)
#' ### rescale into an interval [1,4] as edge width
#' edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
#' ### do visualisation
#' dnet::visNet(g=ig, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
#' }

dcDAGdomainSim <- function (g, domains=NULL, method.domain=c("BM.average","BM.max","BM.complete","average","max"), method.term=c("Resnik","Lin","Schlicker","Jiang","Pesquita"), force=TRUE, fast=TRUE, parallel=TRUE, multicores=NULL, verbose=TRUE)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################

    method.domain <- match.arg(method.domain)
    method.term <- match.arg(method.term)
    
    if(class(g)=="Onto"){
        ig <- dcConverter(g, from='Onto', to='igraph', verbose=F)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'Onto' object.\n")
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
        message(sprintf("First, extract all annotatable domains (%s)...", as.character(Sys.time())), appendLF=T)
    }

    anno <- V(ig)$annotations
    #alldomains <- sort(as.numeric(unique(unlist(anno))))
    alldomains <- sort(unique(unlist(anno)))
    
    ## checking input domains
    domains <- domains[!is.na(domains)]
    if(is.null(domains) || is.na(domains)){
        domains <- alldomains
    }else{
        flag <- domains %in% alldomains
        if(sum(flag)!=0){
            domains <- domains[flag]
        }else{
            domains <- alldomains
        }
    }
    
    ## pre-compute a sparse matrix of input domains x terms
    allterms <- 1:length(anno)
    sGT <- Matrix::Matrix(0, nrow=length(domains), ncol=length(allterms), sparse=T)
    for(j in 1:length(allterms)){
        ind <- match(anno[[j]], domains)
        flag <- ind[!is.na(ind)]
        if(length(flag)!=0){
            sGT[flag,j] <- 1
        }
    }
    colnames(sGT) <- V(ig)$name
    rownames(sGT) <- domains
    
    if(verbose){
        message(sprintf("\tthere are %d input domains amongst %d annotatable domains", length(domains), length(alldomains)), appendLF=T)
    }
    
    ## a list of domains, each containing terms annotated by
    domains2terms <- sapply(1:length(domains), function(x){
        res <- names(which(sGT[x,]==1))
        if(force){
            subg <- dnet::dDAGinduce(ig, nodes_query=res, path.mode="all_paths")
            res <- dnet::dDAGtip(subg)
        }
        return(res)
    })
    names(domains2terms) <- domains
    terms <- unique(unlist(domains2terms))
    
    ## also instore index for terms (in domains2terms)
    domains2terms_index <- sapply(domains2terms, function(x){
        match(x, terms)
    })
    
    if(verbose){
        if(force){
            message(sprintf("Second, pre-compute semantic similarity between %d terms (forced to be the most specific for each domain) using %s method (%s)...", length(terms), method.term, as.character(Sys.time())), appendLF=T)
        }else{
            message(sprintf("Second, pre-compute semantic similarity between %d terms using %s method (%s)...", length(terms), method.term, as.character(Sys.time())), appendLF=T)
        }
    }
    ## pre-compute semantic similarity between terms in subject
    sim.term <- suppressMessages(dnet::dDAGtermSim(ig, terms=terms, method=method.term, parallel=parallel, multicores=multicores, verbose=T))
    
    if(verbose){
        message(sprintf("Last, calculate pair-wise semantic similarity between %d domains using %s method (%s)...", length(domains), method.domain, as.character(Sys.time())), appendLF=T)
    }
    num_domains <- length(domains2terms)
    
    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
    
        flag_parallel <- dnet::dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            if(method.domain=='average'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_domains-1), .inorder=T, .combine=rbind), {
                    ind1 <- domains2terms_index[[i]]
                    progress_indicate(i, num_domains, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_domains
                        ind_js <- domains2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            mean(sim12[,which(new_ind_js==k)])
                        })
                        x <- rep(0, num_domains)
                        x[js] <- res
                        x
                    }
                })
            }else if(method.domain=='max'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_domains-1), .inorder=T, .combine=rbind), {
                    ind1 <- domains2terms_index[[i]]
                    progress_indicate(i, num_domains, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_domains
                        ind_js <- domains2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            max(sim12[,which(new_ind_js==k)])
                        })
                        x <- rep(0, num_domains)
                        x[js] <- res
                        x
                    }
                })
            }else if(method.domain=='BM.average'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_domains-1), .inorder=T, .combine=rbind), {
                    ind1 <- domains2terms_index[[i]]
                    progress_indicate(i, num_domains, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_domains
                        ind_js <- domains2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            x <- as.matrix(sim12[,which(new_ind_js==k)])
                            0.5*(mean(apply(x,1,max)) + mean(apply(x,2,max)))
                        })
                        x <- rep(0, num_domains)
                        x[js] <- res
                        x
                    }
                })
            }else if(method.domain=='BM.max'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_domains-1), .inorder=T, .combine=rbind), {
                    ind1 <- domains2terms_index[[i]]
                    progress_indicate(i, num_domains, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_domains
                        ind_js <- domains2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            x <- as.matrix(sim12[,which(new_ind_js==k)])
                            max(mean(apply(x,1,max)), mean(apply(x,2,max)))
                        })
                        x <- rep(0, num_domains)
                        x[js] <- res
                        x
                    }
                })
            }else if(method.domain=='BM.complete'){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_domains-1), .inorder=T, .combine=rbind), {
                    ind1 <- domains2terms_index[[i]]
                    progress_indicate(i, num_domains, 10, flag=T)
                    fast <- T
                    if(fast){
                        js <- (i+1):num_domains
                        ind_js <- domains2terms_index[js]
                        sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                        new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                        res <- sapply(1:length(ind_js), function(k){
                            x <- as.matrix(sim12[,which(new_ind_js==k)])
                            min(c(apply(x,1,max),apply(x,2,max)))
                        })
                        x <- rep(0, num_domains)
                        x[js] <- res
                        x
                    }
                })
            }

            ## add the last row
            sim <- rbind(sim, rep(0, num_domains))

            sim <- sim + Matrix::t(sim)
            sim <- Matrix::Matrix(sim, sparse=T)
        }
    }
    
    ###### non-parallel computing
    if(flag_parallel==F){
        ## calculate pair-wise semantic similarity between input domains
        sim <- Matrix::Matrix(0, nrow=length(domains), ncol=length(domains), sparse=T)
    
        ## print with possibly greater accuracy:
        ##op <- options(digits.secs = 6)
        ##options(op)
     
        if(method.domain=='average'){
            for(i in 1:(num_domains-1)){
                ind1 <- domains2terms_index[[i]]
                progress_indicate(i, num_domains, 10, flag=T)
                if(fast){
                    js <- (i+1):num_domains
                    ind_js <- domains2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        mean(sim12[,which(new_ind_js==k)])
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_domains){
                        ind2 <- domains2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- mean(sim12)
                    }
                }
            }
        }else if(method.domain=='max'){
            for(i in 1:(num_domains-1)){
                ind1 <- domains2terms_index[[i]]
                progress_indicate(i, num_domains, 10, flag=T)
                if(fast){
                    js <- (i+1):num_domains
                    ind_js <- domains2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        max(sim12[,which(new_ind_js==k)])
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_domains){
                        ind2 <- domains2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- max(sim12)
                    }
                }
            }
        }else if(method.domain=='BM.average'){
            for(i in 1:(num_domains-1)){
                ind1 <- domains2terms_index[[i]]
                progress_indicate(i, num_domains, 10, flag=T)
                if(fast){
                    js <- (i+1):num_domains
                    ind_js <- domains2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        x <- as.matrix(sim12[,which(new_ind_js==k)])
                        0.5*(mean(apply(x,1,max)) + mean(apply(x,2,max)))
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_domains){
                        ind2 <- domains2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- 0.5*(mean(apply(sim12,1,max)) + mean(apply(sim12,2,max)))
                    }
                }
            }
        
        }else if(method.domain=='BM.max'){
            for(i in 1:(num_domains-1)){
                ind1 <- domains2terms_index[[i]]
                progress_indicate(i, num_domains, 10, flag=T)
                if(fast){
                    js <- (i+1):num_domains
                    ind_js <- domains2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        x <- as.matrix(sim12[,which(new_ind_js==k)])
                        max(mean(apply(x,1,max)), mean(apply(x,2,max)))
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_domains){
                        ind2 <- domains2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- max(mean(apply(sim12,1,max)), mean(apply(sim12,2,max)))
                    }
                }
            }
        }else if(method.domain=='BM.complete'){
            for(i in 1:(num_domains-1)){
                ind1 <- domains2terms_index[[i]]
                progress_indicate(i, num_domains, 10, flag=T)
                if(fast){
                    js <- (i+1):num_domains
                    ind_js <- domains2terms_index[js]
                    sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
                    new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
                    res <- sapply(1:length(ind_js), function(k){
                        x <- as.matrix(sim12[,which(new_ind_js==k)])
                        min(c(apply(x,1,max),apply(x,2,max)))
                    })
                    sim[i,js] <- res
                }else{
                    for(j in (i+1):num_domains){
                        ind2 <- domains2terms_index[[j]]
                        ## pairwise similarity between terms
                        sim12 <- as.matrix(sim.term[ind1, ind2])
                        sim[i,j] <- min(c(apply(sim12,1,max),apply(sim12,2,max)))
                    }
                }
            }
        }
        sim <- sim + Matrix::t(sim)
    
    }
    
    rownames(sim) <- colnames(sim) <- domains
    
    ## create an object of class Dnetwork
    dnetwork <- as(sim, "Dnetwork")
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(dnetwork)
}
