#' Function to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data
#'
#' \code{dDAGtermSim} is supposed to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data. Parallel computing is also supported for Linux or Mac operating systems.
#'
#' @param g an object of class "igraph" or "graphNEL". It must contain a vertex attribute called 'annotations' for storing annotation data (see example for howto)
#' @param terms the terms/nodes between which pair-wise semantic similarity is calculated. If NULL, all terms in the input DAG will be used for calcluation, which is very prohibitively expensive!
#' @param method the method used to measure semantic similarity between input terms. It can be "Resnik" for information content (IC) of most informative common ancestor (MICA) (see \url{http://arxiv.org/pdf/cmp-lg/9511007.pdf}), "Lin" for 2*IC at MICA divided by the sum of IC at pairs of terms (see \url{https://www.cse.iitb.ac.in/~cs626-449/Papers/WordSimilarity/3.pdf}), "Schlicker" for weighted version of 'Lin' by the 1-prob(MICA) (see \url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}), "Jiang" for 1 - difference between the sum of IC at pairs of terms and 2*IC at MICA (see \url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}), "Pesquita" for graph information content similarity related to Tanimoto-Jacard index (ie. summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 (see \url{http://www.ncbi.nlm.nih.gov/pubmed/18460186})). By default, it uses "Schlicker" method
#' @param fast logical to indicate whether a vectorised fast computation is used. By default, it sets to true. It is always advisable to use this vectorised fast computation; since the conventional computation is just used for understanding scripts
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a sparse matrix containing pair-wise semantic similarity between input terms. This sparse matrix can be converted to the full matrix via the function \code{as.matrix}
#' @note none
#' @export
#' @import Matrix
#' @seealso \code{\link{dDAGinduce}}, \code{\link{dDAGancestor}}, \code{\link{dDAGgeneSim}}, \code{\link{dCheckParallel}}
#' @include dDAGtermSim.r
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
#' # 4) calculate pair-wise semantic similarity between 5 randomly chosen terms 
#' terms <- sample(V(dag)$name, 5)
#' sim <- dDAGtermSim(g=dag, terms=terms, method="Schlicker", parallel=FALSE)
#' sim
#' }

dDAGtermSim <- function (g, terms=NULL, method=c("Resnik","Lin","Schlicker","Jiang","Pesquita"), fast=T, parallel=TRUE, multicores=NULL, verbose=T)
{
    
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    method <- match.arg(method)
    
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
    IC <- V(ig)$IC
    names(IC) <- V(ig)$name
    
    ## checking input terms
    terms <- terms[!is.na(terms)]
    if(is.null(terms) || is.na(terms)){
        terms <- V(ig)$name
    }else{
        flag <- terms %in% V(ig)$name
        if(sum(flag)!=0){
            terms <- terms[flag]
        }else{
            terms <- V(ig)$name
        }
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
    num_terms <- length(terms)
    if(verbose){
        message(sprintf("Calculate semantic similarity between %d terms using %s method (%s)...", num_terms, method, as.character(Sys.time())), appendLF=T)
    }
    
    ## pre-compute a sparse matrix of children x ancestors
    sCP <- dDAGancestor(ig, term1=terms, term2=NULL, verbose=T)
    
    allterms <- V(ig)$name
    ind <- match(terms, allterms)
    
    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
        flag_parallel <- dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
                
            if(method=="Resnik"){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:num_terms, .inorder=T, .combine=rbind), {
                    ancestor_i <- which(sCP[i,]==1)
                    progress_indicate(i, num_terms, 10, flag=T)
                    fast <- T
                    if(fast){
                        mat <- sCP[i:num_terms,]
                        ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                        flag <- is.element(ancestor_js[,2], ancestor_i)
                        ca_js <- ancestor_js[flag,]
                        ca_js_list <- split(ca_js[,2],ca_js[,1])
                        mica_js <- sapply(ca_js_list, function(x){
                            x[which.max(IC[x])]
                        })
                        js <- as.numeric(names(ca_js_list))+i-1
                        res <- IC[mica_js]
                        x <- rep(0, num_terms)
                        x[js] <- res
                        x
                    }
                })
            }else if(method=="Lin"){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:num_terms, .inorder=T, .combine=rbind), {
                    ancestor_i <- which(sCP[i,]==1)
                    progress_indicate(i, num_terms, 10, flag=T)
                    fast <- T
                    if(fast){
                        mat <- sCP[i:num_terms,]
                        ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                        flag <- is.element(ancestor_js[,2], ancestor_i)
                        ca_js <- ancestor_js[flag,]
                        ca_js_list <- split(ca_js[,2],ca_js[,1])
                        mica_js <- sapply(ca_js_list, function(x){
                            x[which.max(IC[x])]
                        })
                        js <- as.numeric(names(ca_js_list))+i-1
                        res <- 2*IC[mica_js]/(IC[ind[i]]+IC[ind[js]])
                        x <- rep(0, num_terms)
                        x[js] <- res
                        x
                    }
                })
            }else if(method=="Schlicker"){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:num_terms, .inorder=T, .combine=rbind), {
                    ancestor_i <- which(sCP[i,]==1)
                    progress_indicate(i, num_terms, 10, flag=T)
                    fast <- T
                    if(fast){
                        mat <- sCP[i:num_terms,]
                        ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                        flag <- is.element(ancestor_js[,2], ancestor_i)
                        ca_js <- ancestor_js[flag,]
                        ca_js_list <- split(ca_js[,2],ca_js[,1])
                        mica_js <- sapply(ca_js_list, function(x){
                            x[which.max(IC[x])]
                        })
                        js <- as.numeric(names(ca_js_list))+i-1
                        res <- (2*IC[mica_js]/(IC[ind[i]]+IC[ind[js]])) * (1 - 10^(-IC[mica_js]))
                        x <- rep(0, num_terms)
                        x[js] <- res
                        x
                    }
                })
            }else if(method=="Jiang"){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:num_terms, .inorder=T, .combine=rbind), {
                    ancestor_i <- which(sCP[i,]==1)
                    progress_indicate(i, num_terms, 10, flag=T)
                    fast <- T
                    if(fast){
                        mat <- sCP[i:num_terms,]
                        ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                        flag <- is.element(ancestor_js[,2], ancestor_i)
                        ca_js <- ancestor_js[flag,]
                        ca_js_list <- split(ca_js[,2],ca_js[,1])
                        mica_js <- sapply(ca_js_list, function(x){
                            x[which.max(IC[x])]
                        })
                        js <- as.numeric(names(ca_js_list))+i-1
                        tmp <- IC[ind[i]]+IC[ind[js]]-2*IC[mica_js]
                        tmp[tmp>1] <- 1
                        res <- 1- tmp
                        x <- rep(0, num_terms)
                        x[js] <- res
                        x
                    }
                })
            }else if(method=="Pesquita"){
                i <- 1
                sim <- foreach::`%dopar%` (foreach::foreach(i=1:num_terms, .inorder=T, .combine=rbind), {
                    ancestor_i <- which(sCP[i,]==1)
                    progress_indicate(i, num_terms, 10, flag=T)
                    fast <- T
                    if(fast){
                        mat <- sCP[i:num_terms,]
                        ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                        flag <- is.element(ancestor_js[,2], ancestor_i)
                        ca_js <- ancestor_js[flag,]
                        ca_js_list <- split(ca_js[,2],ca_js[,1])
                        mica_js <- sapply(ca_js_list, function(x){
                            x[which.max(IC[x])]
                        })
                        js <- as.numeric(names(ca_js_list))+i-1
                        ## graph information content similarity related to Tanimoto-Jacard index
                        ## summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2
                        ## for all ancestors
                        allan_js_list <- split(ancestor_js[,2],ancestor_js[,1])
                        allan_union <- sapply(allan_js_list, function(x){
                            ux <- union(x, ancestor_i)
                            sum(IC[ux])
                        })
                        ## for all common ancestors
                        allca_union <- sapply(ca_js_list, function(x){
                            sum(IC[x])
                        })
                        ## their ratio
                        res <- allca_union / allan_union
                        x <- rep(0, num_terms)
                        x[js] <- res
                        x
                    }
                })
    
            }

            sim <- sim + Matrix::t(sim)
            sim <- Matrix::Matrix(sim, sparse=T)
        }
    }

    ###### non-parallel computing
    if(flag_parallel==F){
        ## calculate pair-wise semantic similarity between input terms
        sim <- Matrix::Matrix(0, nrow=num_terms, ncol=num_terms, sparse=T)
        
        if(method=="Resnik"){
            for(i in 1:num_terms){
                ancestor_i <- which(sCP[i,]==1)
                progress_indicate(i, num_terms, 10, flag=T)
                if(fast){
                    mat <- sCP[i:num_terms,]
                    ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                    flag <- is.element(ancestor_js[,2], ancestor_i)
                    ca_js <- ancestor_js[flag,]
                    ca_js_list <- split(ca_js[,2],ca_js[,1])
                    mica_js <- sapply(ca_js_list, function(x){
                        x[which.max(IC[x])]
                    })
                    js <- as.numeric(names(ca_js_list))+i-1
                    res <- IC[mica_js]
                    sim[i,js] <- res
                }else{
                    for(j in i:num_terms){
                        ancestor_j <- which(sCP[j,]==1)
                        ## find the common ancestor with the maximal information content: most informative information ancestor (MICA)
                        ancestors <- intersect(ancestor_i, ancestor_j)
                        mica <- ancestors[which.max(IC[ancestors])]
                        res <- IC[mica]
                        sim[i,j] <- res
                    }
                }
            }
        }else if(method=="Lin"){
            for(i in 1:num_terms){
                ancestor_i <- which(sCP[i,]==1)
                progress_indicate(i, num_terms, 10, flag=T)
                if(fast){
                    mat <- sCP[i:num_terms,]
                    ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                    flag <- is.element(ancestor_js[,2], ancestor_i)
                    ca_js <- ancestor_js[flag,]
                    ca_js_list <- split(ca_js[,2],ca_js[,1])
                    mica_js <- sapply(ca_js_list, function(x){
                        x[which.max(IC[x])]
                    })
                    js <- as.numeric(names(ca_js_list))+i-1
                    res <- 2*IC[mica_js]/(IC[ind[i]]+IC[ind[js]])
                    sim[i,js] <- res
                }else{
                    for(j in i:num_terms){
                        ancestor_j <- which(sCP[j,]==1)
                        ## find the common ancestor with the maximal information content: most informative information ancestor (MICA)
                        ancestors <- intersect(ancestor_i, ancestor_j)
                        mica <- ancestors[which.max(IC[ancestors])]
                        res <- 2*IC[mica]/(IC[ind[i]]+IC[ind[j]])
                        sim[i,j] <- res
                    }
                }
            }
        }else if(method=="Schlicker"){
            for(i in 1:num_terms){
                ancestor_i <- which(sCP[i,]==1)
                progress_indicate(i, num_terms, 10, flag=T)
                if(fast){
                    mat <- sCP[i:num_terms,]
                    ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                    flag <- is.element(ancestor_js[,2], ancestor_i)
                    ca_js <- ancestor_js[flag,]
                    ca_js_list <- split(ca_js[,2],ca_js[,1])
                    mica_js <- sapply(ca_js_list, function(x){
                        x[which.max(IC[x])]
                    })
                    js <- as.numeric(names(ca_js_list))+i-1
                    res <- (2*IC[mica_js]/(IC[ind[i]]+IC[ind[js]])) * (1 - 10^(-IC[mica_js]))
                    sim[i,js] <- res
                }else{
                    for(j in i:num_terms){
                        ancestor_j <- which(sCP[j,]==1)
                        ## find the common ancestor with the maximal information content: most informative information ancestor (MICA)
                        ancestors <- intersect(ancestor_i, ancestor_j)
                        mica <- ancestors[which.max(IC[ancestors])]
                        res <- (2*IC[mica]/(IC[ind[i]]+IC[ind[j]])) * (1 - 10^(-IC[mica]))
                        sim[i,j] <- res
                    }
                }
            }
        }else if(method=="Jiang"){
            for(i in 1:num_terms){
                ancestor_i <- which(sCP[i,]==1)
                progress_indicate(i, num_terms, 10, flag=T)
                if(fast){
                    mat <- sCP[i:num_terms,]
                    ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                    flag <- is.element(ancestor_js[,2], ancestor_i)
                    ca_js <- ancestor_js[flag,]
                    ca_js_list <- split(ca_js[,2],ca_js[,1])
                    mica_js <- sapply(ca_js_list, function(x){
                        x[which.max(IC[x])]
                    })
                    js <- as.numeric(names(ca_js_list))+i-1
                    tmp <- IC[ind[i]]+IC[ind[js]]-2*IC[mica_js]
                    tmp[tmp>1] <- 1
                    res <- 1- tmp
                    sim[i,js] <- res
                }else{
                    for(j in i:num_terms){
                        ancestor_j <- which(sCP[j,]==1)
                        ## find the common ancestor with the maximal information content: most informative information ancestor (MICA)
                        ancestors <- intersect(ancestor_i, ancestor_j)
                        mica <- ancestors[which.max(IC[ancestors])]
                        res <- 1- min(1, IC[ind[i]]+IC[ind[j]]-2*IC[mica])
                        sim[i,j] <- res
                    }
                }
            }
        }else if(method=="Pesquita"){
            for(i in 1:num_terms){
                ancestor_i <- which(sCP[i,]==1)
                progress_indicate(i, num_terms, 10, flag=T)
                if(fast){
                    mat <- sCP[i:num_terms,]
                    ancestor_js <- which(matrix(mat==1, nrow=num_terms-i+1), arr.ind=T)
                    flag <- is.element(ancestor_js[,2], ancestor_i)
                    ca_js <- ancestor_js[flag,]
                    ca_js_list <- split(ca_js[,2],ca_js[,1])
                    mica_js <- sapply(ca_js_list, function(x){
                        x[which.max(IC[x])]
                    })
                    js <- as.numeric(names(ca_js_list))+i-1
                    ## graph information content similarity related to Tanimoto-Jacard index
                    ## summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2
                    ## for all ancestors
                    allan_js_list <- split(ancestor_js[,2],ancestor_js[,1])
                    allan_union <- sapply(allan_js_list, function(x){
                        ux <- union(x, ancestor_i)
                        sum(IC[ux])
                    })
                    ## for all common ancestors
                    allca_union <- sapply(ca_js_list, function(x){
                        sum(IC[x])
                    })
                    ## their ratio
                    res <- allca_union / allan_union
                    sim[i,js] <- res
                }else{
                    for(j in i:num_terms){
                        ancestor_j <- which(sCP[j,]==1)
                        ## find the common ancestor with the maximal information content: most informative information ancestor (MICA)
                        ancestors <- intersect(ancestor_i, ancestor_j)
                        mica <- ancestors[which.max(IC[ancestors])]
                        ## graph information content similarity related to Tanimoto-Jacard index
                        ## summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2
                        allancestors <- union(ancestor_i, ancestor_j)
                        res <- sum(IC[ancestors]) / sum(IC[allancestors])
                        sim[i,j] <- res
                    }
                }
            }
        }
        sim <- sim + Matrix::t(sim)
    }
    rownames(sim) <- colnames(sim) <- terms
    
    sim[as.matrix(is.na(sim))] <- 0
    
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
