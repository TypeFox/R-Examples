#' Function to apply dcGO algorithm to infer domain-centric ontology
#'
#' \code{dcAlgo} is supposed to apply dcGO algorithm to infer domain-centric ontology from input files. It requires two input files: 1) an annotation file containing annotations between proteins/genes and ontology terms; 2) an architecture file containing domain architectures for proteins/genes.
#'
#' @param anno.file an annotation file containing annotations between proteins/genes and ontology terms. For example, a file containing annotations between human genes and HP terms can be found in \url{http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt}. As seen in this example, the input file must contain the header (in the first row) and two columns: 1st column for 'SeqID' (actually these IDs can be anything), 2nd column for 'termID' (HP terms). Note: the file should use the tab delimiter as the field separator between columns
#' @param architecture.file an architecture file containing domain architectures (including individual domains) for proteins/genes. For example, a file containing human genes and domain architectures can be found in \url{http://dcgor.r-forge.r-project.org/data/Algo/SCOP_architecture.txt}. As seen in this example, the input file must contain the header (in the first row) and two columns: 1st column for 'SeqID' (actually these IDs can be anything), 2nd column for 'Architecture' (SCOP domain architectures, each represented as comma-separated domains). Note: the file should use the tab delimiter as the field separator between columns
#' @param output.file an output file containing results. If not NULL, a tab-delimited text file will be also written out, with 1st column 'Feature_id' for features/domains, 2nd column 'Term_id' for ontology terms, 3rd column 'Score' for hypergeometric scores (indicative of strength for feature-term associations). Otherwise, there is no output file (by default)
#' @param ontology the ontology identity. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPON" for Human Phenotype ONset and clinical course, "MP" for Mammalian Phenotype, "EC" for Enzyme Commission, "KW" for UniProtKB KeyWords, "UP" for UniProtKB UniPathway. For details on the eligibility for pairs of input domain and ontology, please refer to the online Documentations at \url{http://supfam.org/dcGOR/docs.html}. If NA, then the user has to input a customised RData-formatted file (see \code{RData.ontology.customised} below)
#' @param feature.mode the mode of how to define the features thereof. It can be: "supra" for combinations of one or two successive domains (including individual domains; considering the order), "individual" for individual domains only, and "comb" for all possible combinations (including individual domains; ignoring the order)
#' @param min.overlap the minimum number of overlaps with each term in consideration. By default, it sets to a minimum of 3
#' @param fdr.cutoff the fdr cutoff to call the significant associations between features and terms. By default, it sets to 1e-3
#' @param hscore.type the type of defining the hypergeometric score. It can be: "zscore" for z-score (by default), "fdr" for fdr (after being transformed via \eqn{-1*log_2(fdr)})
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.ontology.customised a file name for RData-formatted file containing an object of S4 class 'Onto' (i.g. ontology). By default, it is NULL. It is only needed when the user wants to perform customised analysis using their own ontology. See \code{\link{dcBuildOnto}} for how to creat this object
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{dcRDataLoader}} for details
#' @return 
#' a data frame containing three columns: 1st column 'Feature_id' for features, 2nd 'Term_id' for terms, and 3rd 'Score' for the hypergeometric score indicative of strength of associations beteen features and terms
#' @note
#' When 'output.file' is specified, a tab-delimited text file is output, with the column names: 1st column 'Feature_id' for features, 2nd 'Term_id' for terms, and 3rd 'Score' for the hypergeometric score indicative of strength of associations beteen features and terms
#' @export
#' @importFrom dnet dDAGannotate dDAGinduce dCheckParallel
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcSplitArch}}, \code{\link{dcConverter}}, \code{\link{dcDuplicated}}, \code{\link{dcAlgoPropagate}}
#' @include dcAlgo.r
#' @examples
#' \dontrun{
#' # 1) Prepare input file: anno.file and architecture.file
#' anno.file <- "http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt"
#' architecture.file <- "http://dcgor.r-forge.r-project.org/data/Algo/SCOP_architecture.txt"
#'
#' # 2) Do inference using built-in ontology
#' res <- dcAlgo(anno.file, architecture.file, ontology="HPPA", feature.mode="supra", parallel=FALSE)
#' res[1:5,]
#' 
#' # 3) Advanced usage: using customised ontology
#' x <- base::load(base::url("http://dcgor.r-forge.r-project.org/data/onto.HPPA.RData"))
#' RData.ontology.customised <- 'onto.HPPA.RData'
#' base::save(list=x, file=RData.ontology.customised)
#' #list.files(pattern='*.RData')
#' ## you will see an RData file 'onto.HPPA.RData' in local directory
#' res <- dcAlgo(anno.file, architecture.file, feature.mode="supra", parallel=FALSE, RData.ontology.customised=RData.ontology.customised)
#' res[1:5,]
#' }

dcAlgo <- function(anno.file, architecture.file, output.file=NULL, ontology=c(NA,"GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPON","MP","EC","KW","UP"), feature.mode=c("supra","individual","comb"), min.overlap=3, fdr.cutoff=1e-3, hscore.type=c("zscore","fdr"), parallel=TRUE, multicores=NULL, verbose=T, RData.ontology.customised=NULL, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/dcGOR")
{

    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    ontology <- match.arg(ontology)
    feature.mode <- match.arg(feature.mode)
    hscore.type <- match.arg(hscore.type)
    
    p.adjust.method <- c("BH","BY","bonferroni","holm","hochberg","hommel")[1]
    min.overlap <- as.integer(min.overlap)
    min.overlap <- ifelse(min.overlap<0, 3, min.overlap)
    fdr.cutoff <- ifelse(fdr.cutoff<=0 | fdr.cutoff>1, 1e-3, fdr.cutoff)
    
    if(is.null(anno.file) | is.na(anno.file)){
        stop("The file 'anno.file' must be provided!\n")
    }
    
    if(is.null(architecture.file) | is.na(architecture.file)){
        stop("The file 'architecture.file' must be provided!\n")
    }
    
    if(!is.na(ontology)){
    
        if(verbose){
            now <- Sys.time()
            message(sprintf("First, load the ontology '%s' (%s) ...", ontology, as.character(now)), appendLF=T)
        }
        
        #########
        ## load ontology information
        g <- dcRDataLoader(paste('onto.', ontology, sep=''), RData.location=RData.location)
        if(class(g)=="Onto"){
            g <- dcConverter(g, from='Onto', to='igraph', verbose=F)
        }
        
    }else if(file.exists(RData.ontology.customised)){
    
        if(verbose){
            now <- Sys.time()
            message(sprintf("First, load customised ontology '%s' (%s)...", RData.ontology.customised, as.character(now)), appendLF=T)
        }
    
        ## load ontology informatio
        g <- ''
        eval(parse(text=paste("g <- get(load('", RData.ontology.customised,"'))", sep="")))
        if(class(g)=="Onto"){
            g <- dcConverter(g, from='Onto', to='igraph', verbose=F)
        }
        ontology <- RData.ontology.customised

    }else{
        stop("There is no input for ontology! Please input one of two parameters ('ontology' and 'RData.ontology.customised').\n")
    }
    
    ########################################################
    ## A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d (%s)", i, B, as.character(Sys.time())), appendLF=T)
            }
        }
    }
    
    ## A function to do prediction
    doEnrichment <- function(genes.group, name.group, anno.sparse, M, N, genes.term.parent.sparse, N_rel, p.adjust.method, min.overlap, fdr.cutoff){
        
        # num of success in sampling
        X <- as.vector(genes.group %*% anno.sparse)
        
        ##########
        # get x_ind to indicate only those having more than min.overlap will be calculated
        x_ind <- which(X>=min.overlap)
        ##########
        
        X <- X[x_ind]
        if(length(X)==0){
            return(NULL)
        }
                
        ###########################
        ## do enrichment analysis
        
        # num of sampling
        K <- sum(genes.group)
        
        M <- M[x_ind]
        
        ## for p-value
        p.value <- stats::phyper(X,M,N-M,K, lower.tail=F, log.p=F)
        ## for adjp
        adjpvals <- stats::p.adjust(p.value, method=p.adjust.method)
        ## Z-score based on theoretical calculation
        x.exp <- K*M/N
        var.exp <- K*M/N*(N-M)/N*(N-K)/(N-1)
        zscore <- suppressWarnings((X-x.exp)/sqrt(var.exp))
        
        ###########################
        ## do relative enrichment analysis
        
        # num of sampling but also in parental terms
        K_rel <- as.vector(genes.group %*% genes.term.parent.sparse)
        
        K_rel <- K_rel[x_ind]
        N_rel <- N_rel[x_ind]
        
        ## for relative p-value
        p.value.rel <- suppressWarnings(stats::phyper(X,M,N_rel-M,K_rel, lower.tail=F, log.p=F))
        p.value.rel[is.na(p.value.rel)] <- 1
        ## for relative adjp
        adjpvals.rel <- stats::p.adjust(p.value.rel, method=p.adjust.method)
        ## for relative Z-score based on theoretical calculation
        x.exp <- K_rel*M/N_rel
        var.exp <- K_rel*M/N_rel*(N_rel-M)/N_rel*(N_rel-K_rel)/(N_rel-1)
        zscore.rel <- suppressWarnings((X-x.exp)/sqrt(var.exp))
        
        ## force those with zscore.rel=NaN to have: zscore.rel=0; p.value.rel=1; adjpvals.rel=1
        flag <- is.na(zscore.rel)
        zscore.rel[flag] <- 0
        adjpvals.rel[flag] <- p.value.rel[flag] <- 1

        #########
        ## get final hscore and fdr
        hscore <- signif(apply(cbind(zscore,zscore.rel), 1, min), digits=2)
        #pval <- apply(cbind(p.value,p.value.rel), 1, max)
        #fdr <- stats::p.adjust(pval, method=p.adjust.method)
        fdr <- apply(cbind(adjpvals,adjpvals.rel), 1, max)

        flag <- fdr==0
        fdr[flag] <- min(fdr[!flag])
        
        res <- data.frame(Feature_id=rep(name.group,length(X)), Term_id=names(M), X, M, K=rep(K,length(X)), N=rep(N,length(X)), zscore, p.value, adjpvals, K_rel, N_rel, zscore.rel, p.value.rel, adjpvals.rel, hscore, fdr, row.names=NULL)
        
        if(hscore.type=='zscore'){

            ind <- which(res$X >= min.overlap & res$fdr < fdr.cutoff)
            if(length(ind)>=1){
                return(res[ind,c(1:2,15)])
            }else{
                return(NULL)
            }
            
        }else if(hscore.type=='fdr'){
            
            ind <- which(res$X >= min.overlap & res$fdr < fdr.cutoff)
            res$fdr <- signif(-1*log10(res$fdr), digits=2)
            if(length(ind)>=1){
                return(res[ind,c(1:2,16)])
            }else{
                return(NULL)
            }
            
        }

    }
    
    ###############################################################################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("Second, import files for annotations '%s' and architectures '%s' (%s) ...", anno.file, architecture.file, as.character(now)), appendLF=T)
    }
    
    ## import annotations
    #tab <- utils::read.delim(anno.file, header=F, sep="\t", nrows=50, skip=1)
    #anno <- utils::read.table(anno.file, header=F, sep="\t", skip=1, colClasses=sapply(tab,class))
    anno <- utils::read.delim(anno.file, header=T, sep="\t", colClasses="character")
    ## import architectures
    #tab <- utils::read.delim(architecture.file, header=F, sep="\t", nrows=50, skip=1)
    #arch <- utils::read.table(architecture.file, header=F, sep="\t", skip=1, colClasses=sapply(tab,class))
    arch <- utils::read.delim(architecture.file, header=T, sep="\t", colClasses="character")
    
    ## replace proteins with internal id
    all <- unique(c(unique(anno[,1]), unique(arch[,1])))
    anno[,1] <- match(anno[,1], all)
    arch[,1] <- match(arch[,1], all)
    
    ###############################################################################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("Third, propagate annotations (%s) ...", as.character(now)), appendLF=T)
    }
    
    ## for annotations, do propagation till the root of the ontology
    if(0){
        ## remove duplicated (if any)
        ind <- dcDuplicated(anno, pattern.wise="row", verbose=verbose)
        anno <- anno[unique(ind),]
    }
    ## convert into a list of terms, each containing annotated proteins/genes
    anno.list <- base::split(x=as.numeric(anno[,1]), f=anno[,2])
    ## do propagation
    if(1){
        dag <- dnet::dDAGannotate(g=g, annotations=anno.list, path.mode="all_paths", verbose=verbose)
        anno.list <- lapply(V(dag)$annotations, function(x) as.numeric(x))
        names(anno.list) <- V(dag)$name
    }
    ## filter based on "min.overlap"
    ind <- which(sapply(anno.list, length) >= min.overlap)
    #anno.list <- anno.list[ind]
    ## parent-child
    dag_sub <- dnet::dDAGinduce(g=dag, nodes_query=names(anno.list), path.mode="all_paths")
    pc <- igraph::get.edgelist(dag_sub, names=TRUE)
    cp_list <- base::split(x=as.character(pc[,1]), f=pc[,2])
    
    if(length(anno.list)==0){
        stop("\tThere is no term being used.\n")
    }else{
        if(verbose){
            now <- Sys.time()
            message(sprintf("\tThere are %d terms used (%s).", length(anno.list), as.character(now)), appendLF=T)
        }
    }
    
    ######################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("Fourth, define groups using feature mode '%s' (%s) ...", feature.mode, as.character(now)), appendLF=T)
    }
    
    ## for architectures, get groups
    if(0){
        ## remove duplicated (if any)
        ind <- dcDuplicated(as.numeric(arch[,1]), pattern.wise="row", verbose=verbose)
        arch <- arch[unique(ind),]
    }
    ## convert into a list of architectures, each containing proteins
    arch.list <- base::split(x=as.numeric(arch[,1]), f=arch[,2])
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\ti) split into features (%s) ...", as.character(now)), appendLF=T)
    }
    data <- as.character(sort(names(arch.list)))
    ## get supradomains: a list of arhcitectures, each containing all possible supradomains
    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
        flag_parallel <- dnet::dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            j <- 1
            supra <- foreach::`%dopar%` (foreach::foreach(j=1:length(data), .inorder=T), {
                progress_indicate(i=j, B=length(data), 10, flag=T)
                suppressMessages(dcSplitArch(da=data[j], feature.mode=feature.mode, sep=",", ignore="_gap_", verbose=verbose))
            })
            names(supra) <- data
        }
    }
    ###### non-parallel computing
    if(flag_parallel==F){
        supra <- lapply(1:length(data),function(j) {
            progress_indicate(i=j, B=length(data), 10, flag=T)
            suppressMessages(dcSplitArch(da=data[j], feature.mode=feature.mode, sep=",", ignore="_gap_", verbose=verbose))
        })
        names(supra) <- data
    }
    
    ################
    if(verbose){
        now <- Sys.time()
        message(sprintf("\tii) obtain feature-based groups (%s) ...", as.character(now)), appendLF=T)
    }
    # get a matrix with two columns: proteins, supradomains thereof
    if(1){
        ###### parallel computing
        flag_parallel <- F
        if(parallel==TRUE){
            flag_parallel <- dnet::dCheckParallel(multicores=multicores, verbose=verbose)
            if(flag_parallel){
                j <- 1
                output_list <- foreach::`%dopar%` (foreach::foreach(j=1:length(supra), .inorder=T), {
                    progress_indicate(i=j, B=length(supra), 10, flag=T)
                    x <- supra[[j]]
                    y <- arch.list[[j]]
                    rep_y <- rep(y, each=length(x))
                    rep_x <- rep(x, length(y))
                    cbind(rep_y, rep_x)
                })
            }
        }
        ###### non-parallel computing
        if(flag_parallel==F){
            output_list <- lapply(1:length(supra),function(j) {
                progress_indicate(i=j, B=length(supra), 10, flag=T)
                x <- supra[[j]]
                y <- arch.list[[j]]
                rep_y <- rep(y, each=length(x))
                rep_x <- rep(x, length(y))
                cbind(rep_y, rep_x)
            })
        }
    }
    protein2supra_mat <- base::do.call(base::rbind, output_list)
    
    ## convert into a list of supradomains (group), each containing proteins
    groups.list <- base::split(x=as.numeric(protein2supra_mat[,1]), f=protein2supra_mat[,2])
    
    ## filter based on "min.overlap"
    ind <- which(sapply(groups.list, length) >= min.overlap)
    groups.list <- groups.list[ind]
    
    if(length(groups.list)==0){
        stop("There is no feature being defined.\n")
    }else{
        if(verbose){
            now <- Sys.time()
            message(sprintf("\tThere are %d features used (%s).", length(groups.list), as.character(now)), appendLF=T)
        }
    }
    
    ##############################################################################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("Finally, estimate associations between %d features and %d terms, with %d min overlaps and %.1e fdr cutoff (%s) ...", length(groups.list), length(anno.list), min.overlap, fdr.cutoff, as.character(now)), appendLF=T)
    }
    
    genes.universe <- unique(unlist(groups.list))
    
    # For a term, get all genes annotated in its all direct parental terms
    terms <- names(anno.list)
    genes.term.parent <- lapply(terms, function(child){
        term.parent <- cp_list[[child]]
        if(!is.null(term.parent)){
            unique(unlist(anno.list[term.parent]))
        }else{
            NULL
        }
    })
    
    #########################################
    # get a sparse matrix of all versuse length(groups.list)
    sparse_list <- lapply(1:length(groups.list), function(j){
        x <- groups.list[[j]]
        tmp_df <- cbind(x, rep(j,length(x)), rep(1, length(x)))
        return(tmp_df)
    })
    sparse_mat <- base::do.call(base::rbind, sparse_list)
    groups.sparse <- Matrix::sparseMatrix(i=sparse_mat[,1], j=sparse_mat[,2], x=sparse_mat[,3], dims=c(length(all),length(groups.list)))
    colnames(groups.sparse) <- names(groups.list)

    # get a sparse matrix of all versuse length(anno.list)
    sparse_list <- lapply(1:length(anno.list), function(j){
        x <- anno.list[[j]]
        tmp_df <- cbind(x, rep(j,length(x)), rep(1, length(x)))
        return(tmp_df)
    })
    sparse_mat <- base::do.call(base::rbind, sparse_list)
    anno.sparse <- Matrix::sparseMatrix(i=sparse_mat[,1], j=sparse_mat[,2], x=sparse_mat[,3], dims=c(length(all),length(anno.list)))
    colnames(anno.sparse) <- names(anno.list)
    
    # get a sparse matrix of all versuse length(genes.term.parent)
    sparse_list <- lapply(1:length(genes.term.parent), function(j){
        x <- genes.term.parent[[j]]
        tmp_df <- cbind(x, rep(j,length(x)), rep(1, length(x)))
        return(tmp_df)
    })
    sparse_mat <- base::do.call(base::rbind, sparse_list)
    genes.term.parent.sparse <- Matrix::sparseMatrix(i=sparse_mat[,1], j=sparse_mat[,2], x=sparse_mat[,3], dims=c(length(all),length(genes.term.parent)))
    #########################################
    
    # num of success in background
    M <- sapply(anno.list, length)
    # num in background
    N <- length(genes.universe)

    # num in parental terms as backgroun
    N_rel <- sapply(genes.term.parent, length)
    
    #########################
    ## do enrichment analysis
    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
        flag_parallel <- dnet::dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            j <- 1
            B <- dim(groups.sparse)[2]
            groups.name <- colnames(groups.sparse)
            res_list <- foreach::`%dopar%` (foreach::foreach(j=1:B, .inorder=T), {
                progress_indicate(i=j, B=B, 10, flag=T)
                genes.group <- groups.sparse[,j]
                name.group <- groups.name[j]
                doEnrichment(genes.group, name.group, anno.sparse, M, N, genes.term.parent.sparse, N_rel, p.adjust.method, min.overlap, fdr.cutoff)
            })
            names(res_list) <- colnames(groups.sparse)[1:B]
        }
    }
    
    ###### non-parallel computing
    if(flag_parallel==F){
        B <- dim(groups.sparse)[2]
        groups.name <- colnames(groups.sparse)
        res_list <- lapply(1:B,function(j) {
            progress_indicate(i=j, B=B, 10, flag=T)
            genes.group <- groups.sparse[,j]
            name.group <- groups.name[j]
            doEnrichment(genes.group, name.group, anno.sparse, M, N, genes.term.parent.sparse, N_rel, p.adjust.method, min.overlap, fdr.cutoff)
        })
        names(res_list) <- colnames(groups.sparse)[1:B]
    }
    
    feature2term_score <- base::do.call(base::rbind, res_list)
    rownames(feature2term_score) <- NULL
    colnames(feature2term_score) <- c("Feature_id", "Term_id", "Score")
    
    # how to deal with non-positive scores
    if(1){
        ## only retain all entries with positive score (otherwise, there will be problematic for prediction)
        feature2term_score <- feature2term_score[as.numeric(feature2term_score[,3])>0, ]
    }
    
    if(!is.null(output.file)){
    
        output <- feature2term_score
        utils::write.table(output, file=output.file, quote=F, row.names=F, sep="\t")
        
        if(file.exists(output.file)){
            message(sprintf("The results have been saved into '%s'.", file.path(getwd(),output.file)), appendLF=T)
        }
    }
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(feature2term_score)
}