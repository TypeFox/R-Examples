#' Function to predict ontology terms given domain architectures (including individual domains)
#'
#' \code{dcAlgoPredict} is supposed to predict ontology terms given domain architectures (including individual domains). It involves 3 steps: 1) splitting an architecture into individual domains and all possible consecutive domain combinations (viewed as component features); 2) merging hscores among component features; 3) scaling merged hscores into predictive scores across terms.
#'
#' @param data an input data vector containing domain architectures. An architecture is represented in the form of comma-separated domains
#' @param RData.HIS RData to load. This RData conveys two bits of information: 1) feature (domain) type; 2) ontology. It stores the hypergeometric scores (hscore) between features (individual domains or consecutive domain combinations) and ontology terms. The RData name tells which domain type and which ontology to use. It can be: SCOP sf domains/combinations (including "Feature2GOBP.sf", "Feature2GOMF.sf", "Feature2GOCC.sf", "Feature2HPPA.sf"), Pfam domains/combinations (including "Feature2GOBP.pfam", "Feature2GOMF.pfam", "Feature2GOCC.pfam", "Feature2HPPA.pfam"), InterPro domains (including "Feature2GOBP.interpro", "Feature2GOMF.interpro", "Feature2GOCC.interpro", "Feature2HPPA.interpro"). If NA, then the user has to input a customised RData-formatted file (see \code{RData.HIS.customised} below)
#' @param merge.method the method used to merge predictions for each component feature (individual domains and their combinations derived from domain architecture). It can be one of "sum" for summing up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} ranked highest hscore 
#' @param scale.method the method used to scale the predictive scores. It can be: "none" for no scaling, "linear" for being linearily scaled into the range between 0 and 1, "log" for the same as "linear" but being first log-transformed before being scaled. The scaling between 0 and 1 is done via: \eqn{\frac{S - S_{min}}{S_{max} - S_{min}}}, where \eqn{S_{min}} and \eqn{S_{max}} are the minimum and maximum values for \eqn{S}
#' @param feature.mode the mode of how to define the features thereof. It can be: "supra" for combinations of one or two successive domains (including individual domains; considering the order), "individual" for individual domains only, and "comb" for all possible combinations (including individual domains; ignoring the order)
#' @param slim.level whether only slim terms are returned. By defaut, it is NULL and all predicted terms will be reported. If it is specified as a vector containing any values from 1 to 4, then only slim terms at these levels will be reported. Here is the meaning of these values: '1' for very general terms, '2' for general terms, '3' for specific terms, and '4' for very specific terms
#' @param max.num whether only top terms per sequence are returned. By defaut, it is NULL and no constraint is imposed. If an integer is specified, then all predicted terms (with scores in a decreasing order) beyond this number will be discarded. Notably, this parameter works after the preceding parameter \code{slim.level} 
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.HIS.customised a file name for RData-formatted file containing an object of S3 class 'HIS'. By default, it is NULL. It is only needed when the user wants to perform customised analysis. See \code{\link{dcAlgoPropagate}} on how this object is created
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{dcRDataLoader}} for details
#' @return 
#' a named list of architectures, each containing predictive scores
#' @note none
#' @export
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcSplitArch}}, \code{\link{dcConverter}}, \code{\link{dcAlgoPropagate}}, \code{\link{dcAlgoPredictMain}}, \code{\link{dcAlgoPredictGenome}}
#' @include dcAlgoPredict.r
#' @examples
#' \dontrun{
#' # 1) randomly generate 5 domains and/or domain architectures
#' x <- dcRDataLoader(RData="Feature2GOMF.sf")
#' data <- sample(names(x$hscore), 5)
#' 
#' # 2) get predictive scores of all predicted terms for this domain architecture
#' ## using 'sequential' method (by default)
#' pscore <- dcAlgoPredict(data=data, RData.HIS="Feature2GOMF.sf", parallel=FALSE)
#' ## using 'max' method
#' pscore_max <- dcAlgoPredict(data=data, RData.HIS="Feature2GOMF.sf", merge.method="max", parallel=FALSE)
#' ## using 'sum' method
#' pscore_sum <- dcAlgoPredict(data=data, RData.HIS="Feature2GOMF.sf", merge.method="sum", parallel=FALSE)
#' 
#' # 3) advanced usage
#' ## a) focus on those terms at the 2nd level (general)
#' pscore <- dcAlgoPredict(data=data, RData.HIS="Feature2GOMF.sf", slim.level=2, parallel=FALSE)
#' ## b) visualise predictive scores in the ontology hierarchy
#' ### load the ontology
#' g <- dcRDataLoader("onto.GOMF", verbose=FALSE)
#' ig <- dcConverter(g, from='Onto', to='igraph', verbose=FALSE)
#' ### do visualisation for the 1st architecture
#' data <- pscore[[1]]
#' subg <- dnet::dDAGinduce(ig, nodes_query=names(data), path.mode="shortest_paths")
#' dnet::visDAG(g=subg, data=data, node.info="term_id")
#' }

dcAlgoPredict <- function(data, RData.HIS=c(NA,"Feature2GOBP.sf","Feature2GOMF.sf","Feature2GOCC.sf","Feature2HPPA.sf","Feature2GOBP.pfam","Feature2GOMF.pfam","Feature2GOCC.pfam","Feature2HPPA.pfam","Feature2GOBP.interpro","Feature2GOMF.interpro","Feature2GOCC.interpro","Feature2HPPA.interpro"), merge.method=c("sum","max","sequential"), scale.method=c("log","linear","none"), feature.mode=c("supra","individual","comb"), slim.level=NULL, max.num=NULL, parallel=TRUE, multicores=NULL, verbose=T, RData.HIS.customised=NULL, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/dcGOR")
{

    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    RData.HIS <- match.arg(RData.HIS)
    merge.method <- match.arg(merge.method)
    scale.method <- match.arg(scale.method)
    feature.mode <- match.arg(feature.mode)
    
    if (is.vector(data)){
        data <- unique(data)
        data <- data[!is.null(data)]
        data <- data[!is.na(data)]
    }else{
        stop("The input data must be a vector.\n")
    }
    
    ## for the interpro, only 'domains' are supported
    if(length(grep("interpro", RData.HIS, perl=T))!=0){
        #feature.mode <- "domains"
    }
    
    ## load RData containing information on 'hscore', 'ic' and 'slim'
    if(!is.na(RData.HIS)){
        if(verbose){
            now <- Sys.time()
            message(sprintf("Load the HIS object '%s' (%s) ...", RData.HIS, as.character(now)), appendLF=T)
        }
        x <- dcRDataLoader(RData=RData.HIS, verbose=verbose, RData.location=RData.location)
    }else if(file.exists(RData.HIS.customised)){
        if(verbose){
            now <- Sys.time()
            message(sprintf("Load the customised HIS object '%s' (%s)...", RData.HIS.customised, as.character(now)), appendLF=T)
        }
        ## load ontology informatio
        x <- ''
        eval(parse(text=paste("x <- get(load('", RData.HIS.customised,"'))", sep="")))
        RData.HIS <- RData.HIS.customised
    }else{
        stop("There is no input for HIS object! Please input one of two parameters ('RData.HIS' or 'RData.HIS.customised').\n")
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
    doPredict <- function(da, hscore, merge.method, scale.method, feature.mode){
        
        res <- suppressMessages(dcSplitArch(da=da, feature.mode=feature.mode, sep=",", ignore="_gap_", verbose=verbose))
        #res <- dcSplitArch(da=da, feature.mode=feature.mode, sep=",", ignore="_gap_", verbose=verbose)
           
        ## no hscores for all features                 
        if(sum(res %in% names(hscore))==0){
            return(NULL)
        }
        
        ## convert feature lists to term lists               
        if(1){
            ind <- match(res, names(hscore))
            flist<- hscore[ind[!is.na(ind)]]
            #flist<- hscore[res]
            if(0){
                tmp <- cbind(names(flist[[1]]), as.numeric(flist[[1]]))
                if(length(flist)>=2){
                    for(i in 2:length(flist)){
                        tmp <- rbind(tmp, cbind(names(flist[[i]]), as.numeric(flist[[i]])))
                    }
                }
                ### split into a list of terms
                #### term, score
                tmp_score <- base::split(x=as.numeric(tmp[,2]), f=tmp[,1])
            }else{
                output_list <- lapply(1:length(flist), function(i){
                    x <- flist[[i]]
                    if(!is.null(x)){
                        tmp_df <- cbind(names(x), as.numeric(x))
                    }
                })
                tmp <- base::do.call(base::rbind, output_list)
                ### split into a list of terms
                #### term, score
                tmp_score <- base::split(x=as.numeric(tmp[,2]), f=tmp[,1])
            }

            
        }else{
            ind <- match(res, names(hscore))
            all <- base::unlist(hscore[ind[!is.na(ind)]])
            tmp_xxx <- sapply(base::strsplit(names(all),"[.]"), function(x){
                x
            })
            tmp <- t(rbind(tmp_xxx, as.numeric(all)))
            ### split into a list of terms
            #### feature, term, score
            tmp_score <- base::split(x=as.numeric(tmp[,3]), f=tmp[,2])
        }
        
        ## get predictive score 
        ### how to combine
        if(merge.method=='max'){
            score <- sapply(tmp_score, function(x) {
                base::max(x)
            })
        }else if(merge.method=='sum'){
            score <- sapply(tmp_score, function(x) {
                base::sum(x)
            })
        }else if(merge.method=='sequential'){
            score <- sapply(tmp_score, function(x) {
                base::sum(base::sort(x, decreasing=T) / (1:length(x)))
            })
        }
        
        ### how to scale 
        if(scale.method=='none'){
            score_scaled <- score
        }else if(scale.method=='linear'){
            score_scaled <- (score-min(score))/(max(score)-min(score))
            #### force min(score_scaled) to be 0.0001
            score_scaled[score_scaled==0] <- 0.0001
            score_scaled <- signif(score_scaled, digits=4)
            ####
        }else if(scale.method=='log'){
            score <- score[score>0] ## make sure that all scores are positive
            score <- log(score)
            score_scaled <- (score-min(score))/(max(score)-min(score))
            #### force min(score_scaled) to be 0.0001
            score_scaled[score_scaled==0] <- 0.0001
            score_scaled <- signif(score_scaled, digits=4)
            ####
        }
        
        return(base::sort(score_scaled, decreasing=T))
    }
    ########################################################
    
    if(verbose){
        message(sprintf("Predictions for %d architectures using '%s' merge method, '%s' scale method and '%s' feature mode (%s)...", length(data), merge.method, scale.method, feature.mode, as.character(Sys.time())), appendLF=T)
    }
    
    hscore <- x$hscore
    
    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
        flag_parallel <- dnet::dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            j <- 1
            pscore <- foreach::`%dopar%` (foreach::foreach(j=1:length(data), .inorder=T), {
                progress_indicate(i=j, B=length(data), 10, flag=T)
                da <- data[j]
                doPredict(da=da, hscore=hscore, merge.method=merge.method, scale.method=scale.method, feature.mode=feature.mode)
            })
            names(pscore) <- data
        }
    }
    
    ###### non-parallel computing
    if(flag_parallel==F){
        pscore <- lapply(1:length(data),function(j) {
            progress_indicate(i=j, B=length(data), 10, flag=T)
            da <- data[j]
            doPredict(da=da, hscore=hscore, merge.method=merge.method, scale.method=scale.method, feature.mode=feature.mode)
        })
        names(pscore) <- data
    }
    
    if(!is.null(slim.level)){
        slim.level <- as.integer(slim.level)
            
        terms <- base::unlist(sapply(x$slim[slim.level], function(x){names(x)}), use.names=F)
        if(!is.null(terms)){
            pscore <- lapply(pscore, function(x){
                ind <- match(names(x), terms)
                base::sort(x[!is.na(ind)], decreasing=T)
            })
            
            if(verbose){
                message(sprintf("Focus on predicted terms at '%s' slim level(s)", paste(as.character(slim.level),collapse=' ')), appendLF=T)
            }
            
        }
    }
    
    if(!is.null(max.num)){
        max.num <- as.integer(max.num)
        
        if(max.num >0){
            pscore <- lapply(pscore, function(x){
                if(length(x) > max.num){
                    val_cf <- x[max.num]
                    x[x>=val_cf]
                }else{
                    x
                }
            })
            
            if(verbose){
                message(sprintf("Focus on top %d predicted terms for each architecture", max.num), appendLF=T)
            }
        }

    }
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(pscore)
}