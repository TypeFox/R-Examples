#' Function to predict ontology terms for genomes with domain architectures (including individual domains)
#'
#' \code{dcAlgoPredictGenome} is supposed to predict ontology terms for genomes with domain architectures (including individual domains). 
#'
#' @param input.file an input file containing genomes and their domain architectures (including individual domains). For example, a file containing Hominidae genomes and their domain architectures can be found in \url{http://dcgor.r-forge.r-project.org/data/Feature/Hominidae.txt}. As seen in this example, the input file must contain the header (in the first row) and two columns: 1st column for 'Genome' (a genome like a container), 2nd column for 'Architecture' (SCOP domain architectures, each represented as comma-separated domains). Alternatively, the input.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param RData.HIS RData to load. This RData conveys two bits of information: 1) feature (domain) type; 2) ontology. It stores the hypergeometric scores (hscore) between features (individual domains or consecutive domain combinations) and ontology terms. The RData name tells which domain type and which ontology to use. It can be: SCOP sf domains/combinations (including "Feature2GOBP.sf", "Feature2GOMF.sf", "Feature2GOCC.sf", "Feature2HPPA.sf"), Pfam domains/combinations (including "Feature2GOBP.pfam", "Feature2GOMF.pfam", "Feature2GOCC.pfam", "Feature2HPPA.pfam"), InterPro domains (including "Feature2GOBP.interpro", "Feature2GOMF.interpro", "Feature2GOCC.interpro", "Feature2HPPA.interpro"). If NA, then the user has to input a customised RData-formatted file (see \code{RData.HIS.customised} below)
#' @param weight.method the method used how to weight predictions. It can be one of "none" (no weighting; by default), "copynum" for weighting copynumber of architectures, and "ic" for weighting information content (ic) of the term, "both" for weighting both copynumber and ic
#' @param merge.method the method used to merge predictions for each component feature (individual domains and their combinations derived from domain architecture). It can be one of "sum" for summing up, "max" for the maximum, and "sequential" for the sequential merging. The sequential merging is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} ranked highest hscore 
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
#' a matrix of terms X genomes, containing the predicted scores (per genome) as a whole
#' @note
#' none
#' @export
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcAlgoPropagate}}, \code{\link{dcAlgoPredict}}
#' @include dcAlgoPredictGenome.r
#' @examples
#' \dontrun{
#' # 1) Prepare an input file containing domain architectures
#' input.file <- "http://dcgor.r-forge.r-project.org/data/Feature/Hominidae.txt"
#'
#' # 2) Do prediction using built-in data
#' output <- dcAlgoPredictGenome(input.file, RData.HIS="Feature2GOMF.sf", parallel=FALSE)
#' dim(output)
#' output[1:10,]
#' 
#' # 3) Advanced usage: using customised data
#' x <- base::load(base::url("http://dcgor.r-forge.r-project.org/data/Feature2GOMF.sf.RData"))
#' RData.HIS.customised <- 'Feature2GOMF.sf.RData'
#' base::save(list=x, file=RData.HIS.customised)
#' #list.files(pattern='*.RData')
#' ## you will see an RData file 'Feature2GOMF.sf.RData' in local directory
#' output <- dcAlgoPredictGenome(input.file, parallel=FALSE, RData.HIS.customised=RData.HIS.customised)
#' dim(output)
#' output[1:10,]
#' }

dcAlgoPredictGenome <- function(input.file, RData.HIS=c(NULL,"Feature2GOBP.sf","Feature2GOMF.sf","Feature2GOCC.sf","Feature2HPPA.sf","Feature2GOBP.pfam","Feature2GOMF.pfam","Feature2GOCC.pfam","Feature2HPPA.pfam","Feature2GOBP.interpro","Feature2GOMF.interpro","Feature2GOCC.interpro","Feature2HPPA.interpro"), weight.method=c("none","copynum","ic","both"), merge.method=c("sum","max","sequential"), scale.method=c("log","linear","none"), feature.mode=c("supra","individual","comb"), slim.level=NULL, max.num=NULL, parallel=TRUE, multicores=NULL, verbose=T, RData.HIS.customised=NULL, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/dcGOR")
{

    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    RData.HIS <- match.arg(RData.HIS)
    merge.method <- match.arg(merge.method)
    scale.method <- match.arg(scale.method)
    weight.method <- match.arg(weight.method)
    feature.mode <- match.arg(feature.mode)
    
    ## for the interpro, only 'individual' are supported
    if(length(grep("interpro", RData.HIS, perl=T))!=0){
        #feature.mode <- "individual"
    }
    
    if(is.matrix(input.file) | is.data.frame(input.file)){
        if(verbose){
            message(sprintf("Load the input file ..."), appendLF=T)
        }
        if(is.data.frame(input.file)){
            input <- cbind(as.character(input.file[,1]), as.character(input.file[,2]))    
        }else{
            input <- input.file
        }
    }else if(is.character(input.file) & input.file!=''){
        if(verbose){
            message(sprintf("Read the input file '%s' ...", input.file), appendLF=T)
        }
        #tab <- read.delim(input.file, header=F, sep="\t", nrows=50, skip=1)
        #input <- read.table(input.file, header=F, sep="\t", skip=1, colClasses=sapply(tab,class))
        input <- utils::read.delim(input.file, header=T, sep="\t", colClasses="character")
    }else{
        return(NULL)
    }
    
    if(nrow(input)==0){
        return(NULL)
    }
    
    # determine the distinct architectures
    data <- sort(unique(input[,2]))

    if(!is.na(RData.HIS)){
        tmp.RData.HIS <- RData.HIS
    }else if(file.exists(RData.HIS.customised)){
        tmp.RData.HIS <- RData.HIS.customised
    }else{
        stop("There is no input for HIS object! Please input one of two parameters ('RData.HIS' or 'RData.HIS.customised').\n")
    }

    if(verbose){
        message(sprintf("Predictions for %d sequences (%d distinct architectures) using '%s' RData, '%s' merge method, '%s' scale method and '%s' feature mode (%s) ...", length(unique(input[,1])), length(data), tmp.RData.HIS, merge.method, scale.method, feature.mode, as.character(Sys.time())), appendLF=T)
        message(paste(c("\n##############################"), collapse=""), appendLF=T)
        message(paste(c("'dcAlgoPredict' is being called..."), collapse=""), appendLF=T)
        message(paste(c("##############################\n"), collapse=""), appendLF=T)
    }
    pscore <- dcAlgoPredict(data=data, RData.HIS=RData.HIS, merge.method=merge.method, scale.method=scale.method, feature.mode=feature.mode, slim.level=slim.level, max.num=max.num, parallel=parallel, multicores=multicores, verbose=verbose, RData.HIS.customised=RData.HIS.customised, RData.location=RData.location)
    
    ########################################################################################################
    if(verbose){
        message(paste(c("##############################"), collapse=""), appendLF=T)
        message(paste(c("'dcAlgoPredict' has been completed!"), collapse=""), appendLF=T)
        message(paste(c("##############################\n"), collapse=""), appendLF=T)
        message(sprintf("A summary in terms of ontology terms using '%s' weight method (%s)...", weight.method, as.character(Sys.time())), appendLF=T)
    }
    
    ##########################
    ## A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d (%s)", i, B, as.character(Sys.time())), appendLF=T)
            }
        }
    }
    
    ## A function to do prediction
    doSummary <- function(x, pscore, weight.method){
        
        # get copynum
        tmp <- base::table(x)
        copynum <- as.numeric(tmp)
        names(copynum) <- names(tmp)
        input_copynum <- as.numeric(copynum)
        
        # only those architectures in this genome
        pscore_tmp <- pscore[names(copynum)]
        output_list <- lapply(1:length(pscore_tmp), function(i){
            x <- pscore_tmp[[i]]
            if(!is.null(x)){
                ## copynum, term, score
                tmp_df <- cbind(rep(input_copynum[i],length(x)), names(x), as.numeric(x))
                return(tmp_df)
            }
        })
        
        ## split into two lists of terms: for the score, and for the copynum
        pscore_mat <- base::do.call(base::rbind, output_list)
        tscore <- base::split(x=as.numeric(pscore_mat[,3]), f=pscore_mat[,2])
        tcopynum <- base::split(x=as.numeric(pscore_mat[,1]), f=pscore_mat[,2])
        
        # give a summary score in terms of ontology terms
        ## whether to weight according to copynumber of architectures
        if(weight.method=="copynum" | weight.method=="both"){
            sscore <- sapply(1:length(tscore), function(i) {
                base::sum(tscore[[i]]*tcopynum[[i]])
            })
        }else{
            sscore <- sapply(1:length(tscore), function(i) {
                base::sum(tscore[[i]])
            })
        }
        names(sscore) <- names(tscore)
        ## whether to weight according to ic of terms
        if(weight.method=="ic" | weight.method=="both"){
            ind <- match(names(sscore), names(ic))
            sscore <- sscore*ic[ind]
        }
        ## scale to the range from 0 to 1
        sscore <- (sscore-min(sscore))/(max(sscore)-min(sscore))
        
        #### force min(sscore) to be 0.0001
        sscore[sscore==0] <- 0.0001
        sscore <- signif(sscore, digits=4)
        ####
        
        return(sscore)
    }
    
    ##########################
        
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
    ic <- x$ic
    
    ### split into a list of genomes
    #### genome, architecture
    tmp_genome <- base::split(x=input[,2], f=input[,1])
    
    ###### parallel computing
    flag_parallel <- F
    if(parallel){
        flag_parallel <- dnet::dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            j <- 1
            gn_sscore <- foreach::`%dopar%` (foreach::foreach(j=1:length(tmp_genome), .inorder=T), {
                progress_indicate(i=j, B=length(tmp_genome), 10, flag=T)
                x <- tmp_genome[[j]]
                doSummary(x=x, pscore=pscore, weight.method=weight.method)
            })
            names(gn_sscore) <- names(tmp_genome)
        }
    }
    
    ###### non-parallel computing
    if(flag_parallel==F){
        gn_sscore <- lapply(1:length(tmp_genome),function(j) {
            progress_indicate(i=j, B=length(tmp_genome), 10, flag=T)
            x <- tmp_genome[[j]]
            doSummary(x=x, pscore=pscore, weight.method=weight.method)
        })
        names(gn_sscore) <- names(tmp_genome)
    }
    
    #####################################################################################    
    #####################################################################################
    ## get a matrix of terms (all in ic) X genomes
    tmp <- ic
    tmp[] <- 0
    res_list <- lapply(gn_sscore, function(x){
        res <- tmp
        res[names(x)] <- as.numeric(x)
        res
    })
    tg_matrix <- base::do.call(base::cbind, res_list)
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(tg_matrix)
}