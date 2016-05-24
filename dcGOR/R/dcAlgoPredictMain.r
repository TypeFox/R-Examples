#' Function to predict ontology terms given an input file containing domain architectures (including individual domains)
#'
#' \code{dcAlgoPredictMain} is supposed to predict ontology terms given an input file containing domain architectures (including individual domains).
#'
#' @param input.file an input file containing domain architectures (including individual domains). For example, a file containing UniProt ID and domain architectures for human proteins can be found in \url{http://dcgor.r-forge.r-project.org/data/Feature/hs.txt}. As seen in this example, the input file must contain the header (in the first row) and two columns: 1st column for 'SeqID' (actually these IDs can be anything), 2nd column for 'Architecture' (SCOP domain architectures, each represented as comma-separated domains). Alternatively, the input.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param output.file an output file containing predicted results. If not NULL, a tab-delimited text file will be also written out; otherwise, there is no output file (by default)
#' @param RData.HIS RData to load. This RData conveys two bits of information: 1) feature (domain) type; 2) ontology. It stores the hypergeometric scores (hscore) between features (individual domains or consecutive domain combinations) and ontology terms. The RData name tells which domain type and which ontology to use. It can be: SCOP sf domains/combinations (including "Feature2GOBP.sf", "Feature2GOMF.sf", "Feature2GOCC.sf", "Feature2HPPA.sf"), Pfam domains/combinations (including "Feature2GOBP.pfam", "Feature2GOMF.pfam", "Feature2GOCC.pfam", "Feature2HPPA.pfam"), InterPro domains (including "Feature2GOBP.interpro", "Feature2GOMF.interpro", "Feature2GOCC.interpro", "Feature2HPPA.interpro"). If NA, then the user has to input a customised RData-formatted file (see \code{RData.HIS.customised} below)
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
#' a data frame containing three columns: 1st column the same as the input file (e.g. 'SeqID'), 2nd for 'Term' (predicted ontology terms), 3rd for 'Score' (along with predicted scores)
#' @note
#' When 'output.file' is specified, a tab-delimited text file is written out, with the column names: 1st column the same as the input file (e.g. 'SeqID'), 2nd for 'Term' (predicted ontology terms), 3rd for 'Score' (along with predicted scores)
#' @export
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcAlgoPropagate}}, \code{\link{dcAlgoPredict}}
#' @include dcAlgoPredictMain.r
#' @examples
#' \dontrun{
#' # 1) Prepare an input file containing domain architectures
#' input.file <- "http://dcgor.r-forge.r-project.org/data/Feature/hs.txt"
#'
#' # 2) Do prediction using built-in data
#' output <- dcAlgoPredictMain(input.file, RData.HIS="Feature2GOMF.sf", parallel=FALSE)
#' output[1:5,]
#' 
#' # 3) Advanced usage: using customised data
#' x <- base::load(base::url("http://dcgor.r-forge.r-project.org/data/Feature2GOMF.sf.RData"))
#' RData.HIS.customised <- 'Feature2GOMF.sf.RData'
#' base::save(list=x, file=RData.HIS.customised)
#' #list.files(pattern='*.RData')
#' ## you will see an RData file 'Feature2GOMF.sf.RData' in local directory
#' output <- dcAlgoPredictMain(input.file, parallel=FALSE, RData.HIS.customised=RData.HIS.customised)
#' output[1:5,]
#' }

dcAlgoPredictMain <- function(input.file, output.file=NULL, RData.HIS=c(NA,"Feature2GOBP.sf","Feature2GOMF.sf","Feature2GOCC.sf","Feature2HPPA.sf","Feature2GOBP.pfam","Feature2GOMF.pfam","Feature2GOCC.pfam","Feature2HPPA.pfam","Feature2GOBP.interpro","Feature2GOMF.interpro","Feature2GOCC.interpro","Feature2HPPA.interpro"), merge.method=c("sum","max","sequential"), scale.method=c("log","linear","none"), feature.mode=c("supra","individual","comb"), slim.level=NULL, max.num=NULL, parallel=TRUE, multicores=NULL, verbose=T, RData.HIS.customised=NULL, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/dcGOR")
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
    tmp <- base::table(input[,2])
    copynum <- as.numeric(tmp)
    names(copynum) <- names(tmp)

    if(!is.na(RData.HIS)){
        tmp.RData.HIS <- RData.HIS
    }else if(file.exists(RData.HIS.customised)){
        tmp.RData.HIS <- RData.HIS.customised
    }else{
        stop("There is no input for HIS object! Please input one of two parameters ('RData.HIS' or 'RData.HIS.customised').\n")
    }

    if(verbose){
        message(sprintf("Predictions for %d sequences (with %d distinct architectures) using '%s' RData, '%s' merge method, '%s' scale method and '%s' feature mode (%s) ...", length(unique(input[,1])), length(copynum), tmp.RData.HIS, merge.method, scale.method, feature.mode, as.character(Sys.time())), appendLF=T)
        
        message(paste(c("\n##############################"), collapse=""), appendLF=T)
        message(paste(c("'dcAlgoPredict' is being called..."), collapse=""), appendLF=T)
        message(paste(c("##############################\n"), collapse=""), appendLF=T)
    }
    pscore <- dcAlgoPredict(data=names(copynum), RData.HIS=RData.HIS, merge.method=merge.method, scale.method=scale.method, feature.mode=feature.mode, slim.level=slim.level, max.num=max.num, parallel=parallel, multicores=multicores, verbose=verbose, RData.HIS.customised=RData.HIS.customised, RData.location=RData.location)
    
    if(verbose){
        message(paste(c("##############################"), collapse=""), appendLF=T)
        message(paste(c("'dcAlgoPredict' has been completed!"), collapse=""), appendLF=T)
        message(paste(c("##############################\n"), collapse=""), appendLF=T)
        
        message(sprintf("Preparations for output (%s)...", as.character(Sys.time())), appendLF=T)
    }
    # prepare the output
    tmp_seq <- input[,1]
    tmp_da <- input[,2]
    ind <- match(tmp_da, names(pscore))
    input_ps <- pscore[ind]
    output_list <- lapply(1:length(input_ps), function(i){
        x <- input_ps[[i]]
        if(!is.null(x)){
            #tmp_df <- cbind(rep(tmp_seq[i],length(x)), rep(tmp_da[i],length(x)), names(x), as.numeric(x))
            tmp_df <- cbind(rep(tmp_seq[i],length(x)), names(x), as.numeric(x))
            return(tmp_df)
        }else{
            return(NULL)
        }
    })
    output <- base::do.call(base::rbind, output_list)
    if(!is.null(output)){
    
        #colnames(output) <- c(colnames(input), "Term", "Score")
        colnames(output) <- c(colnames(input)[1], "Term", "Score")
    
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
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    
    invisible(output)
}