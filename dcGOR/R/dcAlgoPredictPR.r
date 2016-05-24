#' Function to assess the prediction performance via Precision-Recall (PR) analysis
#'
#' \code{dcAlgoPredictPR} is supposed to assess the prediction performance via Precision-Recall (PR) analysis. It requires two input files: 1) a Glod Standard Positive (GSP) file containing known annotations between proteins/genes and ontology terms; 2) a prediction file containing predicted terms for proteins/genes. Note: the known annotations will be recursively propagated towards the root of the ontology.
#'
#' @param GSP.file a Glod Standard Positive (GSP) file containing known annotations between proteins/genes and ontology terms. For example, a file containing annotations between human genes and HP terms can be found in \url{http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt}. As seen in this example, the input file must contain the header (in the first row) and two columns: 1st column for 'SeqID' (actually these IDs can be anything), 2nd column for 'termID' (HP terms). Alternatively, the GSP.file can be a matrix or data frame, assuming that GSP file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param prediction.file a prediction file containing proteins/genes, their predicted terms along with predictive scores. As seen in an example below, this file is usually created via \code{\link{dcAlgoPredictMain}}, containing three columns: 1st column for 'SeqID' (actually these IDs can be anything), 2nd column for 'Term' (ontology terms), 3rd column for 'Score' (predictive score). Alternatively, the prediction.file can be a matrix or data frame, assuming that prediction file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param ontology the ontology identity. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPON" for Human Phenotype ONset and clinical course, "MP" for Mammalian Phenotype, "EC" for Enzyme Commission, "KW" for UniProtKB KeyWords, "UP" for UniProtKB UniPathway. For details on the eligibility for pairs of input domain and ontology, please refer to the online Documentations at \url{http://supfam.org/dcGOR/docs.html}. If NA, then the user has to input a customised RData-formatted file (see \code{RData.ontology.customised} below)
#' @param num.threshold an integer to specify how many PR points (as a function of the score threshold) will be calculated
#' @param bin how to bin the scores. It can be "uniform" for binning scores with equal interval (ie with uniform distribution), and 'quantile' for binning scores with eual frequency (ie with equal number)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.ontology.customised a file name for RData-formatted file containing an object of S4 class 'Onto' (i.g. ontology). By default, it is NULL. It is only needed when the user wants to perform customised analysis using their own ontology. See \code{\link{dcBuildOnto}} for how to creat this object
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{dcRDataLoader}} for details
#' @return 
#' a data frame containing two columns: 1st column 'Precision' for precision, 2nd 'Recall' for recall. The row has the names corresponding to the score threshold.
#' @note
#' Prediction coverage: the ratio between predicted targets in number and GSP targets in number
#' F-measure: the maximum of a harmonic mean between precision and recall along PR curve
#' @export
#' @importFrom dnet dDAGannotate
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcConverter}}, \code{\link{dcDuplicated}}, \code{\link{dcAlgoPredictMain}}
#' @include dcAlgoPredictPR.r
#' @examples
#' \dontrun{
#' # 1) Generate prediction file with HPPA predicitions for human genes
#' architecture.file <- "http://dcgor.r-forge.r-project.org/data/Algo/SCOP_architecture.txt"
#' prediction.file <- "SCOP_architecture.HPPA_predicted.txt"
#' res <- dcAlgoPredictMain(input.file=architecture.file, output.file=prediction.file, RData.HIS="Feature2HPPA.sf", parallel=FALSE)
#'
#' # 2) Calculate Precision and Recall
#' GSP.file <- "http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt"
#' res_PR <- dcAlgoPredictPR(GSP.file=GSP.file, prediction.file=prediction.file, ontology="HPPA")
#' res_PR
#' 
#' # 3) Plot PR-curve
#' plot(res_PR[,2], res_PR[,1], xlim=c(0,1), ylim=c(0,1), type="b", xlab="Recall", ylab="Precision")
#' }

dcAlgoPredictPR <- function(GSP.file, prediction.file, ontology=c(NA,"GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPON","MP","EC","KW","UP"), num.threshold=10, bin=c("uniform","quantile"), verbose=T, RData.ontology.customised=NULL, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/dcGOR")
{

    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    ontology <- match.arg(ontology)
    bin <- match.arg(bin)
    
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
    
    ###############################################################################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("Second, import files for GSP and predictions (%s) ...", as.character(now)), appendLF=T)
    }

    #####################    
    ## import annotations
    if(is.matrix(GSP.file) | is.data.frame(GSP.file)){
        if(is.data.frame(GSP.file)){
            gsp <- cbind(as.character(GSP.file[,1]), as.character(GSP.file[,2]))
        }else{
            gsp <- GSP.file
        }
    }else if(is.character(GSP.file) & GSP.file!='' & !is.null(GSP.file) & !is.na(GSP.file)){
        #tab <- read.delim(GSP.file, header=F, sep="\t", nrows=50, skip=1)
        #gsp <- read.table(GSP.file, header=F, sep="\t", skip=1, colClasses=sapply(tab,class))
        gsp <- utils::read.delim(GSP.file, header=T, sep="\t", colClasses="character")
    }else{
        stop("The file 'GSP.file' must be provided!\n")
    }
    
    #####################
    ## import predictions
    if(is.matrix(prediction.file) | is.data.frame(prediction.file)){
        if(is.data.frame(prediction.file)){
            pred <- cbind(as.character(prediction.file[,1]), as.character(prediction.file[,2]), as.character(prediction.file[,3]))
        }else{
            pred <- prediction.file
        }
    }else if(is.character(prediction.file) & prediction.file!='' & !is.null(prediction.file) & !is.na(prediction.file)){
        #tab <- read.delim(prediction.file, header=F, sep="\t", nrows=50, skip=1)
        #pred <- read.table(prediction.file, header=F, sep="\t", skip=1, colClasses=sapply(tab,class))
        pred <- utils::read.delim(prediction.file, header=T, sep="\t", colClasses="character")
    }else{
        stop("The file 'prediction.file' must be provided!\n")
    }
    pred <- pred[!is.na(pred[,3]),]
    
    ## replace proteins with internal id
    all <- unique(c(unique(gsp[,1]), unique(pred[,1])))
    gsp[,1] <- match(gsp[,1], all)
    pred[,1] <- match(pred[,1], all)
    
    ###############################################################################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("Third, propagate GSP annotations (%s) ...", as.character(now)), appendLF=T)
    }
    
    ## for annotations, do propagation till the root of the ontology
    if(0){
        ## remove duplicated (if any)
        ind <- dcDuplicated(gsp, pattern.wise="row", verbose=verbose)
        gsp <- gsp[unique(ind),]
    }
    ## convert into a list of terms, each containing annotated proteins/genes
    gsp.list <- base::split(x=as.numeric(gsp[,1]), f=gsp[,2])
    ## do propagation
    if(1){
        dag <- dnet::dDAGannotate(g=g, annotations=gsp.list, path.mode="all_paths", verbose=verbose)
        gsp.list <- lapply(V(dag)$annotations, function(x) as.numeric(x))
        names(gsp.list) <- V(dag)$name
    }
    
    ## convert into a list of proteins/genes, each containing terms
    x <- gsp.list
    x_names <- names(x)
    output_list <- lapply(1:length(x), function(i){
        tmp <- x[[i]]
        cbind(C1=rep(x_names[i],length(tmp)), C2=tmp)
    })
    x_mat <- base::do.call(base::rbind, output_list)
    ### split: terms, proteins/genes
    gsp.list.gene <- base::split(x=x_mat[,1], f=as.numeric(x_mat[,2]))
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\tThere are %d genes/proteins in GSP (%s).", length(gsp.list.gene), as.character(now)), appendLF=T)
    }
    
    ######################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("Fourth, process input predictions (%s) ...", as.character(now)), appendLF=T)
    }
    
    ## split into a list of proteins/genes
    ### proteins/genes, terms, score
    tmp_term <- base::split(x=pred[,2], f=as.numeric(pred[,1]))
    tmp_score <- base::split(x=as.numeric(pred[,3]), f=as.numeric(pred[,1]))
    pred.list.gene <- lapply(1:length(tmp_score), function(i) {
        x <- tmp_score[[i]]
        names(x) <- tmp_term[[i]]
        return(x)
    })
    names(pred.list.gene) <- names(tmp_score)
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\tThere are %d genes/proteins in predictions (%s).", length(pred.list.gene), as.character(now)), appendLF=T)
    }
    
    ######################################
    
    gsp.names <- names(gsp.list.gene)
    pred.names <- names(pred.list.gene)
    both.names <- base::intersect(pred.names, gsp.names)
    
    gsp.list.both <- gsp.list.gene[both.names]
    pred.list.both <- pred.list.gene[both.names]
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Fifth, calculate the precision and recall for each of %d predicted and GSP genes/proteins (%s).", length(both.names), as.character(now)), appendLF=T)
    }
    
    ## get all decision threshold
    #tmp <- as.numeric(pred[,3])
    tmp <- unlist(pred.list.both, use.names=F)
    tmp <- unlist(pred.list.gene, use.names=F)
    if(bin=='uniform'){
        max_pred <- base::max(tmp)
        min_pred <- base::min(tmp)
        t <- base::seq(from=max_pred, to=min_pred, length.out=num.threshold+1)
    }else if(bin=='quantile'){
        t <- as.vector( stats::quantile(x=tmp, probs=base::seq(from=1, to=0, length.out=num.threshold+1)) )
    }
    
    ## For each target protein/gene and decision threshold t, calculate the precision and recall
    res_list <- lapply(1:length(both.names), function(j){
        x_gsp <- gsp.list.both[[j]]
        x_pred <- pred.list.both[[j]]
        
        res <- sapply(t, function(x){
            ### a set of predicted terms with score greater than or equal to t
            ind <- which(x_pred>=x)
            ###########################
            called_names <- unique(names(ind)) # in case that one sequence has many different architectures
            ###########################
                        
            callP <- length(called_names)
            ### a set of predicted terms (with score greater than or equal to t) overlapped in GSP
            ind2 <- match(called_names, x_gsp)
            TP <- sum(!is.na(ind2))
            return(rbind(TP,callP))
        })
        
        x_pr <- res[1,] / res[2,]
        x_rc <- res[1,] / length(x_gsp)
        
        ##########
        #x_pr[is.na(x_pr)] <- 0 # in case there are no called
        ##########
        
        return(cbind(precision=x_pr, recall=x_rc))
    })
    names(res_list) <- both.names
    
    ######################################

    if(verbose){
        now <- Sys.time()
        message(sprintf("Finally, calculate the averaged precision and recall (%s).", as.character(now)), appendLF=T)
    }
    
    ## Precision at threshold t is averaged over proteins on which at least one prediction was made above threshold t
    pr_sum <- rep(0, length(t))
    pr_sum_deno <- rep(0, length(t))
    for(i in 1:length(res_list)){
        ###############
        # in case there are no called
        tmp <- res_list[[i]][,1]
        tmp[is.na(tmp)] <- 0
        ###############
        pr_sum <- pr_sum + tmp
        
        pr_sum_deno <- pr_sum_deno + !is.na(res_list[[i]][,1]) # only preditable seq are considered
    }
    pr_ave <- pr_sum / pr_sum_deno
    
    ## Recall at threshold t is averaged over all proteins in GSP
    rc_sum <- rep(0, length(t))
    for(i in 1:length(res_list)){
        rc_sum <- rc_sum + res_list[[i]][,2]
    }
    rc_ave <- rc_sum / length(gsp.names)
    
    ## Prediction coverage: the maximum ratio between predicted proteins in number and GSP proteins in number (over all thresholds t)
    Pcoverage <- max(pr_sum_deno) / length(gsp.names)
    
    ## F-measure: the maximum (over all thresholds t) of a harmonic mean between precision and recall
    Fmeasure <- base::max( (2 * pr_ave * rc_ave) / (pr_ave + rc_ave), na.rm=T)
    
    ##############################################################################################
    
    if(verbose){
        message(sprintf("In summary, Prediction coverage: %.2f (amongst %d targets in GSP), and F-measure: %.2f.", Pcoverage, length(gsp.names), Fmeasure), appendLF=T)
    }
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    res_PR <- data.frame(Precision=pr_ave, Recall=rc_ave, row.names=t)
    
    invisible(res_PR)
}