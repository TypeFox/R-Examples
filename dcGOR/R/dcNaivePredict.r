#' Function to perform naive prediction from input known annotations
#'
#' \code{dcNaivePredict} is supposed to perform naive prediction from input known annotations. For each gene/protein, a term to be predicted are simply the frequency of that term appearing in the known annotations.
#'
#' @param data an input vector containing genes/proteins to be predicted
#' @param GSP.file a Glod Standard Positive (GSP) file containing known annotations between proteins/genes and ontology terms. For example, a file containing annotations between human genes and HP terms can be found in \url{http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt}. As seen in this example, the input file must contain the header (in the first row) and two columns: 1st column for 'SeqID' (actually these IDs can be anything), 2nd column for 'termID' (HP terms). Alternatively, the GSP.file can be a matrix or data frame, assuming that GSP file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param output.file an output file containing predicted results. If not NULL, a tab-delimited text file will be also written out; otherwise, there is no output file (by default)
#' @param ontology the ontology identity. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPON" for Human Phenotype ONset and clinical course, "MP" for Mammalian Phenotype, "EC" for Enzyme Commission, "KW" for UniProtKB KeyWords, "UP" for UniProtKB UniPathway. For details on the eligibility for pairs of input domain and ontology, please refer to the online Documentations at \url{http://supfam.org/dcGOR/docs.html}. If NA, then the user has to input a customised RData-formatted file (see \code{RData.ontology.customised} below)
#' @param max.num an integer to specify how many terms will be predicted for each gene/protein
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.ontology.customised a file name for RData-formatted file containing an object of S4 class 'Onto' (i.g. ontology). By default, it is NULL. It is only needed when the user wants to perform customised analysis using their own ontology. See \code{\link{dcBuildOnto}} for how to creat this object
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{dcRDataLoader}} for details
#' @return 
#' a data frame containing three columns: 1st column the same as the input file (e.g. 'SeqID'), 2nd for 'Term' (predicted ontology terms), 3rd for 'Score' (along with predicted scores)
#' @note
#' When 'output.file' is specified, a tab-delimited text file is written out, with the column names: 1st column the same as the input file (e.g. 'SeqID'), 2nd for 'Term' (predicted ontology terms), 3rd for 'Score' (along with predicted scores).

#' @export
#' @importFrom dnet dDAGannotate
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcAlgoPropagate}}
#' @include dcNaivePredict.r
#' @examples
#' \dontrun{
#' # 1) prepare genes to be predicted
#' input.file <- "http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt"
#' #input.file <- "http://dcgor.r-forge.r-project.org/data/Algo/SCOP_architecture.txt"
#' input <- utils::read.delim(input.file, header=TRUE, sep="\t", colClasses="character")
#' data <- unique(input[,1])
#'
#' # 2) do naive prediction
#' GSP.file <- "http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt"
#' res <- dcNaivePredict(data=data, GSP.file=GSP.file, ontology="HPPA")
#' res[1:10,]
#' 
#' # 3) calculate Precision and Recall
#' res_PR <- dcAlgoPredictPR(GSP.file=GSP.file, prediction.file=res, ontology="HPPA")
#' res_PR
#' 
#' # 4) plot PR-curve
#' plot(res_PR[,2], res_PR[,1], xlim=c(0,1), ylim=c(0,1), type="b", xlab="Recall", ylab="Precision")
#' }

dcNaivePredict <- function(data, GSP.file, output.file=NULL, ontology=c(NA,"GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPON","MP","EC","KW","UP"), max.num=1000, verbose=T, RData.ontology.customised=NULL, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/dcGOR")
{

    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    ontology <- match.arg(ontology)
    
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
        message(sprintf("Second, import files for GSP (%s) ...", as.character(now)), appendLF=T)
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
    
    gsp.counts <- sapply(gsp.list, length)
    gsp.fraction <- sort(gsp.counts/max(gsp.counts), decreasing=T)
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\tThere are %d terms in GSP (%s).", length(gsp.fraction), as.character(now)), appendLF=T)
    }
    
    ######################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("Fourth, do naive predictions for %d genes/proteins (%s) ...", length(data), as.character(now)), appendLF=T)
    }
    
    if(!is.null(max.num)){
        max.num <- as.integer(max.num)
        
        if(max.num >0){
            x <- gsp.fraction
            if(length(x) > max.num){
                val_cf <- x[max.num]
                gsp.fraction <- x[x>=val_cf]
            }
            
            if(verbose){
                message(sprintf("\tFocus on top %d predicted terms for each gene/protein", max.num), appendLF=T)
            }
        }

    }
    
    SeqID <- rep(data, length(gsp.fraction))
    Term <- rep(names(gsp.fraction), length(data))
    Score <- rep(as.numeric(signif(gsp.fraction,digits=4)), length(data))
    output <- cbind(SeqID, Term, Score)
    
    if(!is.null(output)){
   
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