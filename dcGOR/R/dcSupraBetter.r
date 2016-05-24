#' Function to find supra-domains with better scores than their individual domains
#'
#' \code{dcSupraBetter} is supposed to find supra-domains with better scores than their individual domains.
#'
#' @param input.file an input file used to build the object. This input file contains original annotations between domains/features and ontology terms, along with the hypergeometric scores (hscore) in support for their annotations. For example, a file containing original annotations between SCOP domain architectures and GO terms can be found in \url{http://dcgor.r-forge.r-project.org/data/Feature/Feature2GO.sf.txt}. As seen in this example, the input file must contain the header (in the first row) and three columns: 1st column for 'Feature_id' (here SCOP domain architectures), 2nd column for 'Term_id' (GO terms), and 3rd column for 'Score' (hscore)
#' @param output.file an output file containing results. If not NULL, a tab-delimited text file will be also written out, with 1st column 'Feature_id' for features/domains, 2nd column 'Term_id' for ontology terms, 3rd column 'Score' for hypergeometric scores (indicative of strength for feature-term associations). Otherwise, there is no output file (by default)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return 
#' a data frame containing three columns: 1st column 'Feature_id' for features, 2nd 'Term_id' for terms, and 3rd 'Score' for the hypergeometric score indicative of strength of associations beteen features and terms
#' @note
#' When 'output.file' is specified, a tab-delimited text file is output, with the column names: 1st column 'Feature_id' for features, 2nd 'Term_id' for terms, and 3rd 'Score' for the hypergeometric score indicative of strength of associations beteen features and terms
#' @export
#' @seealso \code{\link{dcList2Matrix}}
#' @include dcSupraBetter.r
#' @examples
#' \dontrun{
#' input.file <- "http://dcgor.r-forge.r-project.org/data/Feature/Feature2GO.sf.txt"
#' res <- dcSupraBetter(input.file)
#' res[1:10,]
#' }

dcSupraBetter <- function(input.file, output.file=NULL, verbose=T)
{
    
    if(is.null(input.file) | is.na(input.file)){
        stop("The file 'input.file' must be provided!\n")
    }
    
    if(verbose){
        message(sprintf("Reading the file '%s' ...", input.file), appendLF=T)
    }
    
    #tab <- read.delim(input.file, header=F, sep="\t", nrows=50, skip=1)
    #input <- read.table(input.file, header=F, sep="\t", skip=1, colClasses=sapply(tab,class))
    input <- utils::read.delim(input.file, header=T, sep="\t", colClasses="character")
    
    ## Feature_id, Term_id, Score
    tmp_term <- base::split(x=input[,2], f=input[,1], drop=T)
    tmp_score <- base::split(x=as.numeric(input[,3]), f=input[,1], drop=T)
    flist <- lapply(1:length(tmp_score), function(i){
        x <- tmp_score[[i]]
        names(x) <- tmp_term[[i]]
        return(x)
    })
    names(flist) <- names(tmp_score)
    
    fnames <- names(flist)
    tmp <- base::strsplit(fnames, ',')
    tmp_length <- sapply(tmp, length)
    flist_single <- flist[tmp_length==1]
    flist_multiple <- flist[tmp_length>=2]
    
    if(verbose){
        message(sprintf("There are %d supra-domains and %d individual domains.", length(flist_multiple), length(flist_single)), appendLF=T)
    }
    
    ## find all supra with better score
    names_single <- names(flist_single)
    names_multilpe <- names(flist_multiple)
    res_list <- lapply(1:length(flist_multiple), function(i){
        x <- flist_multiple[[i]]
        x_name <- names_multilpe[i]
        part_names <- unique(unlist(base::strsplit(x_name, ',')))
    
        flag <- match(part_names, names_single)
        tmp <- flag[!is.na(flag)]
        if(length(tmp)>=1){
            part_names <- names_single[tmp]
        }else{
            return(NULL)
        }
    
        part_list <- flist_single[part_names]
        res <- suppressMessages(dcList2Matrix(part_list))
        ### Feature_id, Term_id, Score
        rres <- base::split(x=as.numeric(res[,3]), f=res[,2])
    
        ind_list <- lapply(1:length(x), function(k){
            y <- x[k]
            tmp <- rres[[names(y)]]
            if(!is.null(tmp)){
                #### all has a better score than individuals
                if(base::all(y>tmp)){
                    return(k)
                }
            }else{
                return(k)
            }
        })
    
        ind <- unlist(ind_list)
        if(is.null(ind)){
            return(NULL)
        }else{
            return(x[ind])
        }
    
    })
    names(res_list) <- names_multilpe
    
    supra_mat_additional <- suppressMessages(dcList2Matrix(res_list))
    colnames(supra_mat_additional) <- c("Feature_id","Term_id","Score")
    
    if(verbose){
        message(sprintf("A total of %d annotations (for %d supra-domains) have better scores.", nrow(supra_mat_additional), length(unique(supra_mat_additional[,1]))), appendLF=T)
    }
    
    if(!is.null(output.file)){
    
        output <- supra_mat_additional
        utils::write.table(output, file=output.file, quote=F, row.names=F, sep="\t")
        
        if(file.exists(output.file)){
            message(sprintf("The results have been saved into '%s'.", file.path(getwd(),output.file)), appendLF=T)
        }
    }

    invisible(supra_mat_additional)
}