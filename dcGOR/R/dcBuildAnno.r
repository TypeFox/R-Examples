#' Function to build an object of the S4 class Anno from input files
#'
#' \code{dcBuildAnno} is supposed to build an object of of the S4 class \code{\link{Anno}}, given input files. These input files include 1) a file containing domain information, 2) a file containing term information, and 3) a file containing associations between domains and terms. 
#'
#' @param domain_info.file an input file containing domain information. For example, a file containing InterPro domains (InterPro) can be found in \url{http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt}. As seen in this example, the input file must contain the header (in the first row), and entries in the first column intend to be domain ID (and must be unique). Note: the file should use the tab delimiter as the field separator between columns
#' @param term_info.file an input file containing term information. For example, a file containing Gene Ontology (GO) terms can be found in \url{http://dcgor.r-forge.r-project.org/data/InterPro/GO.txt}. As seen in this example, the input file must contain the header (in the first row) and four columns: 1st column for term ID (must be unique), 2nd column for term name, 3rd column for term namespace, and 4th column for term distance. These four columns must be provided, but the content for the last column can be arbitrary (if it is hard to prepare). Note: the file should use the tab delimiter as the field separator between columns
#' @param association.file an input file containing associations between domains and terms. For example, a file containing associations between InterPro domains and GO Molecular Function (GOMF) terms can be found in \url{http://dcgor.r-forge.r-project.org/data/InterPro/Domain2GOMF.txt}. As seen in this example, the input file must contain the header (in the first row) and two columns: 1st column for domain ID (corresponding to the first column in 'domain_info.file'), 2nd column for term ID (corresponding to the first column in 'term_info.file'). If there are additional columns, these columns will be ignored. Note: the file should use the tab delimiter as the field separator between columns
#' @param output.file an output file used to save the built object as an RData-formatted file. If NULL, this file will be saved into "Anno.RData" in the current working local directory
#' @return 
#' Any use-specified variable that is given on the right side of the assigement sign '<-', which contains the built \code{Anno} object.
#' Also, an RData file specified in "output.file" is saved in the local directory.  
#' @note If there are no use-specified variable that is given on the right side of the assigement sign '<-', then no object will be loaded onto the working environment.
#' @export
#' @seealso \code{\link{Anno}}
#' @include dcBuildAnno.r
#' @examples
#' \dontrun{
#' # build an "Anno" object that contains SCOP domain superfamilies (sf) annotated by GOBP terms
#' InterPro2GOMF <- dcBuildAnno(domain_info.file="http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt", term_info.file="http://dcgor.r-forge.r-project.org/data/InterPro/GO.txt", association.file="http://dcgor.r-forge.r-project.org/data/InterPro/Domain2GOMF.txt", output.file="InterPro2GOMF.RData")
#' InterPro2GOMF
#' }

dcBuildAnno <- function(domain_info.file, term_info.file, association.file, output.file="Anno.RData")
{

    if(is.null(domain_info.file) | is.na(domain_info.file)){
        stop("The file 'domain_info.file' must be provided!\n")
    }
    
    if(is.null(term_info.file) | is.na(term_info.file)){
        stop("The file 'term_info.file' must be provided!\n")
    }
    
    if(is.null(association.file) | is.na(association.file)){
        stop("The file 'association.file' must be provided!\n")
    }
    
    if(is.null(output.file)){
        warnings("Since the output file is not provided, the function will use the default output file 'Anno.RData'!\n")   
        output.file <- "Anno.RData"
    }
    
    
    ###########################
    # function to remove all those no-ASCII
    f <- function(x){
        if(is.numeric(x)){
            x
        }else{
            iconv(x, "latin1", "ASCII", sub="")
        }
    }
    ###########################
    
    ###################
    ## import association
    tab <- utils::read.delim(association.file, header=F, sep="\t", nrows=50, skip=1)
    association <- utils::read.table(association.file, header=F, sep="\t", skip=1, colClasses=sapply(tab,class))
    
    ## sparse matrix
    x <- association[,1]
    y <- association[,2]
    z <- rep(1, nrow(association))
    domains_in_asso <- sort(unique(x))
    terms_in_asso <- sort(unique(y))
    data <- Matrix::sparseMatrix(i=match(x,domains_in_asso), j=match(y,terms_in_asso), x=z)
    rownames(data) <- domains_in_asso
    colnames(data) <- terms_in_asso
    
    ###################
    ## import term info
    term_info <- utils::read.delim(term_info.file,header=T)
    term_info <- as.data.frame(apply(term_info,1:2,f))
    
    colnames(term_info) <- c("ID","Name","Namespace","Distance")
    rownames(term_info) <- term_info[,1]
    term_info <- term_info[order(term_info[,4]),] # order according to the 4th column
    
    ## focus those terms in common
    termID <- intersect(terms_in_asso, rownames(term_info))
    if(sum(is.na(suppressWarnings(as.numeric(termID)))) >= 1){
        termID <- sort(termID)
    }else{
        termID <- sort(as.numeric(termID))
    }
    
    ## define term info
    flag <- match(termID, rownames(term_info))
    term_info <- term_info[flag,]
    term_info <- term_info[order(term_info[,4]),] # order according to the 4th column
    
    ###################
    ## import domain info
    domain_info <- utils::read.delim(domain_info.file,header=T)
    rownames(domain_info) <- domain_info[,1]
    domain_info <- domain_info[order(domain_info[,1]),] # order according to the first column

    ## focus those domains in common
    domainID <- intersect(domains_in_asso, rownames(domain_info))
    if(sum(is.na(suppressWarnings(as.numeric(domainID)))) >= 1){
        domainID <- sort(domainID)
    }else{
        domainID <- sort(as.numeric(domainID))
    }
    
    ## define domain info
    flag <- match(domainID, rownames(domain_info))
    domain_info <- domain_info[flag,]
    domain_info <- domain_info[order(domain_info[,1]),] # order according to the first column
    
    ###################
    f_row <- match(rownames(data), rownames(domain_info))
    f_col <- match(colnames(data), rownames(term_info))
    annoData <- data[f_row, f_col]
    x <- new("Anno", annoData=annoData, termData=as(term_info,"InfoDataFrame"), domainData=as(domain_info,"InfoDataFrame"))
    
    # remove the RData extension 
    output.var <- gsub(".RData$", "", output.file, ignore.case=T, perl=T)
    output.var <- gsub(".RDat$", "", output.var, ignore.case=T, perl=T)
    output.var <- gsub(".RDa$", "", output.var, ignore.case=T, perl=T)
    
    do.call(assign, list(output.var, x))
    save(list=output.var, file=output.file)

    if(file.exists(output.file)){
        message(sprintf("An object of S4 class 'Anno' has been built and saved into '%s'.", file.path(getwd(),output.file)), appendLF=T)
    }

    invisible(x)
}