#' Select the COLOMBOS REST API version to be used for retrieving data
#'
#' @param version positive number 2 or 3 - 3 (current REST API version) as default
#'
#' @references http://colombos.net
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' switchVersion (version = 2) # switch from COLOMBOS REST API 3 to 2
#' }
#'
#' @export
switchVersion <- function(version = 3) {
    if (version==3) options("REST.version"="http://rest.colombos.net/")
    else if (version==2) options("REST.version"="http://rest.legacyv2.colombos.net/")
    else stop("Select COLOMBOS REST API version to use: 2 or 3 (3 default)")
    message(paste("COLOMBOS REST version", version, "REST URL", getOption("REST.version")))
}

#' Returns a character vector corresponding to the currently available organisms.
#'
#' @return A list containing the currently available organisms.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' 
#' library('Rcolombos')
#' listOrganisms()
#' 
#'
listOrganisms <- function () {
    r <- GET(getOption("REST.version"), path = "get_organisms")
    if (r$status_code != 200) {
        stop_for_status(r) # Check the request succeeded
    }
    else {
        tmp <- content(r) # Automatically parse the json output
        response <- data.frame(matrix(unlist(tmp$data, recursive=F), length(tmp$data), 3, byrow=T))
        for (i in 1:3) response[,i] <- sapply(response[,i], as.character)
        colnames(response) <- c("name", "description", "nickname")
        return(response)
    }
}

#' This method takes as parameter a single string, representing an organism,
#' and returns a character vector corresponding to the currently available organisms.
#'
#' @param organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @return A data.frame containing the locustag and description
#' of all the genes for the selected organism.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' listGenes()
#' }
#'
listGenes <- function(organism="ecoli") {
    r <- GET(getOption("REST.version"),path = paste("get_genes/",organism, sep=""))
    if (r$status_code != 200) {
        stop_for_status(r)    # Check the request succeeded
    }
    else {
        tmp <- content(r)
        response <- data.frame(matrix(unlist(tmp$data, recursive=F), length(tmp$data), 2, byrow=T))
        for (i in 1:2) {
            response[,i] = as.character(response[,i])
        }
        colnames(response) <- c("locustag", "gene_name")
        return(response)
    }
}

#' This method takes as parameter a single string, representing an organism,
#' and returns a character vector corresponding to the currently available organisms.
#'
#' @param organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @return A data.frame containing the contrasts and GSM
#' of all the contrasts for the selected organism.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' listContrasts()
#' }
#'
listContrasts <- function(organism="ecoli"){
    
    r <- GET(getOption("REST.version"),path = paste("get_contrasts/",organism, sep=""))
    if (r$status_code!= 200) {
        stop_for_status(r)    # Check the request succeeded
    }
    else {
        # Automatically parse the json output
        tmp <- content(r)
        response <- data.frame(matrix(unlist(tmp$data, recursive=F), length(tmp$data), 2, byrow=T))
        for (i in 1:2) {
            response[,i] = as.character(response[,i])
        }
        colnames(response) <- c("name", "description")
        return(response)
    }
}

#' This method allows to download/import the full compendium for the selected organism
#'
#' @param organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @param path A string indicating the path where the file will be either downloaded or read,
#' if already retrieved
#' 
#' @return A list containing two or three data.frames.
#'  
#' In case \code{\link{switchVersion}} is equal to 2:
#' \item{exprdata}{the full compendium for the selected organism}
#' \item{condannot}{The condition annotation for the selected organism}
#' 
#' In case \code{\link{switchVersion}} is equal to 3:
#' \item{exprdata}{the full compendium for the selected organism}
#' \item{refannot}{The condition annotation for the reference contrasts}
#' \item{testannot}{The condition annotation for the test contrasts}
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' hpylo <- getCompendium("hpylo")
#' }
#'
getCompendium <- function(organism="hpylo", path=NULL){
    if(is.null(path)) path <- getwd() else {}
    destfile <- paste(path,"/",organism, "_compendium_data.zip",sep="")
    if(!file.exists(destfile)){
        r <- GET(getOption("REST.version"),path = paste("get_organism_data/",organism, sep=""))
        if (r$status_code!= 200) {
            stop_for_status(r)
        } 
        else {
            tmp <- content(r)
            download.file( tmp$data, destfile )
            return(parseCompendium(destfile))
        }
    } else {}
    return(parseCompendium(destfile))
}
#' This method allows importing the full compendium for the selected organism from a local file
#'
#' @param destfile A character containing the full path of the downloaded file
#'
#' @return A list containing two or three data.frames.
#'  
#' In case \code{\link{switchVersion}} is equal to 2:
#' \item{exprdata}{the full compendium for the selected organism}
#' \item{condannot}{The condition annotation for the selected organism}
#' 
#' In case \code{\link{switchVersion}} is equal to 3:
#' \item{exprdata}{the full compendium for the selected organism}
#' \item{refannot}{The condition annotation for the reference contrasts}
#' \item{testannot}{The condition annotation for the test contrasts}
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' mtube <- parseCompendium("mtube_compendium_data.zip")
#' }
#'
parseCompendium <- function(destfile){
    out_dir <- strsplit(destfile, "\\.")[[1]][1]
    unzip(destfile, exdir=out_dir) # unzip the files in the proper directory 
    files <- dir(path=out_dir,pattern="colombos_[a-z]+_[a-z]+_[0-9]+.txt")
    temp <- paste(out_dir, files[grep("colombos_[a-z]+_exprdata_[0-9]+.txt", files)], sep="/")
    my_cols <- na.omit(scan(temp, nlines=1, sep="\t", what="c", na.strings="", quiet=TRUE))
    exprdata <- read.csv(temp, row.names=1, skip=7, stringsAsFactors=FALSE, sep="\t", header=FALSE)
    exprdata <- exprdata[,c(2:dim(exprdata)[[2]])] 
    colnames(exprdata) = my_cols; exprdata <- exprdata[,c(2:dim(exprdata)[[2]])]
    ## condition annotations 
    if(getOption("REST.version")=="http://rest.colombos.net/"){
        temp <- paste(out_dir, files[grep("colombos_[a-z]+_refannot_[0-9]+.txt", files)], sep="/")
        refannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
        temp <- paste(out_dir, files[grep("colombos_[a-z]+_testannot_[0-9]+.txt", files)], sep="/")
        testannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
        return( list(exprdata=exprdata, refannot=refannot, testannot=testannot) )

    } else if(getOption("REST.version")=="http://rest.legacyv2.colombos.net/") {
        temp <- paste(out_dir, files[grep("colombos_[a-z]+_condannot_[0-9]+.txt", files)], sep="/")
        condannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T, quote="")
        return( list(exprdata=exprdata, condannot=condannot) )
    } else return(NULL)
}



#' This method takes as parameter a string (the nickname of an organism) and returns a character vector 
#' corresponding to the currently available annotation type for the selected organism.
#'
#' @param organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @return A data.frame containing the name and description of the annotation
#' for the selected organism.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' listAnnotationTypes()
#' }
#'
listAnnotationTypes <- function(organism="ecoli"){
    
    r <- GET(getOption("REST.version"),path = paste("get_annotation_types/",organism, sep=""))
    if (r$status_code!= 200) {
        stop_for_status(r)
    }
    else {
        tmp <- content(r)
        response <- data.frame(matrix(unlist(tmp$data, recursive=F), length(tmp$data), 2, byrow=T))
        for (i in 1:2) {
            response[,i] = as.character(response[,i])
        }
        colnames(response) <- c("name", "description")
        return(response)
    }
}

#' This method takes a string containing the nickname for the selected organism and a string containing the annotation type
#' and return the available entities
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param annotation A character containing the selected annotation type: use \code{\link{listAnnotationTypes}} to display
#' the available types.
#' 
#' @return A vector containing the available entities for the selected annotation type.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  library("Rcolombos")
#'  pathway_entities <- listEntities(organism="bsubt", annotation="Pathway")
#'  Tr_entities <- listEntities("bsubt","Transcriptional regulation")
#' }
#'
listEntities <- function(organism="ecoli", annotation="Pathway"){
    if(is.null(annotation)) stop("Insert a string with the annotation type.\n See listAnnotationTypes for the available types.") else {}
    r <- GET(getOption("REST.version"),path = paste("get_entities",organism, gsub(" ","%20", annotation), sep="/"))
    if (r$status_code!= 200) {
        stop_for_status(r)
    } else {
        tmp <- content(r)
        return(unlist(tmp$data))
    }
}

#' This method allows to retrieve all the annotations for the Reference and Test conditions for a selected organism (nickname) and for a user specified contrast name. Please be aware that only one contrast is allowed in input. It returns a list containing both ReferenceAnnnotation and TestAnnotation.
#' and return the available entities
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param contrast_name annotation A character containing the selected contrast_name type: use \code{\link{listContrasts}} to display the available contrast names.
#' 
#' @return A list of two data.frame, ReferenceAnnnotation and TestAnnotation, containing 2 columns: both the properties and the values for the selected contrast
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  library("Rcolombos")
#'  out <- get_contrast_annotations(organism="bsubt", 
#'  contrast_name="GSM27217.ch2-vs-GSM27217.ch1")
#' }
#'
get_contrast_annotations <- function(organism="bsubt", contrast_name="GSM27217.ch2-vs-GSM27217.ch1"){
    if(is.null(contrast_name)) stop("Insert a string with contrast_name\n See listContrasts for the available contrast names.") else {}
    r <- GET(getOption("REST.version"),path = paste("get_contrast_annotations",organism, contrast_name, sep="/"))
    if (r$status_code!= 200) {
        stop_for_status(r)
    } else {
        tmp <- content(r)
        return(list(ReferenceAnnotation=setNames(do.call(rbind.data.frame, strsplit(unlist(tmp$data$ReferenceAnnotation),":")), c("property", "value")), TestAnnotation=setNames(do.call(rbind.data.frame, strsplit(unlist(tmp$data$TestAnnotation),":")), c("property", "value"))))
    }
}

