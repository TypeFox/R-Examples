#' @title Interface to Colombos Compendia using the Exposed REST API
#' @name Rcolombos
#' @description Provides programmatic access to Colombos, a web based
#'interface for exploring and analyzing comprehensive organism-specific
#'cross-platform expression compendia of bacterial organisms.
#' 
#' @docType package
#' @aliases Rcolombos Rcolombos-package
NULL

#' This method mimics the quick_search functionality of Colombos.
#' It takes a string containg the nickname for the selected organism and a vector of string 
#' representing the genes of interest for the specified organism and returns a list containing 
#' the locustags (gene_names), contrasts and M-values for the current selection.
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param genes A vector of strings representing the genes of interest.
#' @param geneNames boolean if == FALSE (default) returns the locustag otherwise the gene_name for the selected genes.
#' 
#' 
#' @return A data.frame containing locustag (gene_names), 
#' contrasts and M-values for the current organism and genes.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  library("Rcolombos")
#'  my_module <- quick_search(organism="ecoli", 
#'                      genes=c("b0400","b2805","b0567"), 
#'                      geneNames=FALSE)
#'  heatmap(as.matrix(my_module), col=terrain.colors(15))
#' }
#'
quick_search <- function(organism="ecoli", genes, geneNames=FALSE){
    if(is.null(genes)) stop("Insert a character vector with the genes to be imputed.") else {}
    r <- GET(getOption("REST.version"), path = paste("quick_search",organism, paste(genes, collapse=","), sep="/"))
    if (r$status_code != 200) {
        stop_for_status(r)
    } else {
        tmp <- content(r)
        response = suppressWarnings( data.frame(matrix(as.numeric(tmp$data$Mvalues), 
        nrow=length(tmp$data$geneLocustags), ncol=length(tmp$data$contrasts), byrow=T)) )
        colnames(response) <- tmp$data$contrasts; rownames(response) <- tmp$data$geneLocustags
        if (geneNames) rownames(response) <- tmp$data$geneNames
        return(response)
    }
}
#' This method mimics the advanced_search functionality of Colombos.
#' It takes a series of parameters, representing the different settings available on Colombos advanced search
#' and returns a list containing the locustags (gene_names), contrasts and M-values for the current selection.
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param g_ids A vector of strings representing contrast_id, go terms, experiment id or condition id according the search type.
#' @param geneNames boolean if == FALSE (default) return the locustag otherwise the gene_name for the selected genes.
#' @param c_ids A vector of strings representing contrast_id, go terms, experiment id or condition id according the search type.
#' @param by A string eithes genes, contrasts, both allowing the selection by genes entities, contrast entities or both.
#' @param g_search_type A string either genes, go or annotation.
#' @param ann_type A string containing the selected gene_annotation_type: use \code{\link{listEntities}} to display the available entities.
#' @param c_search_type A string either contrast_names. experiment, go, condition use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @return A data.frame containing locustag (gene_names), 
#' contrasts and M-values for the current organism and genes.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  library("Rcolombos")
#'  
#'  # modules by gene entities
#'  g.gn <- advanced_search(organism="bsubt", 
#'                      g_ids=c("cgeB","yfnG"), 
#'                      by="genes", g_search_type="genes")
#'  g.go <- advanced_search(organism="bsubt", 
#'                      g_ids="response to antibiotic, transcription", 
#'                      by="genes", g_search_type="go")
#'  g.anno <- advanced_search(organism="bsubt", 
#'                      g_ids="biotin-carboxyl carrier protein assembly", 
#'                      by="genes", g_search_type="annotation", ann_type="Pathway")
#'  
#'  # modules by contrast entities
#'  c.cn <- advanced_search(organism="bsubt", 
#'                      c_ids=c("GSM27217.ch2-vs-GSM27217.ch1","GSM27218.ch1-vs-GSM27218.ch2"),
#'                      by="contrasts", c_search_type="contrast_names")
#'  c.go <- advanced_search(organism="bsubt", 
#'                      c_ids="response to antibiotic, transcription", 
#'                      by="contrasts", c_search_type="go")
#'  c.exp <- advanced_search(organism="bsubt", 
#'                      c_ids="GSE22296", by="contrasts", c_search_type="experiment")
#'  c.cond <- advanced_search(organism="bsubt", 
#'  c_ids=c("DAPTOMYCIN","H2O2","HPUra","IPTG","MMC","MNCL2","MOENOMYCIN","RAMOPLANIN"),
#'  by="contrasts", c_search_type="condition")
#'  
#'  # modules by both gene and contrast entities
#'  b.go.cn <- advanced_search(organism="bsubt", 
#'                      g_ids="response to antibiotic, transcription", geneNames=F, 
#'                      c_ids=c("GSM27217.ch2-vs-GSM27217.ch1","GSM27218.ch1-vs-GSM27218.ch2"),
#'                      g_search_type="go", c_search_type="contrast_names", by="both")
#'  b.gn.ge <- advanced_search(organism="bsubt", g_ids=c("BSU00020","BSU00100"), 
#'                      geneNames=F, c_ids="GSE22296", g_search_type="genes", 
#'                      c_search_type="experiment", by="both")
#'  b.go.ge <- advanced_search(organism="bsubt", g_ids="response to antibiotic, transcription",
#'                      geneNames=F, c_ids="GSE22296", g_search_type="go", 
#'                      c_search_type="experiment", by="both")
#'  b.gn.cn <- advanced_search(organism="bsubt", 
#'                      g_ids=c("dnaA","dnaN","yaaA","recF","yaaB","gyrB"), geneNames=FALSE, 
#'                      c_ids=c("GSM27217.ch2-vs-GSM27217.ch1","GSM27218.ch1-vs-GSM27218.ch2",
#'                      "GSM27219.ch2-vs-GSM27219.ch1","GSM27278.ch2-vs-GSM27278.ch1",
#'                      "GSM27279.ch1-vs-GSM27279.ch2"), 
#'                      g_search_type="genes", c_search_type="contrast_names", by="both")
#'  heatmap(as.matrix(b.gn.cn), col=terrain.colors(15))
#' }
#'
advanced_search <- function(organism=NULL, g_ids=NULL, geneNames=FALSE, c_ids, by="genes", g_search_type, ann_type, c_search_type){
    if(is.null(organism)) stop("Insert a character vector corresponding to the nickname of the selected organism.") else {}
    if(by=="genes") out <- advanced_search_by_genes(organism, g_ids, geneNames, g_search_type, ann_type)
    else if(by=="contrasts") out <- advanced_search_by_contrasts(organism, c_ids, geneNames, c_search_type)
    else if(by=="both") out <- advanced_search_by_both(organism, g_ids, geneNames, c_ids, g_search_type, ann_type, c_search_type)
    else stop("Wrong by: it should be either genes, contrasts or both!")
}
#' Accessory function allowing the advanced_search by gene_ids, go, annotation
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param ids A vector of strings representing gene_id, go terms or annotation entities according the search type.
#' @param geneNames boolean if == FALSE (default) return the locustag otherwise the gene_name for the selected genes.
#' @param g_search_type A string either genes, go or annotation.
#' @param ann_type A string containing the selected gene_annotation_type: use \code{\link{listEntities}} to display the available entities.
#' 
#' @return A data.frame containing locustag (gene_names), 
#' contrasts and M-values for the current organism and genes.
#'
#' @references http://colombos.net
#'
#' @export
#'
advanced_search_by_genes <- function(organism="bsubt", ids=NULL, geneNames=FALSE, g_search_type="genes", ann_type){
    if(is.null(ids)) stop("Insert the ids for the specific search_type.") else {}
    if(g_search_type=="genes"){
        url_string <- paste(getOption("REST.version"),"advanced_search_by_genes", organism, "genes", paste(ids, collapse=","), sep="/")
    }
    else if(g_search_type=="go"){
        url_string <- paste(getOption("REST.version"),"advanced_search_by_genes", organism, "go", gsub(" ","%20", ids), sep="/")
    }
    else if(g_search_type=="annotation"){
        url_string <- paste(getOption("REST.version"),"advanced_search_by_genes", organism, "annotation", gsub(" ","%20", ann_type), 
                            gsub(" ","%20", ids), sep="/")  
    } else {
        stop("Wrong search_type: it should be either genes, go or annotation!")
    }
    r <- GET(url_string)
    if (r$status_code != 200) {
        stop_for_status(r)
    } else {
        tmp <- content(r)
        response = suppressWarnings( data.frame(matrix(as.numeric(tmp$data$Mvalues), 
        nrow=length(tmp$data$geneLocustags), ncol=length(tmp$data$contrasts), byrow=T)) )
        colnames(response) <- tmp$data$contrasts; rownames(response) <- tmp$data$geneLocustags       
        if (geneNames) rownames(response) <- tmp$data$geneNames
        return(response)
    }
}
#' Accessory function allowing the advanced_search by contrast_ids, go, experiment, condition
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param ids A vector of strings representing contrast_id, go terms, experiment id or condition id according the search type.
#' @param geneNames boolean if == FALSE (default) return the locustag otherwise the gene_name for the selected genes.
#' @param c_search_type A string either contrast_names. experiment, go, condition.
#' 
#' @return A data.frame containing locustag (gene_names), 
#' contrasts and M-values for the current organism and genes.
#'
#' @references http://colombos.net
#'
#' @export
#' 
advanced_search_by_contrasts <- function(organism=NULL, ids=NULL, geneNames=FALSE, c_search_type=NULL){
    if(is.null(ids)) stop("Insert the ids for the specific search_type.") else {}
    if(c_search_type=="contrast_names"){
        url_string <- paste(getOption("REST.version"),"advanced_search_by_contrast", organism, "contrast_names", paste(ids, collapse=","), sep="/")
    }
    else if(c_search_type=="experiment"){
      url_string <- paste(getOption("REST.version"),"advanced_search_by_contrast", organism, "experiment", gsub(" ","%20", ids), sep="/")
    }
    else if(c_search_type=="go"){
        url_string <- paste(getOption("REST.version"),"advanced_search_by_contrast", organism, "go", gsub(" ","%20", ids), sep="/")
    }
    else if(c_search_type=="condition"){
        url_string <- paste(getOption("REST.version"),"advanced_search_by_contrast", organism, "condition", paste(ids, collapse=","), sep="/")
    } else {
        stop("Wrong search_type: it should be either contrast_names, experiment, go or condition!")
    }
    r <- GET(url_string)
    if (r$status_code != 200) {
        stop_for_status(r)
    } else {
        tmp <- content(r)
        response = suppressWarnings( data.frame(matrix(as.numeric(sapply(tmp$data$Mvalues, as.character)), 
        nrow=length(tmp$data$geneLocustags), ncol=length(tmp$data$contrasts), byrow=T)) )
        colnames(response) <- tmp$data$contrasts; rownames(response) <- tmp$data$geneLocustags
        if (geneNames) rownames(response) <- tmp$data$geneNames
        return(response)
    }
}
#' Accessory function allowing the advanced_search by both g_ids and c_ids
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param g_ids A vector of strings representing contrast_id, go terms, experiment id or condition id according the search type.
#' @param geneNames boolean if == FALSE (default) return the locustag otherwise the gene_name for the selected genes.
#' @param c_ids A vector of strings representing contrast_id, go terms, experiment id or condition id according the search type.
#' @param g_search_type A string either genes, go or annotation.
#' @param ann_type A string containing the selected gene_annotation_type: use \code{\link{listEntities}} to display the available entities.
#' @param c_search_type A string either contrast_names. experiment, go, condition.
#' 
#' @return A data.frame containing locustag (gene_names), 
#' contrasts and M-values for the current organism and genes.
#'
#' @references http://colombos.net
#'
#' @export
#' 
advanced_search_by_both <- function(organism, g_ids, geneNames, c_ids, g_search_type, ann_type, c_search_type){
    if(is.null(g_ids)) stop("Insert a character vector with the g_ids to be imputed.") else {}
    if(is.null(c_ids)) stop("Insert a character vector with the c_ids to be imputed.") else {}
    base_url <- paste(getOption("REST.version"),"advanced_search_by_both", sep="/")
    if(g_search_type=="genes"){
        g_string <- paste(base_url, organism, "genes", paste(g_ids, collapse=","), sep="/")
    }
    else if(g_search_type=="go"){
        g_string <- paste(base_url, organism, "go", gsub(" ","%20", g_ids), sep="/")
    }
    else if(g_search_type=="annotation"){
        g_string <- paste(base_url, organism, "annotation", gsub(" ","%20", ann_type), 
                          gsub(" ","%20", g_ids), sep="/")  
    } else {
        stop("Wrong search_type: it should be either genes, go or annotation!")
    }
    #
    if(c_search_type=="contrast_names"){
        url_string <- paste(g_string, "contrast_names", paste(c_ids, collapse=","), sep="/")
    }
    else if(c_search_type=="experiment"){
        url_string <- paste(g_string, "experiment", gsub(" ","%20", c_ids), sep="/")
    }
    else if(c_search_type=="go"){
        url_string <- paste(g_string, "go", gsub(" ","%20", c_ids), sep="/")
    }
    else if(c_search_type=="condition"){
        url_string <- paste(g_string, "condition", paste(c_ids, collapse=","), sep="/")
    } else {
        stop("Wrong search_type: it should be either contrast_names, experiment, go or condition!")
    }
    r <- GET(url_string)
    if (r$status_code != 200) {
        stop_for_status(r)
    } else {
        tmp <- content(r)
        response = suppressWarnings( data.frame(matrix(as.numeric(sapply(tmp$data$Mvalues, as.character)), 
        nrow=length(tmp$data$geneLocustags), ncol=length(tmp$data$contrasts), byrow=T)) )
        colnames(response) <- tmp$data$contrasts; rownames(response) <- tmp$data$geneLocustags
        if (geneNames) rownames(response) <- tmp$data$geneNames
        return(response)
    }
}
