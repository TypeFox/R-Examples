#' Function to view enrichment results of dEnricher
#'
#' \code{dEnricherView} is supposed to view results of enrichment analysis by \code{\link{dEnricher}}. 
#'
#' @param eTerm an object of class "eTerm"
#' @param top_num the maximum number of gene sets (terms) will be viewed
#' @param sortBy which statistics will be used for sorting and viewing gene sets (terms). It can be "adjp" for adjusted p value, "pvalue" for p value, "zscore" for enrichment z-score, "nAnno" for the number of sets (terms), "nOverlap" for the number in overlaps, and "none" for ordering according to ID of gene sets (terms)
#' @param decreasing logical to indicate whether to sort in a decreasing order. If it is null, it would be true for "zscore", "nAnno" or "nOverlap"; otherwise it would be false
#' @param details logical to indicate whether the detailed information of gene sets (terms) is also viewed. By default, it sets to false for no inclusion
#' @return
#' a data frame with following components:
#' \itemize{
#'  \item{\code{setID}: term ID; as rownames}
#'  \item{\code{name}: term name}
#'  \item{\code{nAnno}: number in gene members annotated by a term}
#'  \item{\code{nOverlap}: number in overlaps}
#'  \item{\code{zscore}: enrichment z-score}
#'  \item{\code{pvalue}: nominal p value}
#'  \item{\code{adjp}: adjusted p value}
#'  \item{\code{namespace}: term namespace; optional, it is only appended when "details" is true}
#'  \item{\code{distance}: term distance; optional, it is only appended when "details" is true}
#'  \item{\code{members}: members (represented as Gene Symbols) in overlaps; optional, it is only appended when "details" is true}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dEnricher}}
#' @include dEnricherView.r
#' @examples
#' #dEnricherView(eTerm, top_num=10, sortBy="adjp", decreasing=FALSE, details=TRUE)

dEnricherView <- function(eTerm, top_num=10, sortBy=c("adjp","pvalue","zscore","nAnno","nOverlap","none"), decreasing=NULL, details=F) 
{
    
    if(is.logical(eTerm)){
        stop("There is no enrichment in the 'eTerm' object.\n")
    }
    
    if (class(eTerm) != "eTerm" ){
        stop("The function must apply to a 'eTerm' object.\n")
    }
    
    sortBy <- match.arg(sortBy)
    
    if( is.null(top_num) ){
        top_num <- length(eTerm$set_info$setID)
    }
    if ( top_num > length(eTerm$set_info$setID) ){
        top_num <- length(eTerm$set_info$setID)
    }
    
    if(dim(eTerm$set_info)[1]==1){
        tab <- data.frame( name         = eTerm$set_info$name,
                           nAnno         = sapply(eTerm$gs,length),
                           nOverlap     = sapply(eTerm$overlap,length),
                           zscore       = eTerm$zscore,
                           pvalue       = eTerm$pvalue,
                           adjp         = eTerm$adjp,
                           namespace    = eTerm$set_info$namespace,
                           distance     = eTerm$set_info$distance,
                           members      = sapply(eTerm$overlap, function(x) paste(names(x),collapse=','))
                          )
    }else{
    
        tab <- data.frame( name         = eTerm$set_info$name,
                           nAnno         = sapply(eTerm$gs,length),
                           nOverlap     = sapply(eTerm$overlap,length),
                           zscore       = eTerm$zscore,
                           pvalue       = eTerm$pvalue,
                           adjp         = eTerm$adjp,
                           namespace    = eTerm$set_info$namespace,
                           distance     = eTerm$set_info$distance,
                           members      = sapply(eTerm$overlap, function(x) paste(names(x),collapse=','))
                          )
    }
    
    rownames(tab) <- eTerm$set_info$setID
    
    if(details == T){
        res <- tab[,c(1:9)]
    }else{
        res <- tab[,c(1:6)]
    }
    
    if(is.null(decreasing)){
        if(sortBy=="zscore" | sortBy=="nAnno" | sortBy=="nOverlap"){
            decreasing <- T
        }else{
            decreasing <- F
        }
    }
    
    switch(sortBy, 
        adjp={res <- res[order(res[,6], decreasing=decreasing)[1:top_num],]},
        pvalue={res <- res[order(res[,5], decreasing=decreasing)[1:top_num],]},
        zscore={res <- res[order(res[,4], decreasing=decreasing)[1:top_num],]},
        nAnno={res <- res[order(res[,2], decreasing=decreasing)[1:top_num],]},
        nOverlap={res <- res[order(res[,3], decreasing=decreasing)[1:top_num],]},
        none={res <- res[order(rownames(res), decreasing=decreasing)[1:top_num],]}
    )
    
    if(sortBy=='none'){
    	suppressWarnings(flag <- all(!is.na(as.numeric(rownames(res)))))
    	if(flag){
    		res <- res[order(as.numeric(rownames(res)), decreasing=decreasing)[1:top_num],]
    	}
    }
    
    res
}