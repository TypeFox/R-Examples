#' Function to view enrichment results in a sample-specific manner
#'
#' \code{dGSEAview} is supposed to view results of gene set enrichment analysis but for a specific sample. 
#'
#' @param eTerm an object of class "eTerm"
#' @param which_sample which sample will be viewed
#' @param top_num the maximum number of gene sets will be viewed
#' @param sortBy which statistics will be used for sorting and viewing gene sets. It can be "adjp" for adjusted p value, "gadjp" for globally adjusted p value, "ES" for enrichment score, "nES" for normalised enrichment score, "pvalue" for p value, "FWER" for family-wise error rate, "FDR" for false discovery rate, "qvalue" for q value
#' @param decreasing logical to indicate whether to sort in a decreasing order. If it is null, it would be true for "ES" or "nES"; otherwise it would be false
#' @param details logical to indicate whether the detail information of gene sets is also viewed. By default, it sets to false for no inclusion
#' @return
#' a data frame with following components:
#' \itemize{
#'  \item{\code{setID}: term ID}
#'  \item{\code{ES}: enrichment score}
#'  \item{\code{nES}: normalised enrichment score}
#'  \item{\code{pvalue}: nominal p value}
#'  \item{\code{adjp}: adjusted p value}
#'  \item{\code{gadjp}: globally adjusted p value}
#'  \item{\code{FDR}: false discovery rate}
#'  \item{\code{qvalue}: q value}
#'  \item{\code{setSize}: the number of genes in the set; optional, it is only appended when "details" is true}
#'  \item{\code{name}: term name; optional, it is only appended when "details" is true}
#'  \item{\code{namespace}: term namespace; optional, it is only appended when "details" is true}
#'  \item{\code{distance}: term distance; optional, it is only appended when "details" is true}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dGSEA}}
#' @include dGSEAview.r
#' @examples
#' #dGSEAview(eTerm, which_sample=1, top_num=10, sortBy="adjp", decreasing=FALSE, details=TRUE)

dGSEAview <- function(eTerm, which_sample=1, top_num=10, sortBy=c("adjp","gadjp","ES","nES","pvalue","FWER","FDR","qvalue"), decreasing=NULL, details=F) 
{

    if (class(eTerm) != "eTerm" ){
        stop("The function must apply to a 'eTerm' object.\n")
    }
    
    sortBy <- match.arg(sortBy)
    
    if( is.null(top_num) ){
        top_num <- length(eTerm$set_info$setID)
    }
    if( top_num > length(eTerm$set_info$setID) ){
        top_num <- length(eTerm$set_info$setID)
    }
    
    which_sample <- as.integer(which_sample)
    if(which_sample > ncol(eTerm$adjp)){
        which_sample <- ncol(eTerm$adjp)
    }else if(which_sample < 1){
        which_sample <- 1
    }
    
    gs_size <- sapply(1:length(eTerm$gs), function(x){
        length(eTerm$gs[[x]])
    })
    
    tab <- data.frame(setID         = eTerm$set_info$setID,
                       ES           = eTerm$es[,which_sample],
                       nES          = eTerm$nes[,which_sample],
                       pvalue       = eTerm$pvalue[,which_sample],
                       adjp         = eTerm$adjp[,which_sample],
                       gadjp        = eTerm$gadjp[,which_sample],
                       FWER         = eTerm$fwer[,which_sample],
                       FDR          = eTerm$fdr[,which_sample],
                       qvalue       = eTerm$qvalue[,which_sample],
                       setSize      = gs_size,
                       name         = eTerm$set_info$name,
                       namespace    = eTerm$set_info$namespace,
                       distance     = eTerm$set_info$distance
                      )
    
    
    if(details == T){
        res <- tab
    }else{
        res <- tab[,c(1:9)]
    }
    
    if(is.null(decreasing)){
        if(sortBy=="nES" | sortBy=="ES"){
            decreasing <- T
        }else{
            decreasing <- F
        }
    }
    
    switch(sortBy, 
        ES={res <- res[order(res[,2], decreasing=decreasing)[1:top_num],]},
        nES={res <- res[order(res[,3], decreasing=decreasing)[1:top_num],]},
        pvalue={res <- res[order(res[,4], decreasing=decreasing)[1:top_num],]},
        adjp={res <- res[order(res[,5], decreasing=decreasing)[1:top_num],]},
        gadjp={res <- res[order(res[,6], decreasing=decreasing)[1:top_num],]},
        FWER={res <- res[order(res[,7], decreasing=decreasing)[1:top_num],]},
        FDR={res <- res[order(res[,8], decreasing=decreasing)[1:top_num],]},
        qvalue={res <- res[order(res[,9], decreasing=decreasing)[1:top_num],]}
    )
    
    res
}