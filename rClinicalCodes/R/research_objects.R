#' Associates lists of clinical codes with a ClinicalCodes research object
#' 
#' This function uses the codelist_url slot in a Research object to download codelists
#' 
#' Internally it calls the get_clinicalcodes() function with a url argument
#' @param ro an object of class ResearchObject
#' @export
#' @return a new ResearchObject with associated clinicalcodes dataframes for all codelists
#' @examples \dontrun{
#' # get research object without codes
#' RO <- research_object(article_id = 2, download_codes = FALSE)
#' # associate codes
#' RO <- codelists_for_ro(RO)
#' } 
codelists_for_ro <- function(ro){
    if(class(ro) == "ResearchObject"){
        for(codelist in names(ro$codelists)){
            message(sprintf("Downloading codelist for %s", codelist))
            ro$codelists[[codelist]]$clinicalcodes <- get_ClinicalCodes(url = ro$codelists[[codelist]]$codelist_url)
        }
        ro
    } else stop("Not a valid ResearchObject")
}


#' Downloads Research objects for articles from www.clinicalcodes.org
#' 
#' This function builds an R representation of a ClinicalCodes Research object
#' 
#' The ResearchObject contains metadata describing the article (URI, abstract, ID, title, authors, doi, journal etc.),
#' comments on the article, codelist metadata (associated articles, name, url, number of codes in the list, user field names, comments) 
#' and optional full codelists.
#' The individual codelists can be downloaded directly with the download_codes argument set to TRUE.
#' Otherwise they can be associated later using the codelists_for_ro() function
#' 
#' @param article_ids integer representiong the id of a www.ClinicalCodes.org article
#' @param download_codes logical should the individual code lists be downloaded as part of the research object?
#' @param trim integer How many characters from titles should be included in list names for multiple codelists?
#' @export
#' @return an object of class ResearchObject (R representation of a ClinicalCodes research object)
#' @examples \dontrun{
#' # get research object and codes
#' RO <- research_object(article_ids = 2, download_codes = FALSE)
#' ROs  <- research_object(article_ids = 5:7, download_codes = FALSE)
#' } 
research_object <- function(article_ids, download_codes = FALSE, trim = 50){
    if(length(article_ids) > 1){
        ROs_ <- lapply(article_ids, function(id){
            message(sprintf("Research object for article %d", id))
            ro_url <- sprintf("https://clinicalcodes.rss.mhs.man.ac.uk/medcodes/article/%d/ro", id)
            tryCatch(ro_ <- fromJSON(getURL(ro_url, ssl.verifypeer = 0L, followlocation = 1L)),
                     error = function(e){
                         warning(sprintf("Article %d not found in the repository", id))
                         ro_  <<- NULL
                     })
            if(!is.null(ro_)){
                ro_ <- structure(ro_, class = "ResearchObject")
                if(download_codes) ro_ <- codelists_for_ro(ro_)
            }
            Sys.sleep(1) # going easy on my server!
            ro_
        })
        ROs_  <- ROs_[!sapply(ROs_, is.null)]
        names(ROs_) <- sapply(ROs_, function(x) paste(x$article_ID, strtrim(x$article_title, trim), sep = ": "))
        ROs_
    } else{
        ro_url <- sprintf("https://clinicalcodes.rss.mhs.man.ac.uk/medcodes/article/%d/ro", article_ids)
        tryCatch(ro_ <- fromJSON(getURL(ro_url, ssl.verifypeer = 0L, followlocation = 1L)),
                 error = function(e) stop(sprintf("Article %d not found in the repository", article_ids)))
        ro_ <- structure(ro_, class = "ResearchObject")
        if(download_codes) ro_ <- codelists_for_ro(ro_)
        ro_
    } 
}




