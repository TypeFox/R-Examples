#' Downloads a list of all articles in ClinicalCodes
#' 
#' This function gives a dataframe of all articles listed in the ClinicalCodes repository. 
#' 
#' Article IDs from the list are included and can then be passed to the download functions
#' to import codelists and research objects into R
#' @export
#' @return A dataframe of article information, links and IDs
#' @examples \dontrun{
#' all_articles <- all_ClinicalCodes_articles()
#' }
all_ClinicalCodes_articles <- function(){
    base_url <- "https://clinicalcodes.rss.mhs.man.ac.uk"
    url <- getURL(paste0(base_url, "/medcodes/articles/"), ssl.verifypeer = FALSE)
    doc <- htmlParse(url)
    article_table <- getNodeSet(doc, "//table")[[1]]
    table_data <- readHTMLTable(article_table)
    links <- unlist(xpathSApply(article_table, "//a", xmlGetAttr, name = "href"))
    links <- links[str_detect(links, "/medcodes/article/")]
    table_data$link <- paste0(base_url, links)
    table_data$ID <- as.integer(str_match(links, "([0-9]+)/$")[,2])
    table_data
}

