#' Summarize chromosome counts from API call
#'
#' This function processes and cleans the data returned from the API call for use in downstream analysis.
#'
#' @param counts A 'chrom.counts' object inherited from \code{\link{chrom_counts}}.
#'
#' @details The results from the API call are a bit messy and difficult to use for downstream analyses. This function cleans up the data in three ways. First, it combines aggregates and summarizes all records from each species. Second, many of the counts are combined with text characters (e.g., \code{"#-#"}, \code{"c.#"}, and \code{"#, #, #"}. This function uses regular expressions to pull out all and any numeric values from these strings. Third, some of the records are gametophytic (n) counts and others are from sporophytes (2n); the function simply divides the sporophytic counts in half so that all measurements are on a common scale.
#'
#' IMPORTANT: Use this function with caution. Parsing the counts programmatically may be useful but it may generate erroneous results in some cases if input is in an odd format. For example, if the count is \code{"#+-#"}, the function will return both the first and second \code{#} as valid counts . Given the creativity(?) of researchers in entering data, it is hard to predict all possible ways that the counts may be represented. Therefore, some manual checking will probably be necessary.
#'
#' @return A \code{data.frame} containing the resolved binomial, the count type (gametophytic or sporophytic), the counts, the inferred gametophytic count (for sporophytic records) and the number of records supporting each count.
#'
#' @import dplyr
#' @importFrom data.table rbindlist
#'
#' @export summarize_counts
#'
#' @examples \dontrun{
#'
#' ## Get all counts for genus Castilleja
#' res <- chrom_counts("Castilleja", "genus")
#'
#' ## summarize the results
#' summarize_counts(res)
#'
#' }
summarize_counts <- function(counts){

    if (!inherits(counts, "chrom.counts"))
        stop("Object must be of class 'chrom.counts' returned from chrom_counts()")
    sp <- unique(counts$resolved_binomial)
    df <- lapply(sp, function(x) filter_(counts, ~(resolved_binomial == x)))
    ct <- lapply(df, get_counts)
    tbl_df(rbindlist(ct))
}


## Function for counting sporophytic and gametophytic numbers for individual spp.
get_counts <- function(x){

    cnt <- lapply(c("gametophytic", "sporophytic"),
                  function(y) get_counts_n(x,y))
    df <- data.frame(rbindlist(cnt))
    df$resolved_binomial <- unique(x$resolved_binomial)

    ## reorder
    ord <- c("resolved_binomial", "count_type", "count",
             "inferred_n", "num_records")
    df[,ord]
}

    
## internal function
get_counts_n <- function(x, type){
    cln <- parse_counts(x[,type])
    uni <- sort(unique(cln))
    num <- sapply(uni, function(y) length(cln[cln == y]))

    if (type == "gametophytic"){
        ploidy <- rep("gam", length(uni))
        inf_n  <- uni
    } else {
        ploidy <- rep("spor", length(uni))
        inf_n  <- round(uni/2)
    }

    if (length(num) == 0){
        return(data.frame())
    } else {
        return(data.frame(count_type=ploidy,count=uni,inferred_n=inf_n,
                          num_records=num))
    }
}


parse_counts <- function(x){
    tmp <- na.omit(as.numeric(unlist(strsplit(unlist(x), "[^0-9]+"))))
    tmp[tmp != 0]
}

    

    
    

    
