#' Returns chromosome counts from Chromosome Counts Database API
#'
#' This function calls the Chromsome Counts Database (CCDB) API and returns all counts for specified higher taxa.
#'
#' @param taxa Taxonomic name(s) to query. Can be either a single name, a vector of multiple names or a list. If supplying multiple names, these must all be of the same \code{rank}.
#'
#' @param rank Rank to query.
#'
#' @param full Logical. Whether to return full records. Defaults to \code{FALSE} which returns only partial records. Partial records includes the resolved name as well as the gametophytic (n) and sporophytic (2n) counts.
#'
#' @param foptions additional options to be passed to \code{httr::GET}
#'
#' @details When using the API to query for species, both matched names and resolved names are searched. This means that all records for potential synonyms will be returned as well. Currently species binomials must be specified by either 'genus species' (i.e., space between genus and species) or 'genus_species'.
#'
#' To search for subspecies (subsp.) or varieties (var.), you can use search terms like:
#' 
#' \code{"Solanum acaule var. albicans"}.
#'
#' Searching for \code{"Solanum acaule"} will return all subspecies and varieties.
#'
#' Currently the only acceptable search terms when specifying \code{"majorGroup"} are \code{"Angiosperms"}, \code{"Gymnosperms"}, \code{"Pteridophytes"}, or \code{"Bryophytes"}.
#'
#' @return A \code{data.frame} containing all records matched by query
#'
#' @import dplyr
#' @importFrom httr GET content stop_for_status
#' @importFrom data.table rbindlist
#' 
#' @export chrom_counts
#'
#' @examples \dontrun{
#' 
#' ## Get all counts for genus Castilleja
#' chrom_counts("Castilleja", "genus")
#'
#' ## Get all counts for both Castilleja and Lachemilla
#' chrom_counts(c("Castilleja", "Lachemilla"), "genus")
#'
#' ## Get all counts for Castilleja miniata
#' chrom_counts("Castilleja miniata", "species")
#'
#' ## Get all counts for only Castilleja miniata subsp. elata
#' chrom_counts("Castilleja miniata subsp. elata", "species")
#'
#' ## Note that searching for "Castilleja miniata" will return
#' ## all subspecies and varieties
#'
#' ## Get all counts for the Orobanchaceae
#' chrom_counts("Orobanchaceae", "family")
#' 
#' }
chrom_counts <-  function(taxa,
                          rank=c("species", "genus", "family", "majorGroup"),
                          full=FALSE, foptions=list()){
    
    out <- suppressWarnings(check_ccdb_input(rank, full)) 
    l   <- lapply(taxa, function(x)
                chrom_counts_single(x, rank, out, foptions=foptions))
    res <- tbl_df(rbindlist(l))
    res <- tidy_output(res)
    attr(res, "class") <- c(attr(res,"class"), "chrom.counts")
    res
}

## Internal function to do the individual queries
chrom_counts_single <- function(taxa, rank, out, foptions){

    if (rank == "species")
        taxa <- species_API(taxa)
    
    url <- paste0("http://ccdb.tau.ac.il/services/",
                  out,"/?", rank,"=",taxa,"&format=","json")
    counts_call <- GET(url, foptions)
    stop_for_status(counts_call)
    counts_data_json <- content(counts_call)

    f <- function(x) if (is.null(x)) NA_character_ else x
    counts_data_json <- lapply(counts_data_json, lapply, f)
    counts_data <- data.frame(rbindlist(counts_data_json))
    
    if (length(counts_data) > 0)
        counts_data <- add_binomial(counts_data)

    counts_data
}


## Utility function for checking input
check_ccdb_input <- function(rank, full){
        
    if (length(rank) != 1 | !rank %in% rank_names())
        stop("Specify a single taxonomic rank. \n Options are 'species', 'genus', 'family', and 'majorGroup'.")

    if (full){
        output <- "countsFull"
    } else {
        output <- "countsPartial"
    }
    output
}


## little function returning acceptable rank names
rank_names <- function()
    c("species", "genus", "family", "majorGroup")


## function for removing species from API call
## otherwise API query stalls
species_API <- function(x)
    gsub(" ", "%20", x)


## Function for pulling out species name without the authorities
## Keeping varieties and subspecies
## These are indicated by var. and subsp., respectively
add_binomial <- function(x)
    x %>% rowwise() %>%
    mutate_(resolved_binomial = ~short_species_name(resolved_name))


short_species_name <- function(x){
    tmp <- strsplit(x, split=" ")[[1]]
    ## keep varities and subspecies
    ## Currently depending on this structure
    if ("var." %in% tmp | "subsp." %in% tmp){
        sp <- paste(tmp[1],tmp[2],tmp[3],tmp[4],sep="_")
    } else {
        sp <- paste(tmp[1],tmp[2],sep="_")
    }
    sp
}


tidy_output <- function(x){
    ord_partial <- c("resolved_binomial", "gametophytic", "sporophytic",
                     "resolved_name")

    ord_full    <- c("resolved_binomial", "gametophytic", "sporophytic",
                     "resolved_name", "original_name", "matched_name",
                     "taxonomic_status", "genus", "family", "major_group",
                     "id", "source", "internal_id", "reference", "voucher")

    if (ncol(x) == 0){
        return(x)
    } else if (ncol(x) == 15){
        return(x[,ord_full])
    } else if (ncol(x) == 4){
        return(x[,ord_partial])
    }

}
    
