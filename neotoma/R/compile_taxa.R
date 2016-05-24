
#' Function to convert assemblage taxa to standardized lists.
#'
#' From the assemblage data for the core return assemblage data with the assemblage taxa
#' Currently implemented only for pollen data.
#'
#' The data object uses the smaller pollen subset.  As this package develops we will add the capacity to summarize data output from the translation. Currently we can return only subsets that have been defined in the literature.  These lists include:
#' \itemize{
#'   \item{\code{"P25"} }{ This list is derived from Gavin et al., (2003), and includes 25 pollen taxa.}
#'   \item{\code{"WS64"} }{  This list is derived from Williams and Shuman (2008).}
#'   \item{\code{"WhitmoreFull"} }{  This is the full list associated with the Whitmore et al., (2005) North American Modern Pollen Database.}
#'   \item{\code{"WhitmoreSmall"} }{  As above, but taxa for which both fully resolved and undifferentiated exist these taxa are summed.}
#' }
#'
#' @importFrom stats aggregate
#' @param object A pollen object returned by \code{\link{get_download}}.
#' @param alt.table A user provided table formatted with at least two columns, one called 'taxon' and the other named as in \code{list.name}.
#' @param list.name The taxon compilation list, one of a set of lists from the literature (e.g., \code{"P25"}, \code{"WhitmoreFull"}).  More detail in section Details.
#' @param cf Should taxa listed as *cf*s (*e.g.*, *cf*. *Gilia*) be considered highly resolved?
#' @param type Should taxa listed as types (*e.g.*, *Iva annua*-type) be considered highly resolved?
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return This command returns a list object with the same structure as the parent pollen object returned by \code{\link{get_download}}, or a matrix (or data frame) depending on whether \code{object} is one or the other.  Any pollen taxon not included in the major taxa defined in the pollen gets returned as 'Other'.
#'
#' @examples \dontrun{
#' #  Search for sites with "Thuja" pollen that are older than 8kyr BP and
#' #  that are on the west coast of North America:
#' t8kyr.datasets <- get_dataset(taxonname='Thuja*', loc=c(-150, 20, -100, 60), ageyoung = 8000)
#'
#' #  Returns 3 records (as of 04/04/2013), get dataset for the first record, Gold Lake Bog.
#' GOLDKBG <- get_download(t8kyr.datasets[[1]])
#'
#' gold.p25 <- compile_taxa(GOLDKBG, 'P25')
#'
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#'
#' Gavin DG, Oswald WW, Wahl ER, Williams JW. 2003. A statistical approach to evaluating distance metrics and analog assignments for pollen records. Quaternary Research 60: 356-367.
#'
#' Whitmore J, Gajewski K, Sawada M, Williams JW, Shuman B, Bartlein PJ, Minckley T, Viau AE, Webb III T, Shafer S, Anderson P, Brubaker L. 2005. Modern pollen data from North America and Greenland for multi-scale paleoenvironmental applications. Quaternary Science Reviews 24: 1828-1848.
#'
#' Williams J, Shuman B. 2008. Obtaining accurate and precise environmental reconstructions from the modern analog technique and North American surface pollen dataset. Quaternary Science Reviews. 27:669-687.
#'
#' API Reference:  http://api.neotomadb.org/doc/resources/contacts
#' @keywords utilities
#' @export

compile_taxa <- function(object, list.name, alt.table = NULL, cf = TRUE, type = TRUE){

    if (!inherits(object, c('matrix', 'data.frame', 'download', 'download_list'))) {
        stop(paste('Data object must be a pollen object returned by',
                   'function get_download or a matrix or data frame'))
    }

    if(!is.null(alt.table)) {
        if(!inherits(alt.table, c("matrix", "data.frame"))) {
            stop('The alt.table must be either a matrix or a data frame.')
        }

        pollen.equiv <- alt.table
        avail.lists <- colnames(pollen.equiv)
        if(!list.name %in% avail.lists){
            stop('The list name is not included in your alt.table.')
        }

        if(!'taxon' %in% (avail.lists)){
            stop('The alt.table must contain a column titled taxon.')
        }
        use.list <- which(avail.lists %in% list.name)
    } else {
        avail.lists <- c('P25', 'WS64', 'WhitmoreFull', 'WhitmoreSmall')
        use.list <- which(avail.lists %in% list.name) + 2
    }

    if (cf == FALSE)   list.name <- list.name[is.na(pollen.equiv$cf)]
    if (type == FALSE) list.name <- list.name[is.na(pollen.equiv$type)]

    if (inherits(object, c("download","download_list"))) {

        aggregate.counts <- function(x){
            taxon.matches <- match(colnames(x$counts), pollen.equiv$taxon)

            used.taxa <- pollen.equiv[taxon.matches, ]
            agg.list <- as.vector(used.taxa[, use.list])

            agg.list[is.na(agg.list)] <- 'Other'

            compressed.list <- stats::aggregate(t(x$counts), by = list(agg.list),
                                         sum, na.rm = TRUE)

            compressed.cols <- compressed.list[, 1]

            compressed.list <- t(compressed.list[, -1])
            colnames(compressed.list) <- compressed.cols

            ## We want to make a taxon list like the one returned in get_downloads:
            new.list <- x$taxon.list
            new.list$compressed <- NA

            new.list$compressed <- as.character(pollen.equiv[match(new.list$taxon.name, pollen.equiv$taxon),
                                                             use.list])

            new.list$compressed[is.na(new.list$compressed) & new.list$taxon.name %in% colnames(x$counts)] <- 'Other'
            new.list <- new.list[match(new.list$taxon.name, x$taxon.list$taxon.name),]

            ## Returns a data.frame with taxa in the columns and samples in the rows.
            output <- list(dataset = x$dataset,
                           sample.meta = x$sample.meta,
                           taxon.list = new.list,
                           counts = compressed.list,
                           full.counts = x$counts,
                           lab.data = x$lab.data,
                           chronologies = x$chronologies)

            missed.samples <- as.character(unique(new.list$taxon.name[which(new.list$compressed == 'Other')]))

            if (length(missed.samples)>0){
              warning(paste0('\nThe following taxa could not be found in the existing ',
                             'conversion table:\n', paste(missed.samples, sep = '\n')))
            }

            class(output) <- c('download', 'list')

            output

        }

        if (inherits(object, "download_list")) {
            output <- lapply(object, FUN = aggregate.counts)
            class(output) <- c('download_list', 'list')
        } else {
            output <- aggregate.counts(object)
            class(output) <- c('download', 'list')
        }
    }

    if (inherits(object, c("matrix", "data.frame"))) {
        taxon.matches <- match(colnames(object), pollen.equiv$taxon)
        if (any(is.na(taxon.matches))){
            missed.samples <- colnames(object)[is.na(taxon.matches)]
            warning(paste0('\nThe following taxa could not be found in the existing ',
                           'conversion table:\n', paste(missed.samples, sep = '\n')))
        }

        used.taxa <- pollen.equiv[taxon.matches, ]
        agg.list <- as.vector(used.taxa[, use.list])
        agg.list[is.na(agg.list)] <- 'Other'

        compressed.list <- stats::aggregate(t(object), by = list(agg.list),
                                     sum, na.rm = TRUE)

        compressed.cols <- compressed.list[, 1]

        compressed.list <- t(compressed.list[, -1])
        colnames(compressed.list) <- compressed.cols

        output <- compressed.list

    }
    output
}
