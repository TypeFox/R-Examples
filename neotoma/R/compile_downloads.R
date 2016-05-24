#' Function to convert multiple downloads into a single large table.
#'
#' From the assemblage data for multiple cores, return a single data.frame with columns for site
#' metadata and assemblage data.
#'
#' To support further synthesis and analysis \code{compile_download} works to transform a list
#' returned by \code{\link{get_download}} into a large data frame with columns for site and sample attributes
#' and also with the associated assemblage data at each sample depth.  This function also does the same for
#' single sites.
#'
#' @importFrom plyr ldply
#' @param downloads A download_list as returned by \code{\link{get_download}}, or mutliple downloads joined in a list.
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return This command returns a data frame.
#'
#' @examples \dontrun{
#' #  Search for sites with "Thuja" pollen that are older than 8kyr BP and
#' #  that are on the west coast of North America:
#' t8kyr.datasets <- get_dataset(taxonname='Thuja*', loc=c(-150, 20, -100, 60), ageyoung = 8000)
#'
#' #  Returns 3 records (as of 04/04/2013), get dataset for the first record, Gold Lake Bog.
#' thuja.sites <- get_download(t8kyr.datasets)
#'
#' gold.p25 <- compile_taxa(thuja.sites, 'P25')
#'
#' all.gold <- compile_downloads(gold.p25)
#' 
#' pollen.sums <- rowSums(all.gold[,11:ncol(all.gold)], na.rm=TRUE)
#' 
#' plot(x = all.gold$age, 
#'      y = all.gold$Cupressaceae.Taxaceae / pollen.sums, 
#'      col = all.gold$site.name,
#'      pch = 19)
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

compile_downloads <- function(downloads) {

  if (!('download_list' %in% class(downloads) | all(sapply(downloads, class) == 'download') |
       'download' %in% class(downloads))) {
    stop('compile_datasets can only operate on lists as returned from get_download, or a list of download objects.')
  }

  down.to.df <- function(x) {
    if ('download' %in% class(x)) {
      #  There can be NULL values in the download object.  We'll turn them to NA values:
      if (is.null(x$dataset$site$site.name)) x$dataset$site$site.name <- paste('NoName_ID')
      if (is.null(x$sample.meta$depths)) x$sample.meta$depths <- NA
      if (is.null(x$sample.meta$age)) x$sample.meta$age <- NA
      if (is.null(x$sample.meta$age.older)) x$sample.meta$age.older <- NA
      if (is.null(x$sample.meta$age.younger)) x$sample.meta$age.younger <- NA
      if (is.null(x$dataset$site$lat)) x$dataset$site$lat <- NA
      if (is.null(x$dataset$site$long)) x$dataset$site$long <- NA

      site.info <- data.frame(site.name = x$dataset$site$site.name,
                              depth = x$sample.meta$depth,
                              age = x$sample.meta$age,
                              age.old = x$sample.meta$age.older,
                              age.young = x$sample.meta$age.younger,
                              date.type = x$sample.meta$age.type,
                              lat = x$dataset$site$lat,
                              long = x$dataset$site$long,
                              dataset = x$dataset$dataset.meta$dataset.id,
                              x$counts)
    }
    else{
      #  Dummy data for empty sites.
      site.info <- data.frame(site.name = NA,
                              depth = NA,
                              age = NA,
                              age.old = NA,
                              ageyoung = NA,
                              date.type = NA,
                              lat = NA,
                              long = NA,
                              dataset = NA,
                              unknown = NA)
    }
    site.info
  }


  if (inherits(downloads, 'download_list')) {
      site.info <- plyr::ldply(downloads, down.to.df)
  }
  if (inherits(downloads, 'download')) {
      site.info <- down.to.df(downloads)
  }

  site.info
}
