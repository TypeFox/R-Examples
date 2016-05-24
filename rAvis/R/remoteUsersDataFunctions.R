#' avisContributorsSummary
#' 
#'Returns a table with the 
#'observations aggregated by birdwatcher. 
#'
#' @usage avisContributorsSummary()
#' @note This function does not allow arguments
#' @return This function returns a matrix
#' @export 
#' @examples \dontrun{
#' birdwatchers<- avisContributorsSummary()
#' par (mfrow =c(2,2))
#' plot (birdwatchers[,2],birdwatchers[,3], xlab=colnames (birdwatchers)[2], 
#' ylab=colnames (birdwatchers)[3], pch=19)
#' plot (birdwatchers[,2],birdwatchers[,4], xlab=colnames (birdwatchers)[2], 
#' ylab=colnames (birdwatchers)[4], pch=19)
#' plot (birdwatchers[,2],birdwatchers[,5], xlab=colnames (birdwatchers)[2], 
#' ylab=colnames (birdwatchers)[5], pch=19)
#' plot (birdwatchers[,2],birdwatchers[,6], xlab=colnames (birdwatchers)[2], 
#' ylab=colnames (birdwatchers)[6], pch=19)
#' }

avisContributorsSummary<- function ()
{
  cs <- .avisCacheReturnOrSetup("contributors_summary", ".avisContributorsSummaryInit")

  return(cs)
}

#' avisContributorAggregatedObservations
#' 
#' A function to download the information about the
#'  observations of a birdwatcher.
#'
#' @usage avisContributorAggregatedObservations(contributor_id)
#' @param contributor_id a number setting the id of the birdwatcher (see avisContributorSummary)
#' @return This function returns a dataframe
#' @export 
#' @examples # Explore the contributions of Colectivo Ornitologico Ciguena Negra
#' \dontrun{
#' avisContributorAggregatedObservations (370)
#' }
#' 
#' 
avisContributorAggregatedObservations <- function (contributor_id)
{
  df<- .avisApiFichaUsuario(contributor_id)

  names (df)<- c("SpeciesId", "Observations", "Number", "UTM.10x10", "Birdwatchers")

  return(df)
}

.avisUserNameList <- function()
{
  cached_list = .avisCacheGet('ravis_username_id_list')

  if(is.null(cached_list)){
    # ravis_username_id_list (cached object) is a by-product of avisContributorsSummary process
    avisContributorsSummary()
  }

  .avisCacheGet('ravis_username_id_list')
}

.avisContributorsSummaryInit <- function()
{
  rbs<- .avisApiUsuarios()

  ravis_username_id_list<-unlist(rbs[,2])
  names(ravis_username_id_list)<-unlist(rbs[,1])

  # remove username
  rbs <- rbs[,-(which (colnames(rbs)=="User"))]

  # cache username - id list
  .avisCacheSet('ravis_username_id_list', ravis_username_id_list)
    
  return(rbs)
}
