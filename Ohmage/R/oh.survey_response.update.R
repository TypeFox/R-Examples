#' Update the privacy status of a survey response
#' @param campaign_urn campaign id
#' @param survey_key survey id
#' @param privacy_state Must be one of shared or private.
#' @param ... other arguments passed to oh.call
#' @export

oh.survey_response.update <- function(campaign_urn, survey_key, privacy_state="shared", ...){
	oh.call("/survey_response/update", campaign_urn=campaign_urn, survey_key=survey_key, privacy_state=privacy_state, ...)
	message("Survey ", survey_key, " has been shared.")
}

