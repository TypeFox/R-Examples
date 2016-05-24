#' Return details about the reports you have permission to view.
#' 
#' @inheritParams form_info
#' @inheritParams user_info
#' 
#' @param reportName - the name of your report as displayed in the csv export 
#' URL (which is in lowercase with hyphens replacing spaces of your report name). Do not use 
#' the hash value for this argument. The report should be also public.
#' 
#' @return Name - This is the friendly name you chose when creating this form.
#' @return IsPublic - Indicates whether or not the report is accessible through 
#' the Url by the general public. This value is binary (1 = true, 0 = false).
#' @return Url - This is the URL for your form. Beware using the URL for API or 
#' linking purposes because it changes with the report title.
#' @return Description - Your description of the report.
#' @return DateCreated - A timestamp of when the report was created.
#' @return DateUpdated - A timestamp of when the report was lasted edited in the Wufoo Report Builder.
#' @return Hash - An unchanging hashed value unique to this report on this user's account.
#' 
#' @examples
#' reports_info(showRequestURL = TRUE)
#' reports_info(reportName = "untitled-report")
#' 
#' @export
reports_info <- function(wufoo_name = auth_name(NULL), domain = "wufoo.com", 
                         reportName = NULL, showRequestURL = FALSE, debugConnection = 0L) {
  
  reports_url <- paste0("https://", wufoo_name, ".", domain, "/api/v3/reports.json")
  
  query <- list(reportName = reportName)
  
  executedReportsGetRst <- doRequest(reports_url, query, showURL = showRequestURL, debugConnection = debugConnection)
  
  return(executedReportsGetRst$Reports)
}












