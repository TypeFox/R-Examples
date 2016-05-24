#' @title Return details about the forms you have permission to access.
#'
#' @description This API can be used to create a list of all forms belonging to a user and 
#' dynamically generate a form embed snippet to use in your application.
#'
#' @param formIdentifier - this will give you information about just one form. The call without 
#' the "formIdentifier" will return all forms.
#' @param includeTodayCount - Will give you today's entry count for the form.
#' 
#' @inheritParams user_info
#'
#' @return \url{http://help.wufoo.com/articles/en_US/SurveyMonkeyArticleType/The-Forms-API}
#' 
#' @examples 
#' form_info()
#' 
#' @export
form_info <- function(wufoo_name = auth_name(NULL), formIdentifier = NULL, includeTodayCount = "false", 
                      showRequestURL = FALSE, debugConnection = 0L, domain = "wufoo.com") {
  
  form_url <- paste0("https://", wufoo_name, ".", domain, "/api/v3/forms.json")
  
  query = list(formIdentifier = formIdentifier, includeTodayCount = includeTodayCount)
  
  executedFormGetRst <- doRequest(form_url, query, showURL = showRequestURL, debugConnection = debugConnection)
  df_forms <- executedFormGetRst$Forms
  
  return(df_forms)
}


#' Return responses of your form
#' 
#' @seealso \url{http://help.wufoo.com/articles/en_US/SurveyMonkeyArticleType/The-Entries-GET-API}
#' 
#' @inheritParams form_info
#' @inheritParams user_info
#' 
#' @param systemFields - return system fields. Default: true
#' @param formIdentifier - must be replaced with your form's URL or hash.
#' @param columnNames - a MUST: How should be column names be called. Either "Field1", "Field2" 
#' etc. or "First Name", "Last Name" (tries to make best guess). Default to the second option. 
#' @param sortID - sort on a single ID, as retrieved from the \code{\link{fields_info}}.
#' @param sortDirection - choose to sort your entries ASC (lowest to highest) or DESC (highest to lowest).
#' @param pageStart - the page number you'd like to start from.  Defaults to 0.
#' @param pageSize - the number of entries returned in your page. Defaults to 25; Max = 100.
#' 
#' @description If you have 5 submissions to your form, you'll have 5 elements (rows) in the return.
#' 
#' @return EntryId - This value is the unique identifier for your entry.
#' @return DateCreated - The date that this entry was submitted.
#' @return Created By - The person who created the entry. If submitted through a form, the value 
#' here will be public. If the submission originated in the Entry Manager this value will be the 
#' user name of the submitting user.
#' @return DateUpdated - The date that this entry was edited through the Entry Manager. If the 
#' submission has never been updated, this value will be blank.
#' @return UpdatedBy - The user name of the person who updated the entry in the Entry Manager will 
#' appear in this element.
#' 
#' @examples
#' form_entries(formIdentifier = "z5kqx7h1gtvg4g")
#' form_entries(formIdentifier = "z5kqx7h1gtvg4g", systemFields = "false", showRequestURL = TRUE)
#' 
#' @import dplyr
#' 
#' @export
form_entries <- function(wufoo_name = auth_name(NULL), formIdentifier = NULL, systemFields = "true", 
                         sortID = NULL, sortDirection = NULL, columnNames = FALSE, showRequestURL = FALSE,
                         debugConnection = 0L, domain = "wufoo.com", pageStart = 0, pageSize = 25) {
  
  entries_url <- paste0("https://", wufoo_name, ".", domain, "/api/v3/forms/", formIdentifier ,"/entries.json")
  
  query = list(systemFields = systemFields, sort = sortID, sortDirection = sortDirection, pageStart = pageStart, pageSize = pageSize)
  
  executedEntriesGetRst <- doRequest(entries_url, query, showURL = showRequestURL, debugConnection = debugConnection)
  
  df_entries <- executedEntriesGetRst$Entries
  
  # Make column names more understandable; use names of fields and apply them to the data frame
  if(identical(columnNames, FALSE)) {
    df_entries2 <- data.frame(t(df_entries))
    df_entries2$colNames <- rownames(df_entries2)
    
    fjoined <- fields_info(formIdentifier = formIdentifier) 
    
    # Merge two datasets and later replace names of `df_entries` with those in `df_mergedColNames`
    df_mergedColNames <- left_join(df_entries2, fjoined, by = c("colNames" = "ID"))
    
    colnames(df_entries) <- ifelse(!is.na(df_mergedColNames$Title), df_mergedColNames$Title, 
                                   df_mergedColNames$colNames)
  }
  
  df_entries[df_entries == ""] <- NA 
  
  return(df_entries)
}

#' Return number of responses to your form
#' 
#' @seealso \url{http://help.wufoo.com/articles/en_US/SurveyMonkeyArticleType/The-Entries-GET-API#entrycount}
#' 
#' @inheritParams form_info
#' @inheritParams user_info
#' @inheritParams form_entries
#'   
#' @return EntryCount - number of entries
#' 
#' @examples
#' form_entriesCount(formIdentifier = "z5kqx7h1gtvg4g", showRequestURL = TRUE)
#' 
#' @export
form_entriesCount <- function(wufoo_name = auth_name(NULL), formIdentifier = NULL, showRequestURL = FALSE, 
                              debugConnection = 0L, domain = "wufoo.com") {
  
  entriesCount_url <- paste0("https://", wufoo_name, ".", domain, "/api/v3/forms/", formIdentifier, "/entries/count.json")
  
  executedEntriesCountGetRst <- doRequest(entriesCount_url, showURL = showRequestURL, debugConnection = debugConnection)
  
  return(executedEntriesCountGetRst$EntryCount)
}


#' Return responses of your form, from CSV format
#' 
#' @description This function imports your report csv file as a data frame from the report csv export url (example url below). 
#' The report must be publica and not protected. To verify your report csv export url, 
#' browse to your report, select "Export" then hover over the "Commas (.csv)" option. 
#' Please note that the name of your report will be in lowercase with spaces replaced with hyphens. 
#' For example, the report titled "My Example Report" will be "my-example-report" in the URL as shown below.
#' E.g. \code{https://YourName.wufoo.com/export/reports/manager/NameOfYourReport.csv} 
#' 
#' @inheritParams reports_info
#' @inheritParams user_info
#' @inheritParams form_entries
#'  
#' @examples
#' \dontrun{
#' options(Wufoo_Name = "johnmalc", Wufoo_API = "F1QH-Q64B-BSBI-JASJ")
#' df_csv <- form_entriesFromCSV(reportName = "untitled-report", showRequestURL = F)
#' View(df_csv)
#' }
#' 
#' @importFrom utils read.csv
#'  
#' @export
form_entriesFromCSV <- function(wufoo_name = auth_name(NULL), reportName = NULL, showRequestURL = FALSE, 
                                debugConnection = 0L, domain = "wufoo.com") {
  
  entriesFromCSV_url <- paste0("https://", wufoo_name, ".", domain, "/export/report/manager/", reportName, ".csv")
  
  executedEntriesFromCSVGetRst <- doRequest(entriesFromCSV_url, showURL = showRequestURL, debugConnection = debugConnection)
  
  df_csv <- read.csv(text = executedEntriesFromCSVGetRst, stringsAsFactors = F, header = T, na.strings = c("NA", ""))
  
  return(df_csv)
}





