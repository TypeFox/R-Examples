#' Search for data on documents
#'
#' @export
#'
#' @param document_id A unique id for each document.
#' @param document_type Document types are taken from the House document type code. See Details.
#' @param chamber House or Senate
#' @param committee_id Acronym a committee is associated with the document.
#' @param committee_names Full names of the committees associated with the document.
#' @param congress Session of Congress.
#' @param house_event_id Unique ID for each hearing, assigned by the House.
#' @param hearing_type_code This describes if the meeting is a "markup", "meeting" or "hearing".
#' @param hearing_title Title of the hearing associated with the document.
#' @param published_at Date and time of publication.
#' @param bill_id Bill ID associated with the document.
#' @param description Description of the hearing.
#' @param version_code The short-code for what stage the version of the bill. See GPO for
#' explanations of the version code.
#' @param bioguide_id Unique identifier for a member of Congress if they are associated
#' with the document.
#' @param occurs_at Date and time of a hearing associated with the document.
#' @param urls The original link to the document. The permalink is a link to a copy
#' of the document hosted by the Sunlight Foundation.
#' @param text Extracted text from the document.
#' @param text_preview A preview of the text.
#' @param witness Information about a witness associated with a document.
#' @details
#' \itemize{
#'  \item CV - Committee vote
#'  \item WS - Witness statement
#'  \item WT - Witness truth statement
#'  \item WB - Witness biography
#'  \item CR - Committee report
#'  \item BR - Bill
#'  \item FA - Floor amendment
#'  \item CA - Committee amendment
#'  \item HT - Transcript
#'  \item WD - Witness document
#' }
#' \emph{other} is used for all other documents.
#'
#' @template cg
#' @examples \dontrun{
#' cg_documents()
#' cg_documents(per_page=4)
#' cg_documents(per_page=4)
#' cg_documents(fields=c('document_id','type'))
#'
#' # most parameters are vectorized, pass in more than one value
#' cg_documents(chamber = c("house", "senate"))
#' }

cg_documents <- function(document_id=NULL, document_type=NULL, chamber=NULL, committee_id=NULL,
  committee_names=NULL, congress=NULL, house_event_id=NULL, hearing_type_code=NULL,
  hearing_title=NULL, published_at=NULL, bill_id=NULL, description=NULL, version_code=NULL,
  bioguide_id=NULL, occurs_at=NULL, urls=NULL, text=NULL, text_preview=NULL, witness=NULL,
  fields=NULL, page=1, per_page=20, order=NULL,
  key = NULL, as = 'table', ...) {

  key <- check_key(key)
  fields <- paste0(fields, collapse = ",")
  args <- sc(list(apikey=key, per_page=per_page, page=page, order=order, fields=fields,
                          document_id=document_id, document_type=document_type, chamber=chamber,
                          committee_id=committee_id, committee_names=committee_names,
                          congress=congress, house_event_id=house_event_id,
                          hearing_type_code=hearing_type_code, hearing_title=hearing_title,
                          published_at=published_at, bill_id=bill_id, description=description,
                          version_code=version_code, bioguide_id=bioguide_id, occurs_at=occurs_at,
                          urls=urls, text=text, text_preview=text_preview, witness=witness))
  give_cg(as, cgurl(), "/congressional_documents/search", args, ...)
}
