#' Search BHL across many API methods.
#'
#' @export
#' @inheritParams bhl_getcollections
#' @param method The API method to use.
#' @param pageid The identifier of an individual page in a scanned book.
#' @param ocr return ocr text of the page (logical). Default: FALSE
#' @param names return the names that appear on the page (logical). Default: FALSE
#'
#' @examples \dontrun{
#' bhl_bioherlib(method='GetPageMetadata', pageid=1328690, ocr=TRUE, names=TRUE)
#' bhl_bioherlib(method='GetPageMetadata', pageid=1328690, ocr=TRUE, names=TRUE, as="xml")
#' bhl_bioherlib(method='GetPageMetadata', pageid=1328690, ocr=TRUE, names=TRUE, as="list")
#' }

bhl_bioherlib <- function(method = 'GetPageMetadata', pageid = NULL, ocr = FALSE,
  names = FALSE, as='table', key = NULL, ...) {

  method <- match.arg(method,
    choices = c('GetPageMetadata', 'GetPageOcrText', 'GetPageNames',
     'GetItemMetadata', 'GetItemByIdentifier', 'GetItemPages', 'GetUnpublishedItems',
     'GetTitleMetadata', 'GetTitleItems', 'GetTitleByIdentifier', 'TitleSearchSimple',
     'GetTitleBibTex', 'GetTitleEndNote', 'GetUnpublishedTitles', 'SubjectSearch',
     'GetSubjectTitles', 'AuthorSearch', 'GetAuthorTitles', 'NameCount', 'NameList',
     'NameGetDetail', 'NameSearch', 'GetCollections', 'GetLanguages'))
  args <- bhlc(list(apikey = check_key(key), op = method, pageid = pageid, format = as_f(as),
                       ocr = ocr, names = names))
  bhl_GET(as, args, ...)
}

