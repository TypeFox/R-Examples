#' Get BOLD specimen + sequence data.
#'
#' @export
#' @template args
#' @template otherargs
#' @references \url{http://www.boldsystems.org/index.php/resources/api#combined}
#'
#' @param marker (character) Returns all records containing matching marker codes.
#' @param format (character) One of xml or tsv (default). tsv format gives back a data.frame
#' object. xml gives back parsed xml as a
#' @param sepfasta (logical) If TRUE, the fasta data is separated into a list with names matching
#' the processid's from the data frame
#'
#' @return Either a data.frame, parsed xml, a httr response object, or a list with length two
#' (a data.frame w/o nucleotide data, and a list with nucleotide data)
#'
#' @examples \dontrun{
#' bold_seqspec(taxon='Osmia')
#' bold_seqspec(taxon='Osmia', format='xml')
#' bold_seqspec(taxon='Osmia', response=TRUE)
#' res <- bold_seqspec(taxon='Osmia', sepfasta=TRUE)
#' res$fasta[1:2]
#' res$fasta['GBAH0293-06']
#' 
#' # records that match a marker name
#' res <- bold_seqspec(taxon="Melanogrammus aeglefinus", marker="COI-5P")
#' 
#' # records that match a geographic locality
#' res <- bold_seqspec(taxon="Melanogrammus aeglefinus", geo="Canada")
#'
#' ## curl debugging
#' ### You can do many things, including get verbose output on the curl call, and set a timeout
#' library("httr")
#' head(bold_seqspec(taxon='Osmia', config=verbose()))
#' ## timeout
#' # head(bold_seqspec(taxon='Osmia', config=timeout(1)))
#' ## progress
#' # x <- bold_seqspec(taxon='Osmia', config=progress())
#' }

bold_seqspec <- function(taxon = NULL, ids = NULL, bin = NULL, container = NULL,
  institutions = NULL, researchers = NULL, geo = NULL, marker = NULL, response=FALSE,
  format = 'tsv', sepfasta=FALSE, ...) {
  
  format <- match.arg(format, choices = c('xml','tsv'))
  args <- bc(list(taxon = pipeornull(taxon), geo = pipeornull(geo), 
                  ids = pipeornull(ids), bin = pipeornull(bin), 
                  container = pipeornull(container), 
                  institutions = pipeornull(institutions), 
                  researchers = pipeornull(researchers), 
                  marker = pipeornull(marker), combined_download = format))
  check_args_given_nonempty(args, c('taxon', 'ids', 'bin', 'container',
                                    'institutions', 'researchers',
                                    'geo', 'marker'))
  out <- b_GET(paste0(bbase(), 'API_Public/combined'), args, ...)
  if (response) { 
    out
  } else {
    tt <- rawToChar(content(out, encoding = "UTF-8"))
    if (tt == "") return(NA)
    temp <- switch(format,
           xml = xml2::read_xml(tt),
           tsv = read.delim(text = tt, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    )
    if (!sepfasta) { 
      temp 
    } else {
      if (format == "tsv") {
        fasta <- as.list(temp$nucleotides)
        names(fasta) <- temp$processid
        df <- temp[ , !names(temp) %in% "nucleotides" ]
        list(data = df, fasta = fasta)
      } else { 
        temp 
      }
    }
  }
}
