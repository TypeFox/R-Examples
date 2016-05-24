#' Search BOLD for sequences.
#'
#' Get sequences for a taxonomic name, id, bin, container, institution, researcher, geographic
#' place, or gene.
#'
#' @importFrom stringr str_replace_all str_replace str_split
#' @export
#' @template args
#' @template otherargs
#' @references \url{http://www.boldsystems.org/index.php/resources/api#sequenceParameters}
#'
#' @param marker (character) Returns all records containing matching marker codes.
#'
#' @return A list with each element of length 4 with slots for id, name, gene, and sequence.
#' @examples \dontrun{
#' bold_seq(taxon='Coelioxys')
#' bold_seq(taxon='Aglae')
#' bold_seq(taxon=c('Coelioxys','Osmia'))
#' bold_seq(ids='ACRJP618-11')
#' bold_seq(ids=c('ACRJP618-11','ACRJP619-11'))
#' bold_seq(bin='BOLD:AAA5125')
#' bold_seq(container='ACRJP')
#' bold_seq(researchers='Thibaud Decaens')
#' bold_seq(geo='Ireland')
#' bold_seq(geo=c('Ireland','Denmark'))
#'
#' # Return the httr response object for detailed Curl call response details
#' res <- bold_seq(taxon='Coelioxys', response=TRUE)
#' res$url
#' res$status_code
#' res$headers
#'
#' ## curl debugging
#' ### You can do many things, including get verbose output on the curl call, and set a timeout
#' library("httr")
#' bold_seq(taxon='Coelioxys', config=verbose())[1:2]
#' # bold_seqspec(taxon='Coelioxys', config=timeout(0.1))
#' }

bold_seq <- function(taxon = NULL, ids = NULL, bin = NULL, container = NULL, institutions = NULL,
  researchers = NULL, geo = NULL, marker = NULL, response=FALSE, ...) {
  args <- bc(list(taxon=pipeornull(taxon), geo=pipeornull(geo), ids=pipeornull(ids),
      bin=pipeornull(bin), container=pipeornull(container), institutions=pipeornull(institutions),
      researchers=pipeornull(researchers), marker=pipeornull(marker)))
  check_args_given_nonempty(args, c('taxon','ids','bin','container','institutions','researchers',
                                    'geo','marker'))
  out <- b_GET(paste0(bbase(), 'API_Public/sequence'), args, ...)
  if (response) { 
    out 
  } else {
    tt <- rawToChar(content(out, encoding = "UTF-8"))
    res <- strsplit(tt, ">")[[1]][-1]
    lapply(res, split_fasta)
  }
}
