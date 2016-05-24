#' Check Italian fiscal codes
#' 
#' The function performs fiscal codes check (both regular and temporary
#' ones), computing the last (control) character basing on the others.  The
#' algorithm is intended to be quite "draconian", so you'll better make
#' some pre-cleaning (keep only alphabetic characters and uppercase) if
#' needed. See examples.
#' 
#' @param fc A character vector of fiscal codes.
#' @return The function return \code{TRUE} if the fiscal code is wrong,
#' \code{FALSE} otherwise.
#' @references Law source: D.M. (Ministry of Finance) n. 13813 - 23/12/76 -
#' "Sistemi di codificazione dei soggetti da iscrivere all'anagrafe
#' tributaria". Supp. ord. G.U. 345 29/12/1976.
#' @examples
#' 
#' fc <- c(NA, "qWeASd34D12h 221M   ", " 12312312312 ")
#' wrong_fc(fc)
#' fc <- gsub(" ","", toupper(fc))
#' wrong_fc(fc)
#' 
#' @export 
wrong_fc <- function(fc) 
{
  if( ! (is.character(fc) & is.null(dim(fc))) )
    stop("The input must be a character vector.")
    
  ## Matching patterns and dummy indexes
  reg_ptrn <- "[A-Z]{6}\\d{2}[A-Z]{1}\\d{2}[A-Z]{1}\\d{3}[A-Z]{1}" 
  tmp_ptrn <- "\\d{11}"
  reg_indx <- grepl(reg_ptrn, fc, perl = TRUE)
  tmp_indx <- grepl(tmp_ptrn, fc, perl = TRUE)

  ## Results: wrong until proven to be right (unless NA, handled below)
  fc_error <- rep(TRUE, length(fc))

  ## Check regular fc: to keep C side simple, regular and temporary
  ## codes are checked in two separate functions
  if (any(reg_indx))
    fc_error[reg_indx] <-
      .Call("reg_wrong_fc", fc[reg_indx], PACKAGE="ifctools")

  ## Check temporary fc (extended = FALSE causes more testing to be needed)
  if (any(tmp_indx))
    fc_error[tmp_indx] <-
      .Call("tmp_wrong_fc", fc[tmp_indx], PACKAGE="ifctools")

  ## managing NAs
  fc_error[is.na(fc)] <- NA

  ## Return
  fc_error
}
