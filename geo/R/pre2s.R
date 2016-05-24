#' Read prelude files
#' 
#' Read data files in prelude format, which has 2 line headers and column names
#' seperated from data with dashes.
#' 
#' 
#' @param skr Prelude file name
#' @param rownames Should first column be used as row names? Default FALSE
#' @param dots.in.text Should underscores in column names be replaced with "."?
#' Default TRUE
#' @return A data frame (\code{\link{data.frame}}) containing a representation
#' of the data in the file.
#' @note Call to \code{skipta.texta} could be replaced with a call to
#' \code{\link{chartr}} (as in ROracleUI sql).
#' @seealso Calls \code{\link{skipta.texta}}
#' @keywords ~kwd1
#' @export pre2s
pre2s <-
function(skr, rownames = F, dots.in.text = T)
{
        fields <- count.fields(skr, sep = "\t")
        nrec <- length(fields)
        if(nrec == 2)
                return(NULL)
        collab <- scan(file = skr, what = character(), sep = "\t",
                n = fields[1])
        outp <- read.table(skr, sep = "\t", skip = 2, as.is = T, row.names =
                NULL, na.strings = "")
        names(outp) <- collab
        if(rownames) {
                row.names(outp) <- outp[, 1]
                outp <- outp[, 2:ncol(outp)]
        }
        # change _ in names to .
        if(dots.in.text) names(outp) <- skipta.texta(names(outp))
        return(outp)
}

