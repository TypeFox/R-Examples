#' Compliance Report
#'
#' Generate compliance report by time across treatment groups.
#'
#' @param comply numeric or character vector. Indicator variable for compliance.
#' Should be 1/0 or yes/no.
#' @param treat factor vector. Treatment group for each record.
#' @param time numeric vector. Time for each record.
#' @param times numeric vector. Subset of times to use.
#' @export
#' @examples
#' \dontrun{
#'   complianceReport(rbinom(200, 1, 0.8), as.factor(sample(c('A','B'), 200, replace=TRUE)), sample(10, 200, replace=TRUE))
#' }

complianceReport <- function(comply, treat, time, times=NULL) {
  if(!is.numeric(comply))
    comply <-  1*(comply %in% c('Y','Yes','yes','YES'))
  if(length(times)) {
    s <- time %in% times
    comply <- comply[s]
    treat <- treat[s]
    time <- time[s]
  }
  Compliance <- comply
  latex(summary(Compliance ~ time+stratify(treat)), 
        file=file.path(TexDirName(), 'compliance.tex'), where='hbp!', ncaption=FALSE, ctable=TRUE)
  latex(summary(Compliance ~ time),
        file=file.path(TexDirName(), 'Ocompliance.tex'), where='hbp!', ncaption=FALSE, ctable=TRUE)
  invisible()
}
