#' A collection of wrapper functions aimed at for aiding the authoring of reproducible
#' reports.
#'
#' \pkg{qwraps2} is a collection of helpful functions when
#' working on a varied collection of different analysis reports.  There are two
#' types of functions, helpful data summary functions, formatting results from
#' regression models, and \pkg{ggplot2} wrappers.
#'
#' Several wrappers for \pkg{ggplot2} style graphics, such as ROC, AUC,
#' Bland-Altman, and KM
#' plots are provided.  Named as \code{\link{qroc}}, \code{\link{qacf}}, 
#' \code{\link{qblandaltman}} and
#' \code{\link{qkmplot}} to pay homage to \code{qplot} form \pkg{ggplot2} and
#' the standard names for such plots.
#'
#' Other functions are used to quickly generate meaningful character strings for
#' outputting results in .Rnw, .Rmd, or other similar functions.
#'
#' @section Options: 
#' There are several options which can be set via \code{options} and will be
#' used via \code{getOption}.  The following lists, in alphabetical order the
#' different options which are available and what they control.
#'
#' \itemize{
#'   \item \code{getOptions("qwraps2_alpha", 0.05)} significance level, used for
#'   generating \code{(1 - getOptions("qwraps2_alpha", 0.05)) * 100}\% confidence
#'   intervals, and determining significance for p-value <
#'   \code{getOptions("qwraps2_alpha", 0.05)}.
#'
#'   \item \code{getOptions("qwraps2_frmt_digits", 2)}  Number of digits to the
#'   right of the decimal point for any value other than p-values.
#'
#'   \item \code{getOptions("qwraps2_frmtp_case", "upper")} set to either
#'   'upper' or 'lower' for the case of the 'P' for reporting p-values.
#'
#'   \item \code{getOptions("qwraps2_frmtp_digits", 4)}  Number of digits to the
#'     right of the decimal point to report p-values too.  If 
#'     \code{log10(p-value) < getOptions("qwraps2_frmtp_digits", 4)} then the
#'     output will be "P < 0.01", to however many digits are correct.  Other
#'     options control other parts of the output p-value format.
#'
#'   \item \code{getOptions("qwraps2_frmtp_leading0", TRUE)} to display or not
#'   to display the leading zero in p-values, i.e., if TRUE p-values are
#'   reported as 0.02 versus when FALSE p-values are reported as .02.
#'
#'   \item \code{getOptions("qwraps2_journal", "default")} if a journal has
#'   specific formating for p-values or other statistics, this option will
#'   control the output.  Many other options are ignored if this is any other
#'   than default.  Check the github wiki, or this file,  for current lists of
#'   implemented journal style methods.
#'
#'   \item \code{getOptions("qwraps2_markup", latex)} value set to 'latex' or
#'   to 'markdowm'.  Output is formatted to meet requirements of either markup
#'   language.
#'
#'   \item \code{getOptions("qwraps2_style", "default")} By setting this option
#'     to a specific journal, p-values and other output, will be formatted to
#'     meet journal requirements.
#' }
#'
#' @section Journals with predefined formatting: 
#' 
#' \itemize{
#'   \item Obstetrics \& Gynecology 
#'   \itemize{
#'      \item \url{http://www.editorialmanager.com/ong/default.aspx} 
#'      \item \code{options(qwraps2_journal = "obstetrics_gynecology")}
#'      \item P-value formating as of April 2015:
#' 
#'       Express P values to no more than three decimal places.
#'
#'       Based on observations of published work, leading 0 will be omitted.
#'    }
#'
#'   \item Pediatric Dentistry:  
#'   \itemize{
#'     \item \url{http://www.aapd.org/publications/} 
#'     \item \code{options(qwraps2_journal = "pediatric_dentistry")}
#'     \item P-value formating  as of March 2015.  
#'
#'     If P > .01, the actual value for
#'     P should be expressed to 2 digits.  Non-significant values should not be
#'     expressed as "NS" whether or note P is significant, unless rounding a
#'     significant P-value expressed to 3 digits would make it non significant (ie
#'     P=.049, not P=.05).  If P<.01, it should be express to 3 digits (eg, P=.003,
#'     not P<.05).  Actual P-values should be expressed unless P<.001, in which case
#'     they should be so designated.
#'   }
#'
#' }
#' @docType package
#' @name qwraps2
NULL
