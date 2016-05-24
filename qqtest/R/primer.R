#' Automobile primer paint thickness quality control measurements.
#'
#' Contains process control measurements of thickness of primer applied to
#' automotive body parts in an auto factory.
#' Twice daily, a set of 10 consecutive parts were selected and the thickness in mils (thousandths of an inch)
#' were measured.  For each set of 10 parts, the average (xbar) and the sample standard deviation (s) were also
#' calculated and recorded.  These summaries would be plotted in xbar or s control charts with suitably determined upper and
#' lower control limits.
#' Alternatively, for checking outliers a qqplot (via qqtest) could be used for either xbar or s.
#' 
#' \code{with(primer,qqtest(xbar, main="Averages"))} will effect this plot for xbar.
#' \code{with(primer,qqtest(s,dist="kay", df=9, main ="Standard deviations"))} will effect this plot for s.
#'
#'
#' @format A data frame with 20 rows and 14 variates:
#' \describe{
#'   \item{day}{Day on which the parts were taken and measured.}
#'   \item{batch}{Either the first or second set of 10 consecutive parts taken.}
#'   \item{sel1}{Thickness of primer in mils on the first part sampled in the specified batch of that day.}
#'   \item{sel2}{Thickness of primer in mils on the second part sampled in the specified batch of that day.}
#'   \item{sel3}{Thickness of primer in mils on the third part sampled in the specified batch of that day.}
#'   \item{sel4}{Thickness of primer in mils on the fourth part sampled in the specified batch of that day.}
#'   \item{sel5}{Thickness of primer in mils on the fifth part sampled in the specified batch of that day.}
#'   \item{sel6}{Thickness of primer in mils on the sixth part sampled in the specified batch of that day.}
#'   \item{sel7}{Thickness of primer in mils on the seventh part sampled in the specified batch of that day.}
#'   \item{sel8}{Thickness of primer in mils on the eighth part sampled in the specified batch of that day.}
#'   \item{sel9}{Thickness of primer in mils on the ninth part sampled in the specified batch of that day.}
#'   \item{sel10}{Thickness of primer in mils on the tenth part sampled in the specified batch of that day.}
#'   \item{xbar}{Arithmetic average of the measurements of primer thickness of the 10 parts selected in the specified batch of that day.}
#'   \item{s}{Sample standard deviation of the measurements of primer thickness of the 10 parts selected in the specified batch of that day.}
#' }

#' @source 
#' "Statistical Process Control - SPC", 
#' Automotive Industry Action Group(AIAG), Southfield MI,  (1995), page 64.
#' 
"primer"

primer <- data.frame( day = rep(1:10, each=2),
				   batch = rep(c(1,2), 10),
				   
				   sel1  = c(1.30, 1.01, 1.22, 1.08, 0.98, 1.12, 0.92, 1.04, 1.08, 1.20, 1.25, 1.24, 1.13, 1.08, 
				   1.08, 1.14, 1.06, 1.14, 1.07, 1.13),
				   sel2  = c(1.10, 1.10, 1.05, 1.12, 1.30, 1.30, 1.10, 1.14, 0.92, 1.13, 0.91, 1.34, 1.16, 1.31, 
				   1.26, 1.02, 1.12, 1.22, 1.05, 0.90),
				   sel3  = c(1.20, 1.15, 0.93, 1.11, 1.31, 1.01, 1.13, 1.18, 1.14, 1.19, 0.96, 1.40, 1.12, 1.12, 
				   1.13, 1.14, 0.98, 1.18, 0.97, 1.12), 
				   sel4  = c(1.25, 0.97, 1.08, 1.28, 1.12, 1.20, 1.02, 1.12, 1.20, 1.16, 1.04, 1.26, 1.22, 1.18, 
				   0.94, 0.94, 1.12, 1.27, 1.05, 1.04), 
				   sel5  = c(1.05, 1.25, 1.15, 1.00, 1.08, 1.11, 0.93, 1.00, 1.02, 1.03, 0.93, 1.13, 1.12, 1.15, 
				   1.30, 1.30, 1.20, 1.17, 1.16, 1.40), 
				   sel6  = c(0.95, 1.12, 1.27, 0.95, 1.10, 0.93, 1.17, 1.02, 1.04, 1.25, 1.08, 1.15, 1.07, 1.17, 
				   1.15, 1.08, 1.02, 1.26, 1.02, 1.12), 
				   sel7  = c(1.10, 1.10, 0.95, 1.15, 1.15, 1.02, 1.24, 1.05, 0.94, 1.20, 1.29, 1.08, 1.04, 0.98, 
				   1.07, 0.94, 1.19, 1.15, 1.02, 1.15), 
				   sel8  = c(1.16, 0.90, 1.11, 1.14, 1.35, 1.25, 0.98, 1.34, 1.05, 1.24, 1.42, 1.02, 1.28, 1.05, 
				   1.02, 1.12, 1.03, 1.07, 1.14, 1.01),
				   sel9  = c(1.37, 1.04, 1.12, 1.28, 1.12, 1.05, 1.34, 1.12, 1.12, 1.10, 1.10, 1.05, 1.12, 1.00, 
				   1.22, 1.15, 1.02, 1.02, 1.07, 1.30), 
				   sel10 = c(0.98, 1.08, 1.10, 1.31, 1.26, 1.10, 1.12, 1.05, 1.06, 1.03, 1.00, 1.18, 1.10, 1.26, 
				   1.18, 1.36, 1.09, 1.36, 1.00, 1.14),
					 
				   xbar  = c(1.15, 1.07, 1.10, 1.14, 1.18, 1.11, 1.10, 1.11, 1.06, 1.15, 1.10, 1.19, 1.14, 1.13, 
				   1.14, 1.12, 1.08, 1.18, 1.06, 1.13),
				   s     = c(.136, .098, .106, .120, .121, .115, .136, .101, .086, .079, .170, .125, .070, .107, 
				   .111, .137, .074, .099, .059, .141)
							)
							
	


