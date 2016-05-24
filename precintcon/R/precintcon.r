#' @title Precipitation Intensity, Concentration and Anomaly Analysis
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @author Jonas Teixeira Nery \email{jonas@@ourinhos.unesp.br}
#' @description Contains functions to analyze the precipitation 
#' intensity, concentration and anomaly.
#' @name precintcon
#' @docType package
#' @keywords precipitation rainfall anomaly intensity concentration package
#' @references
#' 
#' Michael J. Hayes (2015). \emph{Comparison of Major Drought Indices}. 
#' National Drought Mitigation Center. \url{http://goo.gl/skHfs9}.
#' 
#' Gibbs, W. J. (1967). \emph{Rainfall deciles as drought indicators}.
#' Bureau of Meteorology Bulletin, n. 48. Commonwealth of Australia. 
#' 
#' Luis, M. D., Gonzolez-Hidalgo, J. C., Brunetti, M., Longares, L. A. (2011). 
#' \emph{Precipitation concentration changes in Spain 1946-2005}. 
#' Natural Hazards and Earth System Science, v. 11, n. 5, p. 1259--1265. 
#' 
#' Martin-Vide, J. (2004). \emph{Spatial distribution of a daily precipitation  
#' concentration index in peninsular Spain}. International Journal of Climatology, 
#' v. 24, n. 8, p. 959--971.
#' 
#' Guttman, N. B. (1999). \emph{Accepting the standardized precipitation index: a 
#' calculation algorithm}. Journal of the American Water Resources Association, v. 35,
#' n. 2, p. 311--322. Blackwell Publishing Ltd.
#' 
#' Barring, L., Hulme, M. (1991). \emph{Filters and approximate confidence intervals for 
#' interpreting rainfall anomaly indices}. Journal of climate, v. 4, n. 8, p. 837--847.
#' @importFrom ggplot2 ggplot aes_string geom_bar scale_x_date
#' scale_y_continuous theme element_text facet_grid geom_line geom_boxplot element_blank
#' ylim ggsave xlab scale_fill_discrete geom_histogram geom_density geom_ribbon scale_x_continuous
#' scale_color_hue scale_fill_hue geom_point
#' @importFrom scales percent_format date_format date_breaks
#' @importFrom graphics par
#' @importFrom stats aggregate pgamma pnorm quantile sd uniroot var
#' @importFrom utils read.table tail
#' @importFrom graphics par
NULL

#' Daily precipitation between 1976 and 2010
#' 
#' This data set gives the daily precipitation (in millimeters) between 1976 and 2010, as compiled by the National Water Agency (ANA), Brazil.
#' 
#' @format A data frame with 420 observations on the following 33 variables:
#' \describe{
#'    \item{\code{year}}{the year.}
#'    \item{\code{month}}{the month.}
#'    \item{\code{d1}}{the precipitation value in millimeters of the day #1 of the month.}
#'    \item{\code{d2}}{the precipitation value in millimeters of the day #2 of the month.}
#'    \item{\code{d3}}{the precipitation value in millimeters of the day #3 of the month.}
#'    \item{\code{d4}}{the precipitation value in millimeters of the day #4 of the month.}
#'    \item{\code{d5}}{the precipitation value in millimeters of the day #5 of the month.}
#'    \item{\code{d6}}{the precipitation value in millimeters of the day #6 of the month.}
#'    \item{\code{d7}}{the precipitation value in millimeters of the day #7 of the month.}
#'    \item{\code{d8}}{the precipitation value in millimeters of the day #8 of the month.}
#'    \item{\code{d9}}{the precipitation value in millimeters of the day #9 of the month.}
#'    \item{\code{d10}}{the precipitation value in millimeters of the day #10 of the month.}
#'    \item{\code{d11}}{the precipitation value in millimeters of the day #11 of the month.}
#'    \item{\code{d12}}{the precipitation value in millimeters of the day #12 of the month.}
#'    \item{\code{d13}}{the precipitation value in millimeters of the day #13 of the month.}
#'    \item{\code{d14}}{the precipitation value in millimeters of the day #14 of the month.}
#'    \item{\code{d15}}{the precipitation value in millimeters of the day #15 of the month.}
#'    \item{\code{d16}}{the precipitation value in millimeters of the day #16 of the month.}
#'    \item{\code{d17}}{the precipitation value in millimeters of the day #17 of the month.}
#'    \item{\code{d18}}{the precipitation value in millimeters of the day #18 of the month.}
#'    \item{\code{d19}}{the precipitation value in millimeters of the day #19 of the month.}
#'    \item{\code{d20}}{the precipitation value in millimeters of the day #20 of the month.}
#'    \item{\code{d21}}{the precipitation value in millimeters of the day #21 of the month.}
#'    \item{\code{d22}}{the precipitation value in millimeters of the day #22 of the month.}
#'    \item{\code{d23}}{the precipitation value in millimeters of the day #23 of the month.}
#'    \item{\code{d24}}{the precipitation value in millimeters of the day #24 of the month.}
#'    \item{\code{d25}}{the precipitation value in millimeters of the day #25 of the month.}
#'    \item{\code{d26}}{the precipitation value in millimeters of the day #26 of the month.}
#'    \item{\code{d27}}{the precipitation value in millimeters of the day #27 of the month.}
#'    \item{\code{d28}}{the precipitation value in millimeters of the day #28 of the month.}
#'    \item{\code{d29}}{the precipitation value in millimeters of the day #29 of the month.}
#'    \item{\code{d30}}{the precipitation value in millimeters of the day #30 of the month.}
#'    \item{\code{d31}}{the precipitation value in millimeters of the day #31 of the month.}
#'  }
#' @name daily
#' @aliases daily
#' @usage data(daily)
#' @docType data
#' @source National Water Agency (ANA) \url{http://www2.ana.gov.br/Paginas/default.aspx}
#' @examples 
#' data(daily)
#' ## maybe pci(daily) ; rai(daily) ...
#' @keywords datasets
NULL

#' Monthly precipitation between 1950 and 1992.
#' 
#' This dataset gives the monthly precipitation (in millimeters) between 1950 and 1992, as compiled by the National Water Agency (ANA), Brazil.
#' 
#' @name monthly
#' @aliases monthly
#' @usage data(monthly)
#' @format 
#' A data frame with 516 observations on the following 3 variables:
#' \describe{
#'    \item{\code{year}}{the year.}
#'    \item{\code{month}}{the month.}
#' 	\item{\code{precipitation}}{the precipitation amount in millimeters.}
#' }
#' @docType data
#' @source National Water Agency (ANA) \url{http://www2.ana.gov.br/Paginas/default.aspx}
#' @examples
#' data(monthly)
#' ## maybe pci(monthly) ; rai(monthly) ...
#' @keywords datasets
NULL