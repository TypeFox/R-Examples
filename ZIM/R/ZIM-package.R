#' @docType package
#' @name ZIM-package
#' @aliases ZIM
#' @title Zero-Inflated Models for Count Time Series with Excess Zeros
#' @description Fits observation-driven and parameter-driven models for count time series with excess zeros. 
#' @details The package \code{ZIM} contains functions to fit statistical models for count time series 
#' with excess zeros (Yang et al., 2013, 2014+). The main function for fitting observation-driven models 
#' is \code{\link{zim}}, and the main function for fitting parameter-driven models is \code{\link{dzim}}. 
#' @note The observation-driven models for zero-inflated count time series can also be fit using the function 
#' \code{\link[pscl]{zeroinfl}} from the \code{pscl} package (Zeileis et al., 2008). 
#' Fitting parameter-driven models is based on sequential Monte Carlo (SMC) methods, which are
#' computer intensive and could take several hours to estimate the model parameters.
#' @author Ming Yang, Gideon K. D. Zamba, and Joeseph E. Cavanaugh
#' @author Maintainer: Ming Yang <mingyang AT hsph.harvard.edu> 
#' @references
#' Yang, M., Cavanaugh, J. E., and Zamba, G. K. D. (2014+). State-space models for count time series
#' with excess zeros. \emph{Statistical Modelling}, To Appear. \cr \cr
#' Yang, M., Zamba, G. K. D., and Cavanaugh, J. E. (2013). Markov regression models for count time series
#' with excess zeros: A partial likelihood approach. \emph{Statistical Methodology}, \bold{14}:26-38. \cr \cr
#' Zeileis, A., Kleiber, C., and Jackman, S. (2008). Regression models for count data in \code{R}.  
#' \emph{Journal of Statistical Software}, \bold{27}(8). 
#' @keywords package
#' @import MASS
NULL


#' @docType data
#' @name injury
#' @title Example: Injury Series from Occupational Health
#' @description Monthly number of injuries in hospitals from July 1988 to October 1995.
#' @source Numbers from Figure 1 of Yau et al. (2004).
#' @references
#' Yau, K. K. W., Lee, A. H. and Carrivick, P. J. W. (2004). Modeling zero-inflated count series with 
#' application to occupational health. \emph{Computer Methods and Programs in Biomedicine}, \bold{74}, 47-52.
#' @examples
#' data(injury)
#' plot(injury, type = "o", pch = 20, xaxt = "n", yaxt = "n", ylab = "Injury Count")
#'   axis(side = 1, at = seq(1, 96, 8))
#'   axis(side = 2, at = 0:9)
#'   abline(v = 57, lty = 2)
#'   mtext("Pre-intervention",  line = 1, at = 25, cex = 1.5)
#'   mtext("Post-intervention", line = 1, at = 80, cex = 1.5)
#' @keywords datasets
NULL


#' @docType data
#' @name syph
#' @title Example: Syphilis Series
#' @description Weekly number of syphilis cases in the United States from 2007 to 2010.
#' @format A data frame with 209 observations on the following 69 variables.
#' \tabular{llll}{
#' \code{year} \tab Year                 \cr
#' \code{week} \tab Week                 \cr
#' \code{a1}   \tab \bold{United States} \cr
#' \code{a2}   \tab \bold{New England}   \cr
#' \code{a3}   \tab Connecticut          \cr
#' \code{a4}   \tab Maine                \cr
#' \code{a5}   \tab Massachusetts        \cr
#' \code{a6}   \tab New Hampshire        \cr
#' \code{a7}   \tab Rhode Island         \cr
#' \code{a8}   \tab Vermont              \cr
#' \code{a9}   \tab \bold{Mid. Atlantic} \cr
#' \code{a10}  \tab New Jersey           \cr
#' \code{a11}  \tab New York (Upstate)   \cr
#' \code{a12}  \tab New York City        \cr
#' \code{a13}  \tab Pennsylvania         \cr 
#' \code{a14}  \tab \bold{E.N. Central}  \cr
#' \code{a15}  \tab Illinois             \cr
#' \code{a16}  \tab Indiana              \cr
#' \code{a17}  \tab Michigan             \cr
#' \code{a18}  \tab Ohio                 \cr
#' \code{a19}  \tab Wisconsin            \cr
#' \code{a20}  \tab \bold{W.N. Central}  \cr
#' \code{a21}  \tab Iowa                 \cr
#' \code{a22}  \tab Kansas               \cr
#' \code{a23}  \tab Minnesota            \cr
#' \code{a24}  \tab Missouri             \cr
#' \code{a25}  \tab Nebraska             \cr
#' \code{a26}  \tab North Dakota         \cr
#' \code{a27}  \tab South Dakota         \cr
#' \code{a28}  \tab \bold{S. Atlantic}   \cr
#' \code{a29}  \tab Delaware             \cr
#' \code{a30}  \tab District of Columbia \cr
#' \code{a31}  \tab Florida              \cr
#' \code{a32}  \tab Georgia              \cr
#' \code{a33}  \tab Maryland             \cr
#' \code{a34}  \tab North Carolina       \cr
#' \code{a35}  \tab South Carolina       \cr
#' \code{a36}  \tab Virginia             \cr
#' \code{a37}  \tab West Virginia        \cr
#' \code{a38}  \tab \bold{E.S. Central}  \cr
#' \code{a39}  \tab Alabama              \cr
#' \code{a40}  \tab Kentucky             \cr
#' \code{a41}  \tab Mississippi          \cr
#' \code{a42}  \tab Tennessee            \cr
#' \code{a43}  \tab \bold{W.S. Central}  \cr
#' \code{a44}  \tab Arkansas             \cr
#' \code{a45}  \tab Louisana             \cr
#' \code{a46}  \tab Oklahoma             \cr
#' \code{a47}  \tab Texas                \cr
#' \code{a48}  \tab \bold{Moutain}       \cr
#' \code{a49}  \tab Arizona              \cr
#' \code{a50}  \tab Colorado             \cr
#' \code{a51}  \tab Idaho                \cr
#' \code{a52}  \tab Montana              \cr
#' \code{a53}  \tab Nevada               \cr
#' \code{a54}  \tab New Mexico           \cr
#' \code{a55}  \tab Utah                 \cr
#' \code{a56}  \tab Wyoming              \cr
#' \code{a57}  \tab \bold{Pacific}       \cr
#' \code{a58}  \tab Alaska               \cr
#' \code{a59}  \tab California           \cr
#' \code{a60}  \tab Hawaii               \cr
#' \code{a61}  \tab Oregon               \cr
#' \code{a62}  \tab Washington           \cr
#' \code{a63}  \tab American Samoa       \cr
#' \code{a64}  \tab C.N.M.I.             \cr
#' \code{a65}  \tab Guam                 \cr
#' \code{a66}  \tab Peurto Rico          \cr
#' \code{a67}  \tab U.S. Virgin Islands  \cr
#' }
#' @note C.N.M.I.: Commonwealth of Northern Mariana Islands.
#' @source CDC Morbidity and Mortality Weekly Report (\url{http://www.cdc.gov/MMWR/}).
#' @examples
#' data(syph)
#' plot(ts(syph$a33), main = "Maryland")
#' @keywords datasets
NULL
