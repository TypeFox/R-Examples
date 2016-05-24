#' Apply the maximum age model to a given De distribution
#'
#' Function to fit the maximum age model to De data. This is a wrapper function
#' that calls calc_MinDose() and applies a similiar approach as described in
#' Olley et al. (2006).
#'
#' \bold{Data transformation} \cr\cr To estimate the maximum dose population
#' and its standard error, the three parameter minimum age model of Galbraith
#' et al. (1999) is adapted. The measured De values are transformed as follows:
#' \cr\cr 1. convert De values to natural logs \cr 2. multiply the logged data
#' to creat a mirror image of the De distribution\cr 3. shift De values along
#' x-axis by the smallest x-value found to obtain only positive values \cr 4.
#' combine in quadrature the measurement error associated with each De value
#' with a relative error specified by sigmab \cr 5. apply the MAM to these data
#' \cr\cr When all calculations are done the results are then converted as
#' follows\cr\cr 1. subtract the x-offset \cr 2. multiply the natural logs by
#' -1 \cr 3. take the exponent to obtain the maximum dose estimate in Gy \cr\cr
#' \bold{Further documentation} \cr\cr Please see \code{\link{calc_MinDose}}.
#'
#' @param data \code{\linkS4class{RLum.Results}} or \link{data.frame}
#' (\bold{required}): for \code{data.frame}: two columns with De
#' \code{(data[,1])} and De error \code{(values[,2])}
#' @param sigmab \code{\link{numeric}} (\bold{required}): spread in De values
#' given as a fraction (e.g. 0.2). This value represents the expected
#' overdispersion in the data should the sample be well-bleached (Cunningham &
#' Walling 2012, p. 100).
#' @param log \code{\link{logical}} (with default): fit the (un-)logged three
#' parameter minimum dose model to De data
#' @param par \code{\link{numeric}} (with default): apply the 3- or
#' 4-parametric minimum age model (\code{par=3} or \code{par=4}).
#' @param bootstrap \code{\link{logical}} (with default): apply the recycled
#' bootstrap approach of Cunningham & Wallinga (2012).
#' @param init.values \code{\link{numeric}} (with default): starting values for
#' gamma, sigma, p0 and mu. Custom values need to be provided in a vector of
#' length three in the form of \code{c(gamma, sigma, p0)}.
#' @param plot \code{\link{logical}} (with default): plot output
#' (\code{TRUE}/\code{FALSE})
#' @param \dots further arguments for bootstrapping (\code{bs.M, bs.N, bs.h,
#' sigmab.sd}).  See details for their usage.
#' @return Please see \code{\link{calc_MinDose}}.
#' @section Function version: 0.3
#' @author Christoph Burow, University of Cologne (Germany) \cr Based on a
#' rewritten S script of Rex Galbraith, 2010 \cr
#' @seealso \code{\link{calc_CentralDose}}, \code{\link{calc_CommonDose}},
#' \code{\link{calc_FiniteMixture}}, \code{\link{calc_FuchsLang2001}},
#' \code{\link{calc_MinDose}}
#' @references Arnold, L.J., Roberts, R.G., Galbraith, R.F. & DeLong, S.B.,
#' 2009. A revised burial dose estimation procedure for optical dating of young
#' and modern-age sediments. Quaternary Geochronology 4, 306-325. \cr\cr
#' Galbraith, R.F. & Laslett, G.M., 1993. Statistical models for mixed fission
#' track ages. Nuclear Tracks Radiation Measurements 4, 459-470. \cr\cr
#' Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M.,
#' 1999. Optical dating of single grains of quartz from Jinmium rock shelter,
#' northern Australia. Part I: experimental design and statistical models.
#' Archaeometry 41, 339-364. \cr\cr Galbraith, R.F., 2005. Statistics for
#' Fission Track Analysis, Chapman & Hall/CRC, Boca Raton. \cr\cr Galbraith,
#' R.F. & Roberts, R.G., 2012. Statistical aspects of equivalent dose and error
#' calculation and display in OSL dating: An overview and some recommendations.
#' Quaternary Geochronology 11, 1-27. \cr\cr Olley, J.M., Roberts, R.G.,
#' Yoshida, H., Bowler, J.M., 2006. Single-grain optical dating of grave-infill
#' associated with human burials at Lake Mungo, Australia. Quaternary Science
#' Reviews 25, 2469-2474.\cr\cr \bold{Further reading} \cr\cr Arnold, L.J. &
#' Roberts, R.G., 2009. Stochastic modelling of multi-grain equivalent dose
#' (De) distributions: Implications for OSL dating of sediment mixtures.
#' Quaternary Geochronology 4, 204-230. \cr\cr Bailey, R.M. & Arnold, L.J.,
#' 2006. Statistical modelling of single grain quartz De distributions and an
#' assessment of procedures for estimating burial dose. Quaternary Science
#' Reviews 25, 2475-2502. \cr\cr Cunningham, A.C. & Wallinga, J., 2012.
#' Realizing the potential of fluvial archives using robust OSL chronologies.
#' Quaternary Geochronology 12, 98-106. \cr\cr Rodnight, H., Duller, G.A.T.,
#' Wintle, A.G. & Tooth, S., 2006. Assessing the reproducibility and accuracy
#' of optical dating of fluvial deposits.  Quaternary Geochronology 1, 109-120.
#' \cr\cr Rodnight, H., 2008. How many equivalent dose values are needed to
#' obtain a reproducible distribution?. Ancient TL 26, 3-10. \cr\cr
#' @examples
#'
#' ## load example data
#' data(ExampleData.DeValues, envir = environment())
#'
#' # apply the maximum dose model
#' calc_MaxDose(ExampleData.DeValues$CA1, sigmab = 0.2, par = 3)
#'
#' @export
calc_MaxDose<- function(
  data,
  sigmab,
  log=TRUE,
  par=3,
  bootstrap=FALSE,
  init.values,
  plot=TRUE,
  ...
){
  res<- calc_MinDose(data, sigmab, log, par, bootstrap, init.values, plot=FALSE, invert=TRUE, ...)
  res@originator<- "calc_MaxDose"
  if (plot) try(plot_RLum.Results(res, ...))

  invisible(res)
}
