#' Basic statistical functions for ff
#' 
#' Basic statistical functions for \code{\link{ff}} vectors and \code{\link{ffdf}} data.frames. 
#' The aim of ffbase is to make working with ff vectors and ffdf data.frame a bit easier.
#'
#' @section Basic operations:
#' \tabular{ll}{
#'    \code{\link{cut.ff}} \tab cut a \code{ff} vector. \cr
#'    \code{\link{c.ff}} \tab concatenate \code{ff} vectors. \cr
#'    \code{\link{unique}} \tab unique for a \code{ff} vector and \code{ffdf}. \cr
#'    \code{\link{duplicated}} \tab duplicated for a \code{ff} vector and \code{ffdf}. \cr
#'    \code{\link{ffmatch}} \tab match for a 2 \code{ff} vectors. \cr
#'    \code{\link{ffdfmatch}} \tab match for 2 \code{ffdf} objects. \cr
#'    \code{\%in\%} \tab \%in\% operator for a \code{ff} vector and \code{ffdf}. \cr
#'    \code{\link{is.na.ff}} \tab is.na for a \code{ff} vector. \cr
#'    \code{+, -, *, /, ^, \%\%, \%/\%} \tab operators for arithmetic on \code{ff} vector. \cr
#'    \code{==, !=, <, <=, >=, >, &, |, !} \tab compare & logic operators for working with \code{ff} vectors. \cr
#'    \code{abs, sign, sqrt, ceiling, floor, trunc, round, signif} \tab Math operators for working on \code{ff} vectors. \cr
#'    \code{log, log10, log2, log1p, exp, expm1} \tab Math operators for working on \code{ff} vectors. \cr
#'    \code{acos, acosh, asin, asinh, atan, atanh} \tab Math operators for working on \code{ff} vectors. \cr
#'    \code{cos, cosh, sin, sinh, tan, tanh} \tab Math operators for working on \code{ff} vectors. \cr
#'    \code{gamma, lgamma, digamma, trigamma} \tab Math operators for working on \code{ff} vectors. \cr
#' }
#' 
#' @section Selections:
#' \tabular{ll}{
#'    \code{\link{subset.ffdf}} \tab subset a \code{ffdf}. \cr
#'    \code{\link{transform.ffdf}} \tab create a new \code{ffdf} based on an existing \code{ffdf} \cr
#'    \code{\link{with.ffdf}} \tab create a ff vector based on columns of an existing \code{ffdf} \cr
#'    \code{\link{within.ffdf}} \tab create a \code{ffdf} data.frame based on columns of an existing \code{ffdf} \cr
#'    \code{\link{ffwhich}} \tab create a \code{ff} integer vector based on a logical expression \cr
#' }
#' 
#' @section Aggregations:
#' \tabular{ll}{
#'    \code{\link{hist.ff}} \tab Calculate a histogram for \code{ff} vector. \cr
#'    \code{\link{quantile.ff}} \tab Get quantiles for \code{ff} vector. \cr
#'    \code{\link{sum.ff}} \tab sum for a \code{ff} vector. \cr
#'    \code{\link{mean.ff}} \tab (trimmed) mean for a \code{ff} vector. \cr
#'    \code{\link{all.ff}} \tab all for logical \code{ff} vector. \cr
#'    \code{\link{min.ff}} \tab min for \code{ff} vector. \cr
#'    \code{\link{max.ff}} \tab max for \code{ff} vector. \cr
#'    \code{\link{cumsum.ff}} \tab cumsum for \code{ff} vector. \cr
#'    \code{\link{cumprod.ff}} \tab cumprod for \code{ff} vector. \cr
#'    \code{\link{range.ff}} \tab range for \code{ff} vector. \cr
#'    \code{\link{table}} \tab table for \code{ff} vectors. \cr
#'    \code{\link{tabulate.ff}} \tab tabulate for \code{ff} vectors. \cr
#'    \code{\link{ffdfdply}} \tab Split, group and aggregate for \code{ffdf} operations. \cr
#' }
#'
#' @section Miscellaneous:
#' \tabular{ll}{
#'    \code{\link{ffordered}} \tab Add a sorted index to a \code{ff} vector.\cr
#'    \code{\link{save.ffdf}} \tab Save a \code{ffdf} in a directory with its containing \code{ff} columns.\cr
#'    \code{\link{load.ffdf}} \tab Loads a \code{ffdf} from a directory\cr
#'    \code{\link{pack.ffdf}} \tab Packs ffdf data.frames into a zip or tar file\cr
#'    \code{\link{unpack.ffdf}} \tab Unpacksdata.frames from a zip or tar file\cr
#'    \code{\link{ffappend}} \tab Append data to a \code{ff} vector.\cr
#'    \code{\link{ffdfappend}} \tab Append data to a \code{ffdf}.\cr
#'    \code{\link{merge.ffdf}} \tab Merge two \code{ffdf} objects. \cr
#'    \code{\link{ffmatch}} \tab match two ff vectors \cr
#'    \code{\link{ffdfmatch}} \tab match two ffdf data.frames \cr
#'    \code{\link{laf_to_ffdf}} \tab Import csv and fixed width files through package LaF. \cr
#' }
#'
#' @example ../examples/pkg.R
#' @name ffbase-package 
#' @aliases ffbase ffbase-package
#' @import bit
#' @importFrom ff ffdf ff is.ff is.ffdf fforder filename as.ff is.factor.ff vmode recodeLevels 
#' @importFrom ff is.open ffdforder appendLevels ffvecapply 'filename<-' maxffmode
#' @importFrom ff arrayIndex2vectorIndex nrow<- .rambytes set.ff hi ffindexordersize ffindexget
#' @importFrom ff ffsort .vNA .vmax .vmin .vmode finalizer<- ffindexorder .vimplemented 
#' @importFrom graphics image
#' @importFrom stats gaussian na.omit runif
#' @importFrom utils head tail tar untar unzip zip
#' @docType package 
{}
