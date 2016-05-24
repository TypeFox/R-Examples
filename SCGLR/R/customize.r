#' Plot customization
#' @name customize
#' @description
#' Parameters used to choose what to plot and how. These parameters are given to \code{\link{plot.SCGLR}} and \code{\link{pairs.SCGLR}}.
#' @details
#' Parameter name can be abbreviated (e.g. \code{pred.col} will be understood as predictors.color).\cr
#' Options can be set globally using \code{options("plot.SCGLR")}. It will then provide default values that can be further overriden by giving explicit parameter value.
#' \tabular{ll}{
#'    \bold{\emph{parameter name}}    \tab \bold{\emph{type (default value). Description.}}\cr
#'    \code{title}                    \tab string (NULL). Main title of plot (override built-in).\cr
#'    \code{labels.auto}              \tab logical (TRUE). Should covariate or predictor labels be aligned with arrows.\cr
#'    \code{labels.offset}            \tab numeric (0.01). Offset by which labels should be moved from tip of arrows.\cr
#'    \code{labels.size}              \tab numeric (1). Relative size for labels. Use it to globally alter label size.\cr
#'    \code{expand}                   \tab numeric (1). Expand factor for windows size. Use it for example to make room for clipped labels.\cr
#'    \code{threshold}                \tab numeric. All covariates and/or predictors whose sum of square correlations with the two components of the plane lower than this threshold will be ignored.\cr
#'
#'    \bold{\code{observations}}      \tab logical (FALSE). Should we draw observations.\cr
#'    \code{observations.size}        \tab numeric (1). Point size. \cr
#'    \code{observations.color}       \tab character ("black"). Point color.\cr
#'    \code{observations.alpha}       \tab numeric (1). Point transparency.\cr
#'    \code{observations.factor}      \tab logical (FALSE). Paint observations according to factor (specify factor).\cr
#'
#'    \bold{\code{predictors}}        \tab logical or array of characters (FALSE). Should we draw predictors and optionally which one (TRUE means all).\cr
#'    \code{predictors.color}         \tab string ("red"). Base color used to draw predictors.\cr
#'    \code{predictors.alpha}         \tab numeric (1). Overall transparency for predictors (0 is transparent, 1 is opaque). \cr
#'    \code{predictors.arrows}        \tab logical (TRUE). Should we draw arrows for predictors.\cr
#'    \code{predictors.arrows.color}  \tab string (predictors.color). Specific color for predictor arrows.\cr
#'    \code{predictors.arrows.alpha}  \tab numeric (predictors.alpha). Transparency for predictor arrows.\cr
#'    \code{predictors.labels}        \tab logical (TRUE). Should we draw labels for predictors.\cr
#'    \code{predictors.labels.color}  \tab string (predictors.color). Specific color for predictor labels.\cr
#'    \code{predictors.labels.alpha}  \tab numeric (predictors.alpha). Transparency for predictor labels.\cr
#'    \code{predictors.labels.size}   \tab numeric (labels.size). Specific size for predictor labels.\cr
#'    \code{predictors.labels.auto}   \tab logical (labels.auto). Should predictor labels be aligned with arrows.\cr
#'
#'    \bold{\code{covariates}}        \tab logical or array of characters (TRUE). Should we draw covariates and optionally which one (TRUE means all).\cr
#'    \code{covariates.color}         \tab string ("black"). Base color used to draw covariates. \cr
#'    \code{covariates.alpha}         \tab numeric (1). Overall transparency for covariates (0 is transparent, 1 is opaque). \cr
#'    \code{covariates.arrows}        \tab logical (TRUE). Should we draw arrows for covariates.\cr
#'    \code{covariates.arrows.color}  \tab string (covariates.color). Specific color for covariate arrows.\cr
#'    \code{covariates.arrows.alpha}  \tab numeric (covariates.alpha). Transparency for covariate arrows.\cr
#'    \code{covariates.labels}        \tab logical (TRUE). Should we draw labels for predictors.\cr
#'    \code{covariates.labels.color}  \tab string (covariates.color). Specific color for predictor labels.\cr
#'    \code{covariates.labels.alpha}  \tab numeric (covariates.alpha). Transparency for covariate labels.\cr
#'    \code{covariates.labels.size}   \tab numeric (labels.size). Specific size for covariate labels.\cr
#'    \code{covariates.labels.auto}   \tab logical (labels.auto). Should covariate labels be aligned with arrows. \cr
#'
#'    \bold{\code{factor}}            \tab logical or character (FALSE). Should we draw a factor chosen among additionnal variables (TRUE mean first one).\cr
#'    \code{factor.points}            \tab logical (TRUE). Should symbol be drawn for factors.\cr
#'    \code{factor.points.size}       \tab numeric (4). Symbol size.\cr
#'    \code{factor.points.shape}      \tab numeric (13). Point shape.\cr
#'    \code{factor.labels}            \tab logical (TRUE). Should factor labels be drawn.\cr
#'    \code{factor.labels.color}      \tab string ("black"). Color used to draw labels.\cr
#'    \code{factor.labels.size}       \tab numeric (labels.size). Specific size for factor labels.\cr
#' }
#' @examples \dontrun{
#' # setting parameters
#' plot(genus.scglr)
#' plot(genus.scglr, predictors=TRUE)
#' plot(genus.scglr, predictors=TRUE, pred.arrows=FALSE)
#' 
#' # setting global style
#' options(plot.SCGLR=list(predictors=TRUE, pred.arrows=FALSE))
#' plot(genus.scglr)
#' 
#' # setting custom style
#' myStyle <- list(predictors=TRUE, pred.arrows=FALSE)
#' plot(genus.scglr, style=myStyle)
#' }
NULL
