#' Studies on the Effectiveness of Writing-to-Learn Interventions
#'
#' Results from 48 studies on the effectiveness of school-based writing-to-learn
#' interventions on academic achievement.
#'
#' @docType data
#'
#' @usage dat.bangertdrowns2004
#'
#' @format A data frame; for documentation, see \code{dat.bangertdrowns2004}
#' in Wolfgang Viechtbauer's R package \code{metafor}.
#'
#'
#' @keywords datasets
#'
#' @details This reproduced dataset and its documentation are credited to Wolfgang
#'  Viechtbauer and his \code{metafor} package (2010). Please see his package for
#'  details.
#'
#' @references Bangert-Drowns, R. L., Hurley, M. M., & Wilkinson, B. (2004). The
#' effects of school-based writing-to-learn interventions on academic achievement:
#' A meta-analysis. Review of Educational Research, 74, 29-58.
#'
#' Viechtbauer, W. (2010). Conducting meta-analysis in R with the metafor package.
#' Journal of Statistical Software, 36(3), 1-48.
#'
#' @source Bangert-Drowns, R. L., Hurley, M. M., & Wilkinson, B. (2004). The
#' effects of school-based writing-to-learn interventions on academic achievement:
#' A meta-analysis. Review of Educational Research, 74, 29-58.
#'
#' @examples
#' \dontrun{
#' dat.bangertdrowns2004
#'
#' # Extracting the effect sizes and sampling variances:
#' effect <- dat.bangertdrowns2004$yi
#' v <- dat.bangertdrowns2004$vi
#'
#' # The weight-function model with no mean model:
#' weightfunct(effect, v)
#'
#' # The weight-function model with a mean model:
#' weightfunct(effect, v, mods=~dat.bangertdrowns2004$info)
#' }
"dat.bangertdrowns2004"
