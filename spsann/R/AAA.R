# Supress R CMD check note 'no visible binding for global variable ...' ########
# Source: http://stackoverflow.com/a/12429344/3365410
# These variables are generated within the functions created with 'autofun'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    names = c("covars.type", "sm", "n_pts", "n_cov", "n_candi", "wp", "conf0", 
              "pre.distri", "pop.prop", "id_fac", "id_num", "probs", "breaks",
              "pcm", "pop_prop", "covars_type", "x_max0", "x.min", "cellsize",
              "y_max0", "y.min", "x_max0", "best_sm", "best_old_sm", "i", 
              "res"))
}

# Import functions from default packages other than `base` #####################
# Source: http://stackoverflow.com/a/31314870/3365410
#' @importFrom stats cor quantile runif terms rbinom
#' @importFrom graphics axis hist legend lines mtext par plot points text
#' @importFrom utils setTxtProgressBar txtProgressBar head tail
#' @importFrom grDevices dev.set dev.prev dev.next dev.new gray
#' @importFrom methods is new slot
#' @importFrom rgeos gUnaryUnion
