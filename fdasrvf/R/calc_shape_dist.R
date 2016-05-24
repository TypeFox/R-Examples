#' Elastic Shape Distance
#'
#' Calculate elastic shape distance between two curves beta1 and beta2
#'
#' @param beta1 array describing curve1 (n,T)
#' @param beta2 array describing curve
#' @param mode Open ("O") or Closed ("C") curves
#' @return d geodesic distance
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' d = calc_shape_dist(beta[,,1,1],beta[,,1,4])
calc_shape_dist <- function(beta1, beta2, mode="O"){
    out = inverse_exp_coord(beta1, beta2, mode)

    return(out$d)
}
