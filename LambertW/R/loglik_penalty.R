#' @rdname loglik-LambertW-utils
#' @description \code{loglik_penalty} computes the penalty for transforming the
#'     data back to the input (see Goerg 2016). This penalty is independent of
#'     the distribution specified by \code{distname}, but only depends on
#'     \eqn{\tau}. If \code{type = "s"} then the penalty term exists if the
#'     distribution is non-negative (see \code{get_distname_family}) and
#'     \code{gamma >= 0}; otherwise, it returns \code{NA}.
#' 
#' @param is.non.negative logical; by default it is set to \code{TRUE} if the
#'     distribution is not a location but a scale family.
#' 
#' @export
loglik_penalty <- function(tau, y, type = c("h", "hh", "s"),
                           is.non.negative = FALSE) {
  
  stopifnot(is.numeric(y),
            is.logical(is.non.negative))
  
  type <- match.arg(type)
  yy <- y
    
  tau <- complete_tau(tau)
  zz <- normalize_by_tau(yy, tau)
  switch(type,
         h = {
           if (tau["delta"] == 0) {
             penalty <- 0
           } else {
             uu <- W_delta(zz, delta = tau["delta"])
             # penalty = sum(log(uu/zz) - log(1 + tau["delta * uu^2))
             penalty <-
                 sum(-tau["delta"]/2 * uu^2 - log(1 + tau["delta"] * uu^2))
           }
         },
         hh = {
           if (all(tau[grepl("delta", names(tau))] == 0)) {
             penalty <- 0
           } else {
             uu <- W_2delta(zz, delta = tau[c("delta_l", "delta_r")])
             ind <- (uu < 0)
             penalty <- sum(-tau["delta_l"]/2 * uu[ind]^2) +
                 sum(-tau["delta_r"]/2 * uu[!ind]^2) - 
                 sum(log(1 + tau["delta_l"] * uu[ind]^2)) -
                 sum(log(1 + tau["delta_r"] * uu[!ind]^2))
           }
         },
         s = {
          if (tau["gamma"] == 0) {
             penalty <- 0
           } else {
             if (is.non.negative) {
               penalty <- sum(log_deriv_W(tau["gamma"] * zz, branch = 0))
             } else {
               penalty <- NA
             }
           }
         })

  return(penalty)
} 
