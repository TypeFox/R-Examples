tpfit <-
function(data, coords, direction, method = "ml", tolerance = pi/8, max.it = 9000, mle = "avg", ...) {
  # Estimation for matrix of transition rates
  #
  #       data vector of data
  #     coords coordinates matrix
  #  direction vector (or versor) of choosen direction
  #     method estimation method c("ml", "ils", "me")
  #  tolerance angle tolerance (in radians)
  #     max.it maximum number of iterations for the optimization (used only for the 'me' method)
  #        mle argument to pass to the function mlen (not used for the 'ils' method)
  #        ... further arguments to pass to tpfit_ils function, such as:
  #    #     * max.dist maximum distance for counting
  #    #     *  mpoints number of lags
  #    #     *        q constant greater than one controlling the growth of rho
  #    #     *     echo logical value to print the optimization output
  #    #     *      ... further arguments to pass to nlminb function
  #    #     *    tpfit tpfit object for a further optimization

  if (method == "ils") return(tpfit_ils(data = data, coords = coords, direction = direction, tolerance = tolerance, ...))
  if (method == "me") return(tpfit_me(data, coords, direction, tolerance, max.it, mle))
  if (method != "ml") warning("Estimation method not recognized. Mean length method (\"ml\") set by default.")

  return(tpfit_ml(data, coords, direction, tolerance, mle))
}
