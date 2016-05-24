multi_tpfit <-
function(data, coords, method = "ml", tolerance = pi/8, rotation = NULL, max.it = 9000, mle = "avg", ...) {
  # Estimation for matrix of transition rates
  #
  #       data vector of data
  #     coords coordinates matrix
  #     method estimation method c("ml", "ils", "me")
  #      scale 
  #  tolerance angle tolerance (in radians)
  #   rotation vector of rotation angles (in radians)
  #     max.it maximum number of iterations for the optimization (used only for the 'me' method)
  #        mle argument to pass to the function mlen (not used for the 'ils' method)
  #        ... further arguments to pass to tpfit_ils function, such as:
  #    #     * max.dist maximum distance for counting
  #    #     *  mpoints number of lags
  #    #     *        q constant greater than one controlling the growth of rho
  #    #     *     echo logical value to print the optimization output
  #    #     *      ... further arguments to pass to nlminb function
  #    #     *    mtpfit mulit_tpfit object for a further optimization
  
  if (method == "ils") return(multi_tpfit_ils(data = data, coords = coords, tolerance = tolerance, rotation = rotation, ...))
  if (method == "me") return(multi_tpfit_me(data, coords, tolerance, max.it, rotation, mle))
  if (method != "ml") warning("Estimation method not recognized. Mean length method (\"ml\") set by default.")

  return(multi_tpfit_ml(data, coords, tolerance, rotation, mle))
}
