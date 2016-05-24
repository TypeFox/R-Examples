##' Displays the probability of non-domination in the variable space. In dimension larger than two, projections in 2D subspaces are displayed.
##' @title Plot uncertainty
##' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
##' @param paretoFront (optional) matrix corresponding to the Pareto front of size \code{[n.pareto x n.obj]},
##' @param type type of uncertainty, for now only the probability of improvement over the Pareto front,
##' @param lower vector of lower bounds for the variables,
##' @param upper vector of upper bounds for the variables,
##' @param resolution grid size (the total number of points is \code{resolution^d}),
##' @param option optional argument (string) for n > 2 variables to define the projection type. The 3 possible values are "mean" (default), "max" and "min",
##' @param nintegpoints number of integration points for computation of mean, max and min values. 
##' @export
##' @importFrom grDevices grey.colors
##' @importFrom graphics image plot.new
##' @details Function inspired by the function \code{\link[KrigInv]{print_uncertainty}} and 
##' \code{\link[KrigInv]{print_uncertainty_nd}} from the package \code{\link[KrigInv]{KrigInv}}.
##' Non-dominated observations are represented with green diamonds, dominated ones by yellow triangles. 
##' @examples
##' \dontrun{ 
##' #---------------------------------------------------------------------------
##' # 2D, bi-objective function
##' #---------------------------------------------------------------------------
##' set.seed(25468)
##' n_var <- 2 
##' fname <- P1
##' lower <- rep(0, n_var)
##' upper <- rep(1, n_var)
##' res1 <- easyGParetoptim(fn=fname, lower=lower, upper=upper, budget=15, 
##' control=list(method="EHI", inneroptim="pso", maxit=20))
##' 
##' plot_uncertainty(res1$history$model, lower = lower, upper = upper)
##' 
##' #---------------------------------------------------------------------------
##' # 4D, bi-objective function
##' #---------------------------------------------------------------------------
##' set.seed(25468)
##' n_var <- 4
##' fname <- DTLZ2
##' lower <- rep(0, n_var)
##' upper <- rep(1, n_var)
##' res <- easyGParetoptim(fn=fname, lower=lower, upper=upper, budget = 40, 
##' control=list(method="EHI", inneroptim="pso", maxit=40))
##' 
##' plot_uncertainty(res$history$model, lower = lower, upper = upper, resolution = 31)
##' } 
plot_uncertainty <- function(model, paretoFront = NULL, type = "pn", lower, upper,
                              resolution = 51, option = "mean", nintegpoints = 400){
  if(is.null(resolution)) resolution <- 51
  n.obj <- length(model)
  observations <- matrix(0, model[[1]]@n, n.obj)
  for (i in 1:n.obj) observations[,i] <- model[[i]]@y
  if (is.null(paretoFront)) paretoFront <- t(nondominated_points(t(observations)))
  
  is.dom <- is_dominated(t(observations))
  
  ## 1D case
  if(model[[1]]@d == 1){
    xgrid <- matrix(seq(lower, upper, length.out = resolution), ncol = 1)
    pn <- prob.of.non.domination(paretoFront = paretoFront, model = model,
                                               integration.points = xgrid)
    plot(xgrid, pn, type = "l")
    points(model[[1]]@X, rep(0, model[[1]]@n), pch = 24, bg = "yellow")
  }
  
  ## 2D case
  if(model[[1]]@d == 2){
    x1 <- seq(lower[1], upper[1], length.out = resolution)
    x2 <- seq(lower[2], upper[2], length.out = resolution)
    Xgrid <- expand.grid(x1, x2)
    pn <- prob.of.non.domination(paretoFront = paretoFront, model = model,
                                 integration.points = Xgrid)
    filled.contour(matrix(pn, resolution), x = x1, y = x2, color.palette = graypalette,
                   main = "Probability of non-domination",
                   plot.axes={axis(1);axis(2);
                     points(model[[1]]@X[is.dom,], pch = 24, bg = "yellow");
                     points(model[[1]]@X[!is.dom,], pch = 23, bg = "green", cex=1.2)}
                   )
  }
  
  ## nD case
  if(model[[1]]@d > 2){
    if(resolution > 40)
      resolution <- 40
    d <- model[[1]]@d
    sub.d <- d - 2
    integration.tmp <- matrix(sobol(n = nintegpoints, dim = sub.d), 
                              ncol = sub.d)
    sbis <- c(1:resolution^2)
    numrow <- resolution^2 * nintegpoints
    integration.base <- matrix(c(0), nrow = numrow, ncol = sub.d)
    for (i in 1:sub.d) {
      col.i <- as.numeric(integration.tmp[, i])
      my.mat <- expand.grid(col.i, sbis)
      integration.base[, i] <- my.mat[, 1]
    }
    s <- seq(from = 0, to = 1, length = resolution)
    sbis <- c(1:nintegpoints)
    s.base <- expand.grid(sbis, s, s)
    s.base <- s.base[, c(2, 3)]
    prediction.points <- matrix(c(0), nrow = numrow, ncol = d)
    par(mfrow = c(d - 1, d - 1), mar=c(4,2,2,2))
    for (d1 in 1:(d - 1)) {
      for (d2 in 1:d) {
        if (d2 != d1) {
          if (d2 < d1) {
            plot.new()
          }
          else {
            myindex <- 1
            for (ind in 1:d) {
              if (ind == d1) {
                prediction.points[, ind] <- lower[ind] + 
                  s.base[, 1] * (upper[ind] - lower[ind])
                scale.x <- lower[ind] + s * (upper[ind] - 
                                               lower[ind])
              }
              else if (ind == d2) {
                prediction.points[, ind] <- lower[ind] + 
                  s.base[, 2] * (upper[ind] - lower[ind])
                scale.y <- lower[ind] + s * (upper[ind] - 
                                               lower[ind])
              }
              else {
                prediction.points[, ind] <- lower[ind] + 
                  integration.base[, myindex] * (upper[ind] - 
                                                   lower[ind])
                myindex <- myindex + 1
              }
            }
            pn <- prob.of.non.domination(paretoFront = paretoFront, model = model,
                                         integration.points = prediction.points)
            myvect <- matrix(pn, nrow = nintegpoints)
            if (option == "mean") {
              myvect <- colMeans(myvect)
            }
            else if (option == "max") {
              myvect <- apply(X = myvect, MARGIN = 2, FUN = max)
            }
            else if (option == "min") {
              myvect <- apply(X = myvect, MARGIN = 2, FUN = min)
            }
            else {
              myvect <- colMeans(myvect)
            }
            
            image(x = scale.x, y = scale.y, z = t(matrix(myvect, resolution)), 
                  col = grey.colors(10), xlab = "", 
                  ylab = "")
            contour(x = scale.x, y = scale.y, z = t(matrix(myvect, resolution)), 
                    add = TRUE)
            points(model[[1]]@X[is.dom,c(d2,d1)], pch = 24, bg = "yellow")
            points(model[[1]]@X[!is.dom,c(d2,d1)], pch = 23, bg = "green", cex=1.2)
          }
        }
      }
    }
    par(mfrow = c(1,1))
  }
}
