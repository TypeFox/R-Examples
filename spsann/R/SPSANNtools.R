#' Auxiliary tools
#' 
#' Auxiliary tools used in the optimization of sample configurations using spatial simulated annealing.
#' 
#' @param osc Object of class \code{OptimizedSampleConfiguration}.
#' 
#' @param at Point of the optimization at which the energy state should be returned. Available options: 
#' \code{"start"}, for the start, and \code{"end"}, for the end of the optimization. Defaults to 
#' \code{at = "end"}.
#' 
#' @param n Number of instances that should be returned. Defaults to \code{n = 1}.
#' 
#' @aliases SPSANNtools objSPSANN
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
# FUNCTION - RETRIEVE THE ENERGY STATE #########################################
#' @export
#' @rdname SPSANNtools
objSPSANN <- function (osc, at = "end", n = 1) {
  UseMethod( "objSPSANN" )
}
#' @export
objSPSANN.OptimizedSampleConfiguration <- 
  function (osc, at = "end", n = 1) {
    
    # Energy state at the start
    if (at == "start") res <- utils::head(osc$objective$energy, n)
    
    # Energy state at the end
    if (at == "end") res <- utils::tail(osc$objective$energy, n)
    
    return(data.frame(res, row.names = ""))
  }
# INTERNAL FUNCTION - COMPUTE ACCEPTANCE PROBABILITY ###########################
.acceptSPSANN <- 
  function (old.energy, new.energy, actual.temp) {
    accept <- min(1, exp((old.energy[[1]] - new.energy[[1]]) / actual.temp))
    accept <- floor(stats::rbinom(n = 1, size = 1, prob = accept))
    return (accept)
  }
# INTERNAL FUNCTION - PLOTTING #################################################
.spSANNplot <-
  function (energy0, energies, k, 
            # acceptance, accept_probs, 
            boundary, new_conf, conf0, y_max0, y.max, x_max0, x.max, best.energy, best.k, MOOP, wp) {
    
    # graphics::par(mfrow = c(1, 2))
    
    # PLOT ENERGY STATES
    grDevices::dev.set(grDevices::dev.prev())
    graphics::par(mar = c(5, 4, 4, 4) + 0.1)

    n <- ncol(energy0)
    l <- colnames(energy0)
    a <- rbind(energy0, energies)
    # col <- c("red", rep("black", n - 1))
    col <- c("red", grDevices::gray(seq(0, 0.5, length.out = n - 1)))
    graphics::plot(
      1, type = 'n', xlim = c(0, k), # ylim = c(0, max(a)), 
      ylim = c(min(a), max(a)), xlab = "jitter", ylab = "energy state")
    graphics::legend("topright", legend = l, lwd = 1, lty = 1:n, col = col)
    for(i in 1:ncol(a)) {
      graphics::lines(a[, i] ~ c(1:k), type = "l", lty = i, col = col[i])
    }
    graphics::lines(x = c(-k, 0), y = rep(energy0[1], 2), col = "red")
    graphics::lines(x = rep(best.k, 2), y = c(-500, best.energy[1]), col = "blue")
    
    # PLOT SAMPLE CONFIGURATION
    grDevices::dev.set(grDevices::dev.next())
    bb <- sp::bbox(boundary)
    if (methods::is(boundary, "SpatialPoints")) {
      sp::plot(x = boundary, pch = 20, cex = 0.1)
      # plot(x = boundary, pch = 20, cex = 0.1)
      # plot(boundary@coords, pch = 20, cex = 0.1)
    } else {
      sp::plot(x = boundary)
    }
    graphics::points(conf0[, 1], conf0[, 2], pch = 1, cex = 0.5, col = "lightgray")
    graphics::points(new_conf[, 1], new_conf[, 2], pch = 20, cex = 0.5)
    if (!missing(wp)) {
      graphics::points(new_conf[wp, 1], new_conf[wp, 2], pch = 20, cex = 1, col = "red")
    }
    
    # plot maximum shift in the x and y coordinates
    x <- c(bb[1, 1], bb[1, 2])
    y <- rep(bb[2, 1], 2) - 0.02 * y_max0
    graphics::lines(x = x, y = y, col = "gray", lwd = 12)
    
    y <- c(bb[2, 1], bb[2, 2])
    x <- rep(bb[1, 1], 2) - 0.02 * x_max0
    graphics::lines(y = y, x = x, col = "gray", lwd = 12)
    
    x <- c(bb[1, 1], bb[1, 1] + x.max)
    y <- rep(bb[2, 1], 2) - 0.02 * y_max0
    graphics::lines(x = x, y = y, col = "orange", lwd = 12)
    
    x <- rep(bb[1, 1], 2) - 0.02 * x_max0
    y <- c(bb[2, 1], bb[2, 1] + y.max)
    graphics::lines(y = y, x = x, col = "orange", lwd = 12)
    
    # plot labels for maximum shift in the x and y coordinates
    x <- bb[1, 1] + (bb[1, 2] - bb[1, 1]) / 2
    y <- bb[2, 1] - 0.02 * y_max0
    graphics::text(x = x, y = y, labels = "maximum shift in the X axis")
    
    x <- bb[1, 1] - 0.02 * x_max0
    y <- bb[2, 1] + (bb[2, 2] - bb[2, 1]) / 2
    graphics::text(y = y, x = x, srt = 90, labels = "maximum shift in the Y axis")
  }
# THE ORIGINAL spSANN FUNCTION #################################################
# .energyState <- 
#   function (fun, points, ...) {
#     if (missing(fun) || missing(points)) {
#       stop ("'fun' and 'points' are mandatory arguments")
#     }
#     return (do.call(fun, list(points, ...)))
#   }
# spatial simulated annealing
# spSANN <-
#   function (points, candi, x.max, x.min, y.max, y.min, fun, ...,
#             iterations = 10000, plotit = TRUE, boundary,
#             acceptance = list(initial = 0.99, cooling = iterations / 10),
#             stopping = list(max.count = 200), progress = TRUE, 
#             verbose = TRUE) {
#     if (plotit){
#       par0 <- graphics::par()
#     }
#     n_pts             <- dim(points)[1]
#     conf0       <- points
#     old_sys_config    <- conf0
#     energy0     <- .energyState(fun = fun, points = old_sys_config, ...)
#     old_energy_state  <- energy0
#     count             <- 0
#     best_energy_state <- Inf
#     energies     <- vector()
#     accept_probs      <- vector()
#     x_max             <- vector()
#     y_max             <- vector()
#     x_max0            <- x.max
#     y_max0            <- y.max
#     if (progress) {
#       pb <- txtProgressBar(min = 1, max = iterations, style = 3)
#     }
#     time0             <- proc.time()
#     for (k in 1:iterations) {
#       id <- sample(c(1:n_pts), 1)
#       new_conf <- spJitterFinite(old_sys_config, candi = candi,
#                                        x.max = x.max, x.min = x.min, 
#                                        y.max = y.max, y.min = y.min,
#                                        which.pts = id)
#       x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
#       x_max[k] <- x.max
#       y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
#       y_max[k] <- y.max
#       new_energy_state <- .energyState(fun=fun, points = new_conf, ...)
#       random_prob <- runif(1)
#       actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
#       accept_probs[k] <- actual_prob
#       if (new_energy_state <= old_energy_state) {
#         old_sys_config <- new_conf
#         old_energy_state <- new_energy_state
#         count <- 0
#       } else {
#         if (new_energy_state > old_energy_state & random_prob<=actual_prob) {
#           old_sys_config <- new_conf
#           old_energy_state <- new_energy_state
#           count <- count + 1
#           if (verbose) {
#             if (count == 1) {
#               cat("\n", count, "iteration with no improvement... p = ", 
#                   random_prob, "\n")
#             } else {
#               cat("\n", count, "iterations with no improvement... p = ", 
#                   random_prob, "\n")
#             }
#           }
#         } else {
#           new_energy_state <- old_energy_state
#           new_conf <- old_sys_config
#           count <- count + 1
#           if (verbose) {
#             if (count == 1) {
#               cat("\n", count, "iteration with no improvement... stops at",
#                   stopping$max.count, "\n")
#             } else {
#               cat("\n", count, "iterations with no improvement... stops at",
#                   stopping$max.count, "\n")
#             }
#           }
#         }
#       }
#       energies[k] <- new_energy_state
#       if (new_energy_state < best_energy_state / 1.0000001) {
#         best_k <- k
#         best_sys_config <- new_conf
#         best_energy_state <- new_energy_state
#         best_old_energy_state <- old_energy_state
#         old_sys_config <- old_sys_config
#       }
#       if (any(round(seq(1, iterations, 10)) == k)) {
#         if (plotit){
#           .spSANNplot(energy0, energies,k,acceptance,accept_probs, 
#                       boundary, new_conf, conf0, y_max0, y_max, 
#                       x_max0, x_max)
#         } 
#       }
#       if (count == stopping[[1]]) {
#         if (new_energy_state > best_energy_state * 1.000001) {
#           old_sys_config <- old_sys_config
#           new_conf <- best_sys_config
#           old_energy_state <- best_old_energy_state
#           new_energy_state <- best_energy_state
#           count <- 0
#           cat("\n", "reached maximum count with suboptimal system configuration\n")
#           cat("\n", "restarting with previously best system configuration\n")
#           if (count == 1) {
#             cat("\n", count, "iteration with no improvement... stops at",
#                 stopping[[1]], "\n")
#           } else {
#             cat("\n", count, "iterations with no improvement... stops at",
#                 stopping[[1]], "\n")
#           }
#         } else {
#           break
#         }
#       }
#       if (progress) {
#         setTxtProgressBar(pb, k)
#       }      
#     }
#     if (progress) {
#       close(pb)
#     }
#     if (plotit){
#       graphics::par(par0)
#     }
#     res <- new_conf
#     criterion <- c(energy0, energies)
#     a <- attributes(res)
#     a$energy.state <- criterion
#     running_time <- (proc.time() - time0) / 60
#     a$running.time <- running_time
#     attributes(res) <- a
#     cat("running time = ", round(running_time[3], 2), " minutes", sep = "")
#     return (res)
#   }
# End!
