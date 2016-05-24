#' Generate hysteresis plots for theta parameter estimates
#'
#' @param hysteresis_output The list object output from the hysteresis function.
#' @param ... Additional arguments currently not supported.
#' @export
hysteresis_plot <- function(hysteresis_output,
                ...){
  #define colors
  UMASS_BLUE <- rgb(51,51,153,155,maxColorValue = 255)
  UMASS_RED <- rgb(153,0,51,155,maxColorValue = 255)
  num_thetas <- length(hysteresis_output[[1]]$theta_values)
  middle <- floor(num_thetas/2) + 1
  Theta = SE = Density = Order = NULL
  possible_structural_terms <- c("out2stars",
                                 "in2stars",
                                 "ctriads",
                                 "mutual",
                                 "ttriads",
                                 "edges")


  pretty_structural_terms <- c("Out Two-Stars",
                                 "In Two-Stars",
                                 "Cyclic Triads",
                                 "Mutual Dyads",
                                 "Transitive Triads",
                                 "Network Density")

  #loop over thetas
  for(i in 1:length(hysteresis_output)){
    modelFrame <- data.frame(Theta = c(hysteresis_output[[i]]$theta_values,
                                       rev(hysteresis_output[[i]]$theta_values)),
                             Density = hysteresis_output[[i]]$mean_densities,
                             SE = apply(hysteresis_output[[i]]$network_densities,
                                        2,sd),
                             Order = c(rep("Ascending",num_thetas),
                                       rep("Descending",num_thetas))
    )
    index <- which(possible_structural_terms == hysteresis_output[[i]]$term)

    data <- data.frame(modelFrame)

    zp1 <- ggplot2::ggplot(data, ggplot2::aes(colour = Order, group = Order)) +
      ggplot2::scale_color_manual(values = c(UMASS_BLUE,UMASS_RED))

    zp1 <- zp1 + ggplot2::geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)

    zp1 <- zp1 + ggplot2::geom_ribbon(data = data,
                 ggplot2::aes(x = Theta,
                     ymin = Density - SE*(-qnorm((1 - 0.95)/2)),
                     ymax = Density + SE*(-qnorm((1 - 0.95)/2))),
                alpha = 0.3,
                fill = "grey90") +

      ggplot2::geom_point(ggplot2::aes(x = Theta,
                                       y = Density,
                                       shape = Order,
                                       colour = Order),
                 size = 3) +
      ggplot2::scale_shape_manual(values = c( 16, 17))

    zp1 <- zp1  + ggplot2::theme_bw() +
      ggplot2::theme(legend.justification = c(1, 0),
                     legend.position = c(1, 0),
                     legend.key = ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"),
                     legend.title = ggplot2::element_blank()) +
      ggplot2::geom_hline(yintercept = hysteresis_output[[1]]$observed_density) +
      ggplot2::geom_vline(xintercept = hysteresis_output[[i]]$theta_values[middle])+
      ggplot2::ggtitle(pretty_structural_terms[index])

    print(zp1)
    Sys.sleep(1)
  }

}

# old confidence intervals
#     zp1 <- zp1 + ggplot2::geom_linerange(
#       ggplot2::aes(x = Theta,
#                    y = Density,
#                    ymin = Density - SE*(-qnorm((1-0.9)/2)),
#                    ymax = Density + SE*(-qnorm((1-0.9)/2))),
#                    lwd = 1)
#     zp1 <- zp1 + ggplot2::geom_linerange(
#       ggplot2::aes(x = Theta,
#                    y = Density,
#                    ymin = Density - SE*(-qnorm((1-0.95)/2)),
#                    ymax = Density + SE*(-qnorm((1-0.95)/2))),
#                    lwd = 0.5) +
