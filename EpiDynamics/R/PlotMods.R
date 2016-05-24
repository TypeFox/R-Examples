#' Plot results of capm model functions
#' @description Plot results of EpiDynamics' functions.
#' @param model.out output of aEpiDynamics' function.
#' @param variables column index for the variables in model.out to be plotted.
#' @param x.label string with the name of x axis.
#' @param y.label string with the name of y axis.
#' @param legend.title string with the legend title.
#' @param line.size scalar to define the thick of the lines (points for bifurcations) to be plotted.
#' @param text.size scalar to define the size of axis texts and titles.
#' @param grid logical to indicate if each variable must be plotted in a separated panel.
#' @param bifur logical to indicate if \code{model.out} represent a bifurcation.
#' @export
#' @examples
#' # Parameters and initial conditions.
#' parameters <- list(beta0 = 17 / 13, beta1 = 0.1, gamma = 1 / 13,
#'                    omega = 2 * pi / 365, mu = 1 / (50 * 365))
#' 
#' initials <- c(S = 1 / 17, I = 1e-4, 
#'               R = 1 - 1 / 17 - 1e-4)
#' 
#' # Solve the system.
#' sir.sinusoidal.forcing <- SIRSinusoidalForcing(pars = parameters, 
#'                                                init = initials, 
#'                                                time = 0:(60 * 365))
#' PlotMods(sir.sinusoidal.forcing)                                          
#'                                                
#' # Solve bifurcation dynamics for 20 years.
#' # If max(time) < 3650, bifurcation dynamics are solved for 3650 time-steps.
#' parameters2 <- list(beta0 = 17 / 13, beta1 = seq(0.001, 0.251, by = 0.001),
#'                    gamma = 1 / 13, omega = 2 * pi / 365, mu = 1 / (50 * 365))
#' # Uncomment the following lines:
#' # bifur <- SIRSinusoidalForcing(pars = parameters2, 
#' #                               init = initials,
#' #                               time = 0:(20 * 365))
#' # PlotMods(bifur, bifur = TRUE)
#' 
PlotMods <- function(model.out = NULL, variables = NULL, x.label = NULL, y.label = NULL, legend.title = 'variable', line.size = 1, text.size = 14,  grid = TRUE, bifur = FALSE) {
  value <- variable <- x <- y <- NULL
  if (bifur == FALSE) {
    if (is.null(variables)) {
      if (is.data.frame(model.out)) {
        model <- model.out
      } else {
        model <- model.out$results
      }
      variables <- 1:ncol(model)
    } else {
      if (is.character(variables)) {
        variables = c('time', variables)
      } else {
        variables = c(1, variables)
      }
      if (is.data.frame(model.out)) {
        model <- model.out
      } else {
        model <- model.out$results
      }
    }
    if (is.null(x.label)) {
      x.label <- 'time'
    }
    if (is.null(y.label)) {
      y.label <- 'value'
    }
    if (grid) {
      vplayout <- function(x, y) {
        viewport(layout.pos.row = x, layout.pos.col = y)
      } 
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(length(variables) - 1, 1)))
      for (i in 1:(length(variables) - 1)) {
        tmp <- model[ , c(variables[1], variables[-1][i])]
        var <- names(tmp)[2]
        names(tmp)[2] <- 'value'
        print(ggplot(tmp, aes(time, value)) +
                geom_line(size = line.size, color = i) +
                ylab(var) + xlab(x.label),
              vp = vplayout(i, 1))  
      }
    } else {
      model <- melt(model[ , variables], id = 'time')
      print(ggplot(model, aes(time, value, color = variable)) +
        geom_line(size = line.size) +
        scale_colour_discrete(name = legend.title) +
        xlab(x.label) + ylab(y.label))
    }
  }
  if (bifur) {
    if (is.null(x.label)) {
      x.label <- names(model.out)[1]
    }
    if (is.null(y.label)) {
      y.label <- 'Level of infection'
    }
    plt <- data.frame(model.out[ , 1],
                      matrix(as.matrix(model.out[ , -1]), ncol = 1))
    names(plt) <- c('x', 'y')
    if (is.null(x.label)) {
      x.label <- names(model.out)[1]
    }
    ggplot(plt, aes(x, y)) +
      geom_point(size = line.size + 2) +
      xlab(x.label) + ylab(y.label) +
    theme(axis.text.x = element_text(size = text.size),
          axis.text.y = element_text(size = text.size),
          axis.title.x = element_text(size = text.size + 1),
          axis.title.y = element_text(size = text.size + 1))
    
  }
}
