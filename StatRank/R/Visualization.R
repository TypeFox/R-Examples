# this is to get around issues with ggplot2 code not passing package checks
if(getRversion() >= "2.15.1")  utils::globalVariables(c("x", "y", "group", "column", "prob", "p", "m", "r"))

################################################################################
# Visualize MultiType
################################################################################

#' Multitype Random Utility visualizer
#' 
#' @param multitype.output output from a multitype fitter
#' @param min left boundary of graphed x-axis
#' @param max right boundary of graphed x-axis
#' @param names names of alternatives
#' @param ncol number of columns in final output
#' @return none
#' @export
#' @examples
#' library(ggplot2)
#' library(gridExtra)
#' Data.Tiny <- matrix(c(1, 2, 3, 3, 2, 1, 1, 2, 3), ncol = 3, byrow = TRUE)
#' multitype.output <- Estimation.RUM.MultiType.MLE(Data.Tiny, iter = 1, dist = "norm", ratio = .5)
#' names <- 1:3
#' #run the following code to make plots
#' #plots <- Visualization.MultiType(multitype.output, -2, 2, names, 3)
Visualization.MultiType <- function(multitype.output, min, max, names, ncol) {
  m <- dim(multitype.output$Mean)[2]

  process.MultiType <- function(output.of.multitype, m) {
    params <- list()
    for(i in 1:m) {
      params[[i]] <- list()
      params[[i]]$Mean <- output.of.multitype$Mean[,i] - mean(output.of.multitype$Mean)
      params[[i]]$SD <- sqrt(output.of.multitype$SD[,i])
      params[[i]]$Gamma <- output.of.multitype$Gamma[1,]
    }
    params
  }
  
  parameters <- process.MultiType(multitype.output, m)
  
  plots <- list()
  linesize <- 1
  for(i in 1:m) {
    means <- parameters[[i]]$Mean
    sds <- parameters[[i]]$SD
    gammas <- parameters[[i]]$Gamma
    get.density <- function(x, weights) sum(dnorm(c(x, x), mean = means, sd = sds) * gammas * weights)
    xs <- seq(min, max, by = 0.01)
    df <- data.frame(x = rep(xs, 3))
    df$y <- c(sapply(xs, function(x) get.density(x, c(1, 1))), sapply(xs, function(x) get.density(x, c(1, 0))), sapply(xs, function(x) get.density(x, c(0, 1))))
    df$group <- rep(c(1, gammas), rep(length(xs), 3))
    plots[[i]] <- ggplot(df, aes(x = x, y = y, linetype = factor(group), color = factor(group))) + 
      geom_line(size = linesize) + labs(title = names[i], x = NULL, y = NULL) + 
      geom_vline(xintercept = means[1], linetype = 2, size = linesize, color = "blue") + 
      geom_vline(xintercept = means[2], linetype = 6, size = linesize, color = "red") + 
      scale_color_manual(values=c("red", "blue", "black")) +
      scale_linetype_manual(values=c(6, 2, 1)) +
      theme(legend.position = "none")
  }
  do.call(gridExtra::grid.arrange, c(plots, ncol=ncol))
}

################################################################################
# Generating the Graphs
################################################################################

Generate.Pairwise.Matrix.Plot <- function(C.matrix, m, title = "", minprob = 0, maxprob = 0, reordering = NULL, labels = NULL) {
  base_size <- 15
  
  if(length(reordering) == 0) {
    reordering <- 1:m
  }
  
  if(length(labels) == 0) {
    labels = reordering
  }
  
  ggplot(C.matrix, aes(x = column, y = row)) +
    geom_tile(aes(fill = prob), colour = "white") + geom_text(aes(label = paste0(round(prob * 100), "%")), size = base_size * 0.3, color = "white") +
    scale_fill_gradient(limits = c(minprob, maxprob), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9), low = "darkblue", high = "darkred", labels = function (x) paste0(round(x * 100), "%"), guide = "legend") +
    labs(title = title, x = NULL, y = NULL, fill=NULL, size = base_size) + 
    scale_x_discrete(expand = c(0, 0), breaks = 2:m, labels = labels[-1]) +
    scale_y_reverse(breaks = 1:(m-1), labels = labels[-m]) +
    theme_bw(base_size = base_size) + 
    theme(legend.position = c(.15, .35), 
          axis.ticks = element_blank(), 
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(color = "gray50"),
          axis.text.x = element_text(color = "gray50", angle = 45, hjust=1),
          legend.text = element_text(color = "gray50")
    ) +
    coord_fixed() + 
    guides(fill = guide_legend(keywidth = 2.5, keyheight = 2.5))
}

#' Creates pairwise matrices to compare inference results with the empirical pairwise probabilities
#' 
#' @param Data.pairs datas broken into pairs
#' @param Parameters The Parameter element of a result from an Estimation function
#' @param get.pairwise.prob function that we use to generate the pairwise probability of beating
#' @param name.of.method names of the alternatives
#' @return none
#' @export
#' @examples
#' library(ggplot2)
#' library(gridExtra)
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' Parameters <- Estimation.PL.GMM(Data.Test.pairs, 5)$Parameters
#' PL.Pairwise.Prob <- function(a, b) a$Mean / (a$Mean + b$Mean)
#' Visualization.Pairwise.Probabilities(Data.Test.pairs, Parameters, PL.Pairwise.Prob, "PL")
Visualization.Pairwise.Probabilities <- function(Data.pairs, Parameters, get.pairwise.prob, name.of.method) {
  
  m <- max(Data.pairs[,c(1,2)])
  
  means <- rep(NA, m)
  
  # where it's not multitype
  if(is.null(Parameters[["Gamma"]])) for(i in 1:m) means[i] <- Parameters[[i]]$Mean
  
  # where it is multitype
  else for(i in 1:m) means[i] <- sum(Parameters[[i]]$Mean * Parameters[[i]]$Gamma)
  reordering <- order(-means)
  
  # calculate the empirical differences
  C.matrix.empirical <- generateC(Data.pairs, m)
  C.matrix.empirical.reordered <- C.matrix.empirical[reordering, reordering]
  C.empirical <- turn_matrix_into_table(C.matrix.empirical.reordered)
  
  # calculate the model differences
  C.matrix.model <- matrix(0, nrow = m, ncol = m)
  for(i in 1:m) for(j in 1:m) if(i != j) {
    C.matrix.model[i, j] <- get.pairwise.prob(Parameters[[i]], Parameters[[j]])
  }
  C.matrix.model.reordered <- C.matrix.model[reordering, reordering]
  C.model <- turn_matrix_into_table(C.matrix.model.reordered)
  
  
  minprob <- min(C.empirical['prob'], C.model['prob'])
  maxprob <- max(C.empirical['prob'], C.model['prob'])
  
  p1 <- Generate.Pairwise.Matrix.Plot(C.empirical, m, "Empirical", minprob, maxprob, reordering)
  p2 <- Generate.Pairwise.Matrix.Plot(C.model, m, "Model", minprob, maxprob, reordering)
  
  P <- C.empirical['prob'] / sum(C.empirical['prob'])
  Q <- C.model['prob'] / sum(C.model['prob'])
  
  kl.divergence <- round(sum(log(P/Q)*P), 10)
  
  name.of.data <- tail(strsplit(deparse(substitute(Data.pairs)), "\\.")[[1]], n=1)
  
  
  gridExtra::grid.arrange(p1, p2, nrow = 1)
  
}

makewireframe <- function(df, ...) {
  lattice::wireframe(z ~ x * y, data = df,
            drape = TRUE,
            colorkey = TRUE,
            scales = list(arrows = FALSE, draw = FALSE),
            xlab = NULL,
            ylab = NULL,
            zlab = NULL, 
            zlim = c(.4, 1),
            screen = list(z = -45, x = -60),
            par.settings = list(axis.line = list(col = "transparent")),
            ...)
}

Visualization.Surfaceplots <- function(Data.pairs, m, Estimate, pairwise.prob = NA, prior = 0, ...) {
  
  C.model <- generateC.model(Estimate, pairwise.prob, ...)[Estimate$order, Estimate$order]  
  C.empirical <- generateC(Data.pairs, m, prior = prior)[Estimate$order, Estimate$order]
  
  uppertriangle <- row(C.model) < col(C.model)
  C.model[!uppertriangle] <- .5
  C.empirical[!uppertriangle] <- .5
  
  C.model.df <- data.frame(z = as.numeric(C.model), x = as.numeric(row(C.model)), y = as.numeric(row(C.model)))
  C.empirical.df <- data.frame(z = as.numeric(C.empirical), x = as.numeric(row(C.empirical)), y = as.numeric(row(C.empirical)))
  C.diff.df <- data.frame(z = as.numeric(C.empirical - C.model), x = as.numeric(row(C.empirical)), y = as.numeric(row(C.empirical)))
  
   if (requireNamespace("lattice", quietly = TRUE)) {
     p1 <- makewireframe(C.model.df, main = "Model         ", col.regions = rainbow(600), at = seq(.4, 1, by = .001), zlim = c(min(C.model, C.empirical), 1))
     p2 <- makewireframe(C.empirical.df, main = "Empirical         ", col.regions = rainbow(600), at = seq(.4, 1, by = .001), zlim = c(min(C.model, C.empirical), 1))
     p3 <- makewireframe(C.diff.df, main = "Difference         ", col.regions = c(cm.colors(500)), zlim = c(-.25, .25), at = seq(-.25, .25, by = .01))
     gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
   }

}

#' RUMplot visualization
#' 
#' Creates marginal random utility density plots for each alternatives given an 
#' Estimation object for a PL or Nonparameteric model
#' 
#' @param RUM choice of Exponential, Gumbel, or Nonparametric
#' @param Estimate fitted RUM object
#' @param min minimum x value to display
#' @param max maximum x value to display
#' @param ncol number of columns in the visualization
#' @param names names of alternatives
#' @export
#' @examples
#' library(ggplot2)
#' library(gridExtra)
#' Data.Tiny <- matrix(c(1, 2, 3, 3, 2, 1, 1, 2, 3), ncol = 3, byrow = TRUE)
#' Estimate <- Estimation.PL.GMM(Breaking(Data.Tiny, method = "full"), m = 3)
#' Visualization.RUMplots("Exponential", Estimate, names = 1:3)
Visualization.RUMplots <- function(RUM = "Exponential", Estimate = NA, min = -5, max = 5, ncol= 5, names = NA) {
  if(RUM == "Exponential") {
    Visualization.RUMplots.Exponential(Estimate, min, max, ncol, names)
  } else if(RUM == "Gumbel") {
    Visualization.RUMplots.Gumbel(Estimate, min, max, ncol, names)
  } else if(RUM == "Nonparametric") {
    Visualization.RUMplots.Nonparametric(Estimate, max, ncol, names)
  } else {
    stop(paste("RUM", RUM, "not recognized.", sep = " "))
  }
}
  
Visualization.RUMplots.Exponential <- function(Estimate, min, max, ncol, names) {
  m <- Estimate$m
  x <- seq(min, max, by = (max - min) / 500)
  plots <- llply(1:m, function(i) qplot(x, dexp(x, rate = 1/Estimate$Mean[i]), geom = "line") + labs(y = NULL, x = NULL, title = names[i]) + geom_vline(xintercept = Estimate$Mean[i], linetype = "dashed")) #+ scale_y_continuous(limits = c(0, 20))) 
  do.call(gridExtra::grid.arrange, c(plots, ncol=ncol))
}

Visualization.RUMplots.Gumbel <- function(Estimate, min, max, ncol, names) {
  m <- Estimate$m
  x <- seq(min, max, by = (max - min) / 500)
  means <- log(1/Estimate$Mean)
  plots <- llply(1:m, function(i) qplot(x, dgumbel(x, location = means[i]), geom = "line", size = I(1)) + labs(y = NULL, x = NULL, title = names[i]) + geom_vline(xintercept = means[i], linetype = "dashed")) #+ scale_y_continuous(limits = c(0, 20))) 
  do.call(gridExtra::grid.arrange, c(plots, ncol=ncol))

}

Visualization.RUMplots.Nonparametric <- function(Estimate, ymax, ncol, names) {
  x.star <- seq(0, 1, len = 512)
  m <- Estimate$m
  plots <- llply(1:m, function(i) qplot(x.star, Estimate$rum.densities[[i]], geom = "line", size = I(1) ) + labs(y = NULL, x = NULL, title = names[i]) + scale_y_continuous(limits = c(0, ymax)))
  do.call(gridExtra::grid.arrange, c(plots, ncol=ncol))
}

#' RPD Visualization
#' 
#' Creates histograms of the empriical rank position distribution for each alternative
#' in rank data
#' 
#' @param Data full, top partial, or sub partial data
#' @param ymax maximum value of density to show on graph
#' @param ncol number of columns visualization is displayed in
#' @param names names of alternatives
#' @export
#' @examples
#' library(ggplot2)
#' library(gridExtra)
#' data(Data.Test)
#' Visualization.Empirical(Data.Test, 0.5)
Visualization.Empirical <- function(Data, ymax, ncol = 5, names = NA) {
  m <- ncol(Data)
  empirical.plots <- list()
  for(i in 1:m) {
    temp <- data.frame(x = 1:m, y = colSums(Data==i))
    temp$p <- temp$y / sum(temp$y)
    empirical.plots[[i]] <- ggplot(temp, aes(x, p)) + geom_bar(stat = "identity") +
      scale_x_reverse(breaks = seq(1, m, by = 1), labels = seq(1, m, by = 1)) +
      ggtitle(names[i]) + labs(x = NULL, y = NULL) + scale_y_continuous(limits = c(0, ymax))
  }
  
  do.call(gridExtra::grid.arrange, c(empirical.plots, ncol=ncol))
}