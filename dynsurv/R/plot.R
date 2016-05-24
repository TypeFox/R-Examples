##############################################################################
##
##   R package dynsurv by Xiaojing Wang, Jun Yan, and Ming-Hui Chen
##   Copyright (C) 2011
##
##   This file is part of the R package dynsurv.
##
##   The R package dynsurv is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package dynsurv is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package dynsurv. If not, see <http://www.gnu.org/licenses/>.
##
##############################################################################

##############################################################################
# Plot coefficient data frame
##############################################################################
plotCoef <- function(object, smooth=FALSE, ...) {

    p <- ggplot(data=object, aes_string(x="Time"))

    if (!smooth)
        p <- p + geom_step(aes_string(y="Mid"), direction="vh") +
            geom_step(aes_string(y="High"), direction="vh", linetype=2) +
            geom_step(aes_string(y="Low"), direction="vh", linetype=2)
    else
       p <- p + geom_line(aes_string(y="Mid")) +
           geom_line(aes_string(y="High"), linetype=2) +
           geom_line(aes_string(y="Low"), linetype=2)

    if (length(levels(factor(object$Model))) == 1)
        p <- p + facet_wrap(~ Cov, scales="free_y")
    else
        p <- p + facet_grid(Cov ~ Model, scales="free_y")

    p <- p + ylab("Coefficient") +
        theme(plot.margin=unit(rep(0, 4), "lines"))

    p
}


##############################################################################
# Plot iteration jump data frame returned by jump.bayesCox
##############################################################################
plotJumpTrace <- function(object, ...) {
    p <- ggplot(data=object, aes_string(x="Iter", y="Count")) +
        geom_line(size=0.1, alpha=0.6) +
        facet_wrap(~ Cov) +
        xlab("Iteration") + ylab("Pieces of Coefficient") +
        theme(plot.margin=unit(rep(0, 4), "lines"))

    p
}

plotJumpHist <- function(object, ...) {
  ## p <- ggplot(data=object, aes(x=factor(Count))) +
  p <- ggplot(data=object, aes_string(x="Count")) +
    ##   stat_bin(aes(y =..count../sum(..count..))) +
    stat_bin(aes_string(y = "..density..")) +
      facet_wrap(~ Cov) +
        xlab("Pieces of Coefficient") + ylab("Relative Frequency") +
          theme(plot.margin=unit(rep(0, 4), "lines"))
  
  p
}

##############################################################################
# Plot the latent variance nu from the bayesCox model
##############################################################################
plotNu <- function(object, ...) {
  cnt <- "..density.."
  ## p <- ggplot(data=object, aes_string(x="Value")) +
  p <- ggplot(data=object, aes_string(x="Value")) +
    ##  stat_bin(aes(y =..count../sum(..count..))) +
    stat_bin(aes_string(y = cnt)) +
      xlab("Nu") + ylab("Relative Frequency") +
        theme(plot.margin=unit(rep(0, 4), "lines"))
  
  if (length(levels(factor(object$Model))) == 1)
    p <- p + facet_wrap(~ Cov)
  else
    p <- p + facet_grid(Cov ~ Model)

  p
}
