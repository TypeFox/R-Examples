##//////////////////////////////////////////////////////////////////////////////
##//zzz.R
##//////////////////////////////////////////////////////////////////////////////
##
##==============================================================================
##author: R Luminescence Package Team
##organisation:
##version.: 0.2.1
##date: 2013-11-10
##==============================================================================
# Set namespace .LuminescenceEnv ------------------------------------------
.LuminescenceEnv <- new.env(parent = emptyenv())

# Assign variables to Namespace -------------------------------------------
##variable col to define colours in the functions for output
assign("col",
       unlist(colors())[c(261,552,51,62,76,151,451,474,654,657,100,513,23,612,129,27,551,393,80,652,555)],
       pos = ".LuminescenceEnv",
       envir = .LuminescenceEnv)



##==============================================================================
##on Attach
.onAttach <- function(libname,pkgname){

  ##set startup message
  try(packageStartupMessage(paste("Welcome to the R package Luminescence version ",
                              packageDescription(pkg="Luminescence")$Version,
                              " [Built: ",
                              strsplit(packageDescription(pkg="Luminescence")$Packaged, ";")[[1]][1],
                             "]", sep=""),
                            "\n",
                            get_Quote()), silent=TRUE)
}

##==============================================================================
# DO NOT TOUCH! -----------------------------------------------------------



#' sTeve - sophisticated tool for efficient data validation and evaluation
#'
#' This function provides a sophisticated routine for comprehensive
#' luminescence dating data analysis.
#'
#' This amazing sophisticated function validates your data seriously.
#'
#' @param n_frames \code{\link{integer}} (with default): n frames
#' @param t_animation \code{\link{integer}} (with default): t animation
#' @param n.tree \code{\link{integer}} (with default): How many trees do you
#' want to cut?
#' @param type \code{\link{integer}} (optional): Make a decision: 1, 2 or 3
#' @return Validates your data.
#' @note This function should not be taken too seriously.
#' @author R Luminescence Team, 2012-2013
#' @seealso \link{plot_KDE}
#' @references #
#' @keywords manip
#' @examples
#'
#' ##no example available
#'
#' @export
sTeve<- function(n_frames = 10, t_animation = 2, n.tree = 7, type) {

  ## allow new overlay plot
  par(new = TRUE)

  ## infer month of year
  month <- as.numeric(strsplit(x = as.character(Sys.Date()), split = "-")[[1]][2])

  ## select showtime item based on month or user-defined type
  if(missing(type) == TRUE) {
    if(month >= 1 & month <= 3) {
      type <- 1
    } else if(month >3 & month <= 11) {
      type <- 2
    } else if(month > 11 & month <= 12) {
      type <- 3
    }
  }



  if(type == 1) {
    ## SHOWTIME OPTION 1
    Sys.sleep(5)
    shape::emptyplot()
    shape::filledrectangle(wx = 0.9, wy = 0.4,
                    mid = c(0.5, 0.5),
                    lcol ="red",
                    lwd=1,
                    col=0,
                    angle = 45)

    text(x=0.5, y=0.5,
         labels="NOT FUNNY",
         cex=2,
         col="red",
         font=2,
         srt=45)
  } else if(type == 2) {

    ## SHOWTIME OPTION 2
    plot(NA, xlim = c(0, 10),
         ylim = c(0, 10),
         main = "",
         xlab = "",
         ylab = "",
         axes = FALSE,
         frame.plot = FALSE)

    n_frames <- n_frames
    t_animation <- t_animation

    dt <- t_animation / n_frames
    x1 <- seq(0, 10, length.out = n_frames)
    y1 <- rep(1.5, n_frames)
    r1 <- 0.5

    x2 <- seq(0, 16, length.out = n_frames)
    y2 <- rep(8.5, n_frames)
    r2 <- 0.5

    x4 <- seq(11, 0, length.out = n_frames)
    y4 <- rep(5, n_frames)
    r4 <- 0.5

    # set angles for each step of mouth opening
    angles_mouth <- rep(c(0.01, 0.25, 0.5, 0.25),
                        length.out = n_frames)

    for(i in 1:n_frames){
      # define pacman circles
      shape::filledcircle(r1 = r1,
                   r2 = 0.00001,
                   mid = c(x1[i], y1[i]),
                   from = angles_mouth[i],
                   to = 2 * pi - angles_mouth[i],
                   col = "yellow")
      shape::filledcircle(r1 = r2,
                   r2 = 0.00001,
                   mid = c(x2[i], y2[i]),
                   from = angles_mouth[i],
                   to = 2 * pi - angles_mouth[i],
                   col = "yellow")
      shape::filledcircle(r1 = r4,
                   r2 = 0.00001,
                   mid = c(x4[i], y4[i]),
                   from = angles_mouth[i] + 3,
                   to = 2 * pi - angles_mouth[i] + 3,
                   col = "yellow")

      # dinfine eyes for pacman
      points(x1[i] + 0.2, y1[i] + 0.75, pch = 21, bg = 1, cex = 0.7)
      points(x2[i] + 0.2, y2[i] + 0.75, pch = 21, bg = 1, cex = 0.7)
      points(x4[i] - 0.05, y4[i] + 0.75, pch = 21, bg = 1, cex = 0.7)

      Sys.sleep(dt)

      shape::plotcircle(r = 1.1 * r1,
                 mid = c(x1[i], y1[i]),
                 col = "white",
                 lcol = "white")
      shape::plotcircle(r = 1.1 * r2,
                 mid = c(x2[i], y2[i]),
                 col = "white",
                 lcol = "white")
      shape::plotcircle(r = 1.1 * r4,
                 mid = c(x4[i], y4[i]),
                 col = "white",
                 lcol = "white")
    }
  } else if(type == 3) {
    ## calculate display ratio
    f <- par()$pin[2] / par()$pin[1]

    ## create new overlay plot
    plot(NA,
         xlim = c(0, 100),
         ylim = c(0, 100),
         axes = F,
         frame.plot = FALSE,
         xlab = "",
         ylab = "")

    ## create semi-transparent layer
    polygon(x = c(-100, -100, 200, 200),
            y = c(-100, 200, 200, -100),
            col = rgb(1,1,1, 0.8),
            lty = 0)

    ## draw christmas trees
    n = n.tree
    tree.x <- runif(n, 10, 90)
    tree.y <- runif(n, 10, 90)
    tree.size <- runif(n, 0.3, 1.5)

    for(i in 1:n) {
      ## stem
      polygon(x = c(tree.x[i] - 1.5 * tree.size[i],
                    tree.x[i] - 1.5 * tree.size[i],
                    tree.x[i] + 1.5 * tree.size[i],
                    tree.x[i] + 1.5 * tree.size[i]) ,
              y = c(tree.y[i] - 12 * tree.size[i],
                    tree.y[i] - 1 * tree.size[i],
                    tree.y[i] - 1 * tree.size[i],
                    tree.y[i] - 12* tree.size[i]),
              col = "rosybrown4",
              lty = 0)

      ## branch one
      shape::filledellipse(rx1 = 10 * tree.size[i],
                    rx2 = 0.00001,
                    mid = c(tree.x[i], tree.y[i] + 3 * tree.size[i]),
                    col = "darkgreen",
                    from = 4.0143,
                    to = 5.41052)

      ## branch two
      shape::filledellipse(rx1 = 8 * tree.size[i],
                    rx2 = 0.00001,
                    mid = c(tree.x[i], tree.y[i] + 7 * tree.size[i]),
                    col = "darkgreen",
                    from = 4.0143,
                    to = 5.41052)

      ## branch three
      shape::filledellipse(rx1 = 6 * tree.size[i],
                    rx2 = 0.00001,
                    mid = c(tree.x[i], tree.y[i] + 9 * tree.size[i]),
                    col = "darkgreen",
                    from = 4.0143,
                    to = 5.41052)

      ## branch four
      shape::filledellipse(rx1 = 4 * tree.size[i],
                    rx2 = 0.00001,
                    mid = c(tree.x[i], tree.y[i] + 11 * tree.size[i]),
                    col = "darkgreen",
                    from = 4.0143,
                    to = 5.41052)

      ## sphere one
      shape::filledellipse(rx1 = 1 * f * tree.size[i],
                    ry1 = 1 * tree.size[i],
                    mid = c(tree.x[i] + 2 * tree.size[i],
                            tree.y[i] + 5 * tree.size[i]),
                    col = shape::shadepalette(n = 20, endcol = "darkred"))

      ## sphere two
      shape::filledellipse(rx1 = 0.8 * f * tree.size[i],
                    ry1 = 0.8 * tree.size[i],
                    mid = c(tree.x[i] - 1 * tree.size[i],
                            tree.y[i] + -3 * tree.size[i]),
                    col = shape::shadepalette(n = 20, endcol = "orange"))

      ## sphere three
      shape::filledellipse(rx1 = 1.2 * f * tree.size[i],
                    ry1 = 1.2 * tree.size[i],
                    mid = c(tree.x[i] - 1.7 * tree.size[i],
                            tree.y[i] + 2 * tree.size[i]),
                    col = shape::shadepalette(n = 20, endcol = "yellow3"))

      ## sphere four
      shape::filledellipse(rx1 = 1 * f * tree.size[i],
                    ry1 = 1 * tree.size[i],
                    mid = c(tree.x[i] + 3 * tree.size[i],
                            tree.y[i] - 4 * tree.size[i]),
                    col = shape::shadepalette(n = 20, endcol = "darkblue"))

      Sys.sleep(0.1)
    }

    ## add snow
    points(runif(300, 0, 100), runif(300, 0, 100), pch = 8, col = "lightgrey")
  }
}#end function
