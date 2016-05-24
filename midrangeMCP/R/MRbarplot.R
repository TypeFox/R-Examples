#' Plotting the result of the multiple comparison procedures
#'
#' \code{MRbarplot} creates a bar plot with vertical or horizontal bars
#'     to compare the mean treatments by the tests:
#'     Skott-Knott midrange, Skott-Knott range, Student-Newman-Keuls
#'      and Tukey midrange.
#' @param x An object of the \code{MRtest} function
#' @param MCP Allows choosing the multiple comparison test.
#'     The \emph{defaut} is "all". This option will perform all tests
#'     available in the \code{MRtest} object.
#' @param col A specification for the plotting color.
#'     The \emph{defaut} is \code{heat.colors(10)}.
#' @param horiz a logical value. If \code{FALSE}, the bars are drawn
#'     vertically with the first bar to the left. If \code{TRUE},
#'     the bars are drawn horizontally with the first at the bottom.
#' @param ... Parameters of the \code{\link{barplot}} function
#' @return \code{MRbarplot} return the bar plot of the tests chosen
#'     ("SKM", "SKR", "SNKM" and "TM")
#'     to evaluate the treatment means.
#' @details The \code{MCP} argument allows choosing several tests
#'     of multiple comparisons from the
#'     \code{MRtest} object. For plots in papers, use
#'     \code{col = gray.colors(10)}. For details, see
#'     \code{\link[grDevices]{colors}} function.
#' @examples
#' # Simulated data (completely randomized design)
#'
#' rv <- c(100.08, 105.66, 97.64, 100.11, 102.60, 121.29, 100.80,
#'         99.11, 104.43, 122.18, 119.49, 124.37, 123.19, 134.16,
#'         125.67, 128.88, 148.07, 134.27, 151.53, 127.31)
#'
#' # Treatments
#' treat <- factor(rep(LETTERS[1:5], each = 4))
#'
#' # Anova
#' res     <- aov(rv~treat)
#'
#' # Loading the midrangeMCP package
#' library(midrangeMCP)
#'
#' # Choosing tests
#' results <- MRtest(y = res, trt = "treat", alpha = 0.05,
#'                    main = "Multiple Comparison Procedures",
#'                    MCP = c("SKM", "TM"))
#'
#' MRbarplot(results, MCP = "all") # It will be shown two
#'                              # graphs. First, for the
#'                              # results of \code{'SKM'}
#'                              # and the second for the
#'                              # results of \code{'TM'}.
#'
#' MRbarplot(results, MCP = "SKM") # It will be shown
#'                                  # only the graph
#'                                  # for the result of
#'                                  # \code{'SKM'}
#'
#' # Plot for papers
#' MRbarplot(results, MCP = "all", col = gray.colors(10))
#' @export
MRbarplot <- function (x,
                       MCP = "all",
                       col = heat.colors(10),
                       horiz = FALSE, ...) {
  mcps  <- c("SKM", "SKR", "SNKM", "TM")
  tests <- sort(x$Tests)
  if (all(MCP == "all")) {
    MCP = tests
  }
  if (length(MCP) == 1) {
    par(mfrow = c(1, 1))
  }
  if (length(MCP) == 2) {
    par(mfrow = c(2, 1))
  }
  if (length(MCP) > 2) {
    par(mfrow = c(2,2))
  }
  if (length(tests) > 1){
    writest <- tests[1]
    for (i in 2:length(tests)) {
      writest <- paste(writest, tests[i], sep = " ")
    }
  } else {
    writest <- tests
  }
  MCP   <- sort(MCP)
  namcp <- pmatch(MCP, tests)
  if (any(is.na(namcp))) {
    stop("The choice of the tests in the MCP argument must be in accordance with the tests chosen in the x argument \n Options: ", writest, call. = FALSE)
  }
  nmcps <- sort(pmatch(tests, mcps))
  means <- x$Groups[[1]][[1]]
  names(means) <- rownames(x$Groups[[1]])
  n     <- length(means)
  lim1  <- max(means)/8
  lim2  <- max(means) + mean(means) * 0.8
  lim3  <- max(means) + mean(means) * 0.4

  if (n > 8) {
    par(las=2) # make label text perpendicular to axis
  }

  # Skott-Knott Midrange plot
  if (any(MCP == "SKM")) {
    groups <- x$Groups[[1]][[2]]
    if (horiz) {
      ind    <- barplot (means, horiz = horiz, col = col,
                         # Color of the bar borders
                         border="black",
                         # (color for x and y labels) Times New Roman : italic
                         font.lab=8,
                         # Font (font for axis annotation) TT Times New Roman : plain
                         font.axis = 6,
                         # Range for x label
                         xlim = c(0, lim3), # Range for x label
                         cex.axis = 0.6,
                         cex.names = 0.6,
                        ...)
      grid(nx = NA, ny = NULL)
      barplot (means,
               horiz = horiz,
               col = col,
               border = "black", # Color of the bar borders
               font.lab = 8, # Font (font for x and y labels ) TT Times New Roman : italic
               font.axis = 6, # Font (font for axis annotation ) TT Times New Roman : plain
               xlim = c(0, lim2), # Range for x label
               cex.axis = 0.6,
               cex.names = 0.6,
               add = TRUE, ...)
      title(mcps[1], font.main = 7, line = 0.5) # font.main = 7 - Font TT Times New Roman : bold
    } else {
      ind    <- barplot(means,
                        horiz = horiz,
                        col = col,
                        border = "black", # Color of the bar borders
                        font.lab = 8, # Font TT (color for x and y labels) Times New Roman : italic
                        font.axis = 6, # Font (font for axis annotation) TT Times New Roman : plain
                        ylim = c(0, lim2), # Range for y label
                        cex.axis = 0.6,
                        cex.names = 0.6,
                        ...)
      grid(nx = NA, ny = NULL)
      barplot(means,
              horiz = horiz,
              col = col,
              border = "black", # Color of the bar borders
              font.lab = 8, # Font (font for x and y labels ) TT Times New Roman : italic
              font.axis = 6, # Font (font for axis annotation ) TT Times New Roman : plain
              ylim = c(0, lim2), # Range for y label
              cex.axis = 0.6,
              cex.names = 0.6,
              add = TRUE, ...)
      title(mcps[1], font.main = 7, line = 0.5) # font.main = 7 - Font TT Times New Roman : bold
    }

    for (k in 1:n) {
      if (horiz) {
        text(means[k] + means[k]/12, ind[k], groups[k], srt = 0, vfont = c("serif", "plain"))
      }
      else {
        text(ind[k], means[k] + lim1, groups[k], srt = 30, vfont = c("serif", "plain"))
      }
    }
  }

  # Skott-Knott Range plot
  if (any(MCP == "SKR")){
    cont <- nmcps <= 1
    cont <- length(cont[cont == TRUE])
    groups <- x$Groups[[cont + 1]][[2]]
    if (horiz) {
      ind    <- barplot(means,
                        horiz = horiz,
                        border="black", # Color of the bar borders
                        font.lab=8, # Font TT (color for x and y labels) Times New Roman : italic
                        font.axis = 6, # Font (font for axis annotation) TT Times New Roman : plain
                        xlim = c(0,lim3), # Range for x label
                        cex.axis = 0.6,
                        cex.names = 0.6,
                        ...)
      grid(nx = NA, ny = NULL)
      barplot(means,
              horiz = horiz,
              col = col,
              border = "black", # Color of the bar borders
              font.lab = 8, # Font (font for x and y labels ) TT Times New Roman : italic
              font.axis = 6, # Font (font for axis annotation ) TT Times New Roman : plain
              xlim = c(0, lim3),
              cex.axis = 0.6,
              cex.names = 0.6,
              add = TRUE, ...)
      title(mcps[2], font.main = 7, line = 0.5) # font.main = 7 - Font TT Times New Roman : bold
    } else{
      ind    <- barplot(means,
                        horiz = horiz,
                        border="black", # Color of the bar borders
                        font.lab=8, # Font TT (color for x and y labels) Times New Roman : italic
                        font.axis = 6, # Font (font for axis annotation) TT Times New Roman : plain
                        ylim = c(0,lim2), # Range for y label
                        cex.axis = 0.6,
                        cex.names = 0.6,
                        ...)
      grid(nx = NA, ny = NULL)
      barplot(means,
              horiz = horiz,
              col = col,
              border = "black", # Color of the bar borders
              font.lab = 8, # Font (font for x and y labels ) TT Times New Roman : italic
              font.axis = 6, # Font (font for axis annotation ) TT Times New Roman : plain
              ylim = c(0, lim2), # Range for y label
              cex.axis = 0.6,
              cex.names = 0.6,
              add = TRUE, ...)
      title(mcps[2], font.main = 7, line = 0.5) # font.main = 7 - Font TT Times New Roman : bold
    }
    for (k in 1:n) {
      if (horiz){
        text(means[k] + means[k]/12, ind[k], groups[k], srt = 0, vfont = c("serif", "plain"))
      } else {
        text(ind[k], means[k] + lim1, groups[k], srt = 30, vfont = c("serif", "plain"))
      }
    }
  }

  # SNK Midrange plot
  if (any(MCP == "SNKM")){
    cont <- nmcps <= 2
    cont <- length(cont[cont == TRUE])
    groups <- x$Groups[[cont + 1]][[2]]
    if (horiz) {
      ind    <- barplot(means,
                        horiz = horiz,
                        border = "black", # Color of the bar borders
                        font.lab = 8, # Font TT (color for x and y labels) Times New Roman : italic
                        font.axis = 6, # Font (font for axis annotation) TT Times New Roman : plain
                        xlim = c(0, lim3), # Range for x label
                        font.main = 7, # Font TT Times New Roman : bold
                        cex.axis = 0.6,
                        cex.names = 0.6,
                        ...)
      grid(nx = NA, ny = NULL)
      barplot(means,
              horiz = horiz,
              col = col,
              border = "black", # Color of the bar borders
              font.lab = 8, # Font (font for x and y labels ) TT Times New Roman : italic
              font.axis = 6, # Font (font for axis annotation ) TT Times New Roman : plain
              xlim = c(0, lim3), # Range for x label
              cex.axis = 0.6,
              cex.names = 0.6,
              add = TRUE, ...)
      title(mcps[3], font.main = 7, line = 0.5) # font.main = 7 - Font TT Times New Roman : bold
    } else{
      ind    <- barplot(means,
                        horiz = horiz,
                        border = "black", # Color of the bar borders
                        font.lab = 8, # Font TT (color for x and y labels) Times New Roman : italic
                        font.axis = 6, # Font (font for axis annotation) TT Times New Roman : plain
                        ylim = c(0, lim2), # Range for y label
                        cex.axis = 0.6,
                        cex.names = 0.6,
                        ...)
      grid(nx = NA, ny = NULL)
      barplot(means,
              horiz = horiz,
              col = col,
              border = "black", # Color of the bar borders
              font.lab = 8, # Font (font for x and y labels ) TT Times New Roman : italic
              font.axis = 6, # Font (font for axis annotation ) TT Times New Roman : plain
              ylim = c(0, lim2), # Range for y label
              cex.axis = 0.6,
              cex.names = 0.6,
              add = TRUE, ...)
      title(mcps[3], font.main = 7, line = 0.5) # font.main = 7 - Font TT Times New Roman : bold
    }
    for (k in 1:n) {
      if (horiz){
        text(means[k] + means[k]/12, ind[k], groups[k], srt = 0, vfont = c("serif", "plain"))
      } else {
        text(ind[k], means[k] + lim1, groups[k], srt = 30, vfont = c("serif", "plain"))
      }
    }
  }

  # Tukey Midrange plot
  if (any(MCP == "TM")){
    cont <- nmcps <= 3
    cont <- length(cont[cont == TRUE])
    groups <- x$Groups[[cont + 1]][[2]]
    if (horiz) {
      ind    <- barplot(means,
                        horiz = horiz,
                        col = col,
                        border = "black", # Color of the bar borders
                        font.lab = 8, # Font TT (color for x and y labels) Times New Roman : italic
                        font.axis = 6, # Font (font for axis annotation) TT Times New Roman : plain
                        xlim = c(0,lim3), # Range for x label
                        cex.axis = 0.6,
                        cex.names = 0.6,
                        ...)
      grid(nx = NA, ny = NULL)
      barplot(means,
              horiz = horiz,
              col = col,
              border = "black", # Color of the bar borders
              font.lab = 8, # Font (font for x and y labels ) TT Times New Roman : italic
              font.axis = 6, # Font (font for axis annotation ) TT Times New Roman : plain
              xlim = c(0, lim3), # Range for x label
              cex.axis = 0.6,
              cex.names = 0.6,
              add = TRUE, ...)
      title(mcps[4], font.main = 7, line = 0.5) # font.main = 7 - Font TT Times New Roman : bold
    } else{
      ind    <- barplot(means,
                        horiz = horiz,
                        col = col,
                        border = "black", # Color of the bar borders
                        font.lab = 8, # Font TT (color for x and y labels) Times New Roman : italic
                        font.axis = 6, # Font (font for axis annotation) TT Times New Roman : plain
                        ylim = c(0,lim2), # Range for y label
                        cex.axis = 0.6,
                        cex.names = 0.6,
                        ...)
      grid(nx = NA, ny = NULL)
      barplot(means,
              horiz = horiz,
              col = col,
              border = "black", # Color of the bar borders
              font.lab = 8, # Font (font for x and y labels ) TT Times New Roman : italic
              font.axis = 6, # Font (font for axis annotation ) TT Times New Roman : plain
              ylim = c(0, lim2), # Range for y label
              cex.axis = 0.6,
              cex.names = 0.6,
              add = TRUE, ...)
      title(mcps[4], font.main = 7, line = 0.5) # font.main = 7 - Font TT Times New Roman : bold
    }

    for (k in 1:n) {
      if (horiz){
        text(means[k] + means[k]/12, ind[k], groups[k], srt = 0, vfont = c("serif", "plain"))
      } else {
        text(ind[k], means[k] + lim1, groups[k], srt = 30, vfont = c("serif", "plain"))
      }
    }
  }
  par(las=0)
}
