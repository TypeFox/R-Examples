#' Plot figural analogies.
#'
#' \code{plot_fa} plots the information of figural analogies.
#' This information is previously generated with function
#' \code{\link{build_fa}}.
#'
#' The \code{plot_fa} function should not be assigned to an
#' object, except when willing to use the optional arguments
#' \code{switch.from} and \code{switch.to}. In the latter
#' case, the object name should be the same as the name of the
#' object of class \code{'fa_items'}, and argument
#' \code{which} must designate the isomorph to be affected
#' even if there is only one isomorph available. For example:
#' \preformatted{object <- plot_fa(object, which = 1,
#' switch.from = 1, switch.to = 2)}
#'
#' @param items An object of class \code{'fa_items'} generated with function \code{build_fa}. No default.
#' @param which A numeric vector designating which isomorph(s) to plot. Plot all by default.
#' @param mode A character string designating plot mode "A", "B" or "C". Plot mode "A" by default.
#' @param language A character string designating English ("E"), German ("D") or Spanish ("S") language. Default is "E".
#' @param language.dir A character string designating language for output files. "A" by default selects all languages.
#' @param form.int A character string designating the form from "A" to "D" of the internal main shape, or "R" for random. Default is "A".
#' @param form.ext A character string designating the form from "A" to "D" of the trapezium, or "R" for random. Default is "A".
#' @param size.shape A number designating the size of every shape. Default is 1.
#' @param size.dot A number designating the size of every shape dot. Default is 2.
#' @param size.line A number designating the thickness of every shape. Default is 1.
#' @param size.q A number designating the size of the question mark. Default is 3.5.
#' @param size.word A number designating the size of the verbal options. Default is 1.2.
#' @param info Should rules applied and correct answers be informed? True by default.
#' @param directory Character string designating a file in your PC where to store the isomorphs.
#' @param switch.from Number 'p' designating an option from 1 to 8 to switch with 'q'.
#' @param switch.to Number 'q' designating an option from 1 to 8 to switch with 'p'.
#' @return A data frame containing rules applied and right answers when \code{info=T} by default, or an object of class \code{'fa_items'} when \code{which} has length 1, its value is greater than 0 and both \code{switch.from} and \code{switch.to} are greater than 0.
#' @author Diego Blum \email{<blumworx@gmail.com>}
#' @references Blum, D., Holling, H., Galibert, M.S., & Forthmann, B. (2016). Task difficulty prediction of figural analogies. \emph{Intelligence, 56,} 72-81.
#' @seealso \code{\link{IMak}}
#' @examples
#' ## Create two isomorphs with one rule, setting the correct answer to 1:
#' one <- build_fa(isomorphs = 2, dot.mov = c(1, 2), correct = 1)
#' ## Plot them:
#' plot_fa(one)
#' ## Change the correct answer of item 2 from position 1 to position 2:
#' one <- plot_fa(one, which = 2, switch.from = 1, switch.to = 2)
#' ## Choose a directory and save the items:
#' # dir1 <- "enter your new directory here"
#' # plot_fa(one, directory = dir1)
#'
#' ## Create four isomorphs with two rules:
#' two <- build_fa(isomorphs = 4, mirror = 1, trap.rot = c(90, 45))
#' ## Plot them in German language:
#' plot_fa(two, language = "D")
#' ## Plot only items 2 and 3 in Spanish and choose form "B" for the internal main shape:
#' plot_fa(two, language = "S", form.int = "B", which = c(2, 3))
#' ## Choose a different directory and save these two items by keeping the latter configuration:
#' # dir2 <- "enter your new directory here"
#' # plot_fa(two, which = c(2, 3), language.dir = "S", form.int = "B", directory = dir2)
#'
#' ## Create 20 isomorphs with three rules. Set automatic=F and affect the options:
#' three <- build_fa(isomorphs = 20, mirror = 1, trap.rot = c(90, 45), dot.mov = c(1, 2),
#' automatic = FALSE, al.mirror = c(0, 1), al.trap.rot = -45, al.dot.mov = 1)
#' ## Plot them:
#' plot_fa(three)
#' ## Plot each individual shape of item 13 in German language only:
#' plot_fa(three, which = 13, mode = "C", language = "D")
#' ## Choose a different directory and save the item parts:
#' # dir3 <- "enter your new directory here"
#' # plot_fa(three, which = 13, mode = "C", language.dir = "D", directory = dir3)
#' @importFrom "grDevices" "dev.off" "png"
#' @importFrom "graphics" "layout" "lines" "par" "plot" "points" "text"
#' @importFrom "utils" "write.table"
#' @export
plot_fa <- function(
  items,
  which = 0,
  mode = "A",
  language = "E",
  language.dir = "A",
  form.int = "A",
  form.ext = "A",
  size.shape = 1,
  size.dot = 2,
  size.line = 1,
  size.q = 3.5,
  size.word = 1.2,
  info = T,
  directory = F,
  switch.from = 0,
  switch.to = 0)

{

  #Forbidden:
  if (!mode %in% c("A", "B", "C") | !language %in% c("E", "D", "S") |
        !language.dir %in% c("A", "E", "D", "S") | !form.int %in% c("A", "B", "C", "D", "R") |
        !form.ext %in% c("A", "B", "C", "D", "R") | !info %in% c(0, 1) |
        is.numeric(size.shape) == F | is.numeric(size.dot) == F | is.numeric(size.line) == F |
        is.numeric(size.q) == F | is.numeric(size.word) == F)
    stop("Incorrect data input.")
  if (!class(items) == "fa_items")
    stop("Object of class 'fa_items' not found.")
  if (switch.from > 0 & switch.to > 0 & (length(which) != 1 | sum(which) == 0))
    stop("Item to switch options is not correctly designated.")
  if ((switch.from > 0 | switch.to > 0) & (!switch.from %in% 1:8 | !switch.to %in% 1:8))
    stop("Incorrect data input for argument 'switch.from' or 'switch.to'.")

  # BASIC CODE FOR FEATURES

  # Rotation function:
  rot <- function(x, angle=0) {
    cbind(cos(angle*0.0174532925)*x[,1] - sin(angle*0.0174532925)*x[,2],
    sin(angle*0.0174532925)*x[,1] + cos(angle*0.0174532925)*x[,2])
  }

  # Reflection function:
  mir <- function(x) {
    cbind(-1*x[,1], x[,2])
  }

  # Broken circle coordinates:
  t <- seq(0, 2*pi, length=200)
  bcirc <- cbind(3.5*sin(t[25:200]), 3.5*cos(t[25:200]))

  # Number of isomorphs:
  isomorphs <- length(items)

  # a, b, c vector:
  abc <- c("a", "b", "c")

  # Items to plot are all or some?
  all_or_some <- 0
  if (sum(which %in% 0:isomorphs) == length(which) &
        (sum(which > 0) == 0 | sum(which > 0) == length(which))) {
    all_or_some <- if (sum(which == 0) == length(which))
      "all" else if (sum(which > 0) == length(which))
        "some"
  } else stop("Incorrect data imput for 'which'.")

  # Creating a data frame with item characteristics:
  information <- data.frame(NA, NA, NA, NA, NA, NA)
  names(information) <- c("Isomorphs", "Correct", "      ",
                          "General rules", "      ", "Total")
  information[,6] <- if (all_or_some == "all")
    isomorphs else length(which)
  information[,4] <- as.character()
  rulenames <- 0
  for (i in 1:length(items[[1]][[12]][[2]])) {
    rulenames[i] <- switch(items[[1]][[12]][[2]][i],
                           "1" = "Main shape rotation",
                           "2" = "Reflection",
                           "3" = "Trapezium rotation",
                           "4" = "Subtraction",
                           "5" = "Dot edge movement")
  }
  for (i in 1:length(rulenames)) {
    information[i,4] <- rulenames[i]
  }
  correct_ones <- 0
  correct_ones <- as.data.frame(correct_ones, ncol=2)
  for (i in 1:isomorphs) {
    correct_ones[i,1] <- i
    correct_ones[i,2] <- paste(items[[i]][[12]][1])
  }
  if (all_or_some == "some") {
    correct_ones <- correct_ones[which,]
  }
  for (i in 1:length(correct_ones[,1])) {
    information[i,1:2] <- ""
    information[i,1] <- correct_ones[i,1]
    information[i,2] <- correct_ones[i,2]
  }
  for (i in 1:length(information[,1])) {
    information[i,3] <- ""
    information[i,5] <- ""
  }
  if (length(rulenames) < length(information[,1]))
    for (i in (length(rulenames) + 1):length(information[,1])) {
      information[i,4] <- ""
    }
  if (length(rulenames) > length(correct_ones[,1]))
    for (i in (length(correct_ones[,1]) + 1):length(rulenames)) {
      information[i,1] <- ""
      information[i,2] <- ""
    }
  if (length(rulenames) != 1 | length(correct_ones[,1]) != 1)
    for (i in 2:length(information[,1])) {
      information[i,6] <- ""
    }

  # PLOTTING FUNCTION

  # Function to plot figures (it will be assigned to a list):
  plot.figure <- function(x) {

    # Identify reflection:
    bcirc_mir_or_not <- if (x$mirror == 0) bcirc else mir(bcirc)
    boot_mir_or_not <- if (x$mirror == 0) boot else mir(boot)

    # Once reflection is identified, apply rotation:
    rot_bcirc <- rot(bcirc_mir_or_not, angle=x$rotation)
    rot_boot <- rot(boot_mir_or_not, angle=x$rotation)
    rot_hammer <- rot(hammer, angle=x$hampos)

    # How much to add when plotting in a different way?
    increment <- c(0, 0, 0)

    # Plot parameters:
    if(directory != F & mode == "C") increment[1] <- 0.5
    if(directory == F & mode == "C") increment[1] <- 0.5
    limits <- 6-size.shape-increment[1]
    par(mai = c(0, 0, 0, 0), pty="s")
    plot(xlim=c(-(limits), limits), ylim=c(-(limits), limits),
         0, type="l", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

    # Drawing the plot:
    if(directory != F) {
      increment[2] <- 10
      increment[3] <- 23
    }
    if(directory != F & mode == "C") {
      increment[2] <- 13
      increment[3] <- 17
    }
    if(directory == F & mode == "C") {
      increment[2] <- 5.5
      increment[3] <- 4.5
    }
    lines(rot_bcirc, lwd=size.line + increment[2])
    for (i in 1:5) {
      if (x$bootlines[i] == 1) lines(rot_boot[i:(i + 1),], lwd=size.line + increment[2])
    }
    point <- cbind(rot_boot[x$dotpos,][1], rot_boot[x$dotpos,][2])
    points(point, pch=16, cex=size.dot + increment[3])
    lines(rot_hammer, lwd=size.line + increment[2])}

  # Setting the number of items to plot:
  items_to_plot <- 0
  if (all_or_some == "all")
    items_to_plot <- 1:isomorphs else
    items_to_plot <- which

  # Setting the directory (if any) to save the items and change language:
  if(directory != F) {
    setwd(directory)
    language <- if(language.dir == "A") 1:3 else
                if(language.dir == "E") 1 else
                if(language.dir == "D") 2 else
                if(language.dir == "S") 3
  } else
    language <- if(language == "E") 1 else
                if(language == "D") 2 else
                if(language == "S") 3
  langname <- 0
  for(i in 1:length(language)) {
    langname[i] <- if(language[i] == 1) "E" else
                   if(language[i] == 2) "D" else
                   if(language[i] == 3) "S"
  }
  langname <- as.vector(langname)

  # Properties of every item to be plotted:
  for (m in items_to_plot){

    # Change places of options when choosing a single item:
    if (length(which) == 1 & which[1] > 0 & switch.from %in% 1:8 &
        switch.to %in% 1:8 & directory == F) {
      if (switch.from == items[[which[1]]][[12]][[1]] |
          switch.to == items[[which[1]]][[12]][[1]]) {
          if (information[1,2] == switch.from) {
            items[[which[1]]][[12]][[1]] <- switch.to
            information[1,2] <- switch.to
          } else {
            items[[which[1]]][[12]][[1]] <- switch.from
            information[1,2] <- switch.from
          }
        warning(paste("Correct answer of item ", which, " changed to position ", information[1,2], ".", sep=""))
      }
      items_to_change <- 0
      items_to_change <- as.list(items_to_change)
      items_to_change[[1]] <- items[[which[1]]][[switch.from]]
      items_to_change[[2]] <- items[[which[1]]][[switch.to]]
      items[[which[1]]][[switch.to]] <- items_to_change[[1]]
      items[[which[1]]][[switch.from]] <- items_to_change[[2]]
    }

    for (h in 1:length(language)){

      if(directory != F & mode != "C")
        png(paste("item", m, langname[h], ".png", sep=""),
            width = switch(mode,A=4800,B=6000),
            height = switch(mode, A=3600, B=1500))

      switch (mode,
              A = {layout(matrix(c(0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0,
                                 0, 0, 0, 0, 3, 3, 4, 4, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 0,
                                 0, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 0),
                                 5, 12, byrow=TRUE), heights=c(1.5, 1.5, 1, 1.5, 1.5))},
              B = {layout(matrix(c(1, 2, 0, 5, 6, 7, 8, 9,
                                 3, 4, 0, 10, 11, 12, 13, 14),
                                 2, 8, byrow=TRUE), widths=c(1.5, 1.5, 1, 1.5, 1.5, 1.5, 1.5, 1.5))},
              C = layout(matrix(1)))
      par(mar = rep(1, 4)) # Shrink margins.

      # Correct plotting when plotting each figure separately:
      if (directory != F) {
        if (mode == "C") {
          plot.fa.size.q <- size.q + 30
          plot.fa.size.word <- size.word + 12
        } else {
          plot.fa.size.q <- size.q + 41
          plot.fa.size.word <- size.word + 16
        }
      } else
          if (mode == "C") {
            plot.fa.size.q <- size.q + 7
            plot.fa.size.word <- size.word + 2.5
          } else {
            plot.fa.size.q <- size.q
            plot.fa.size.word <- size.word
          }

      # Width of every single frame when mode is C:
      width.single <- 800

      # Internal boot coordinates subject to change:
      coord.x <- 3.5*sin(t[25])
      coord.y <- 3.5*cos(t[25])
      randomboot <- list(
        A = matrix(byrow=T, ncol=2, c(0, 3.5, 0, .875, -2, -1.125, coord.x - 1.425,
                                  -1.125, coord.x, .3, coord.x, coord.y)),
        B = matrix(byrow=T, ncol=2, c(0, 3.5, 0, .5, -2.2, .5, .5, -2.2, coord.x,
                                  -2.2 + (coord.x - .5), coord.x, coord.y)),
        C = matrix(byrow=T, ncol=2, c(0, 3.5, -2, 1.5, -2, -1.125, coord.x - 1.425,
                                  -1.125, coord.x - 1.425, coord.y - 1.425, coord.x, coord.y)),
        D = matrix(byrow=T, ncol=2, c(0, 3.5, -2, 1.5, -2, -1.125, .5, 1.375,
                                  coord.x, 1.375 - (coord.x - .5), coord.x, coord.y)))
      boot <- if (form.int == "R") randomboot[[sample(1:4, 1)]] else
        switch(form.int,
               "A" = randomboot[[1]],
               "B" = randomboot[[2]],
               "C" = randomboot[[3]],
               "D" = randomboot[[4]])

      # External trapezium (aka: hammer) coordinates subject to change:
      randomhammer <- list(
        A= cbind(c(1, -1, 0, 1, 1), c(3.5, 3.5, 4.5, 4.5, 3.5)),
        B= cbind(c(1, -1, -1, 0, 1), c(3.5, 3.5, 4.5, 4.5, 3.5)),
        C= cbind(c(0, -1, -1, 1, 0), c(3.5, 3.5, 4.5, 4.5, 3.5)),
        D= cbind(c(1, 0, -1, 1, 1), c(3.5, 3.5, 4.5, 4.5, 3.5)))
      hammer <- if (form.ext == "R")
        randomhammer[[sample(1:4, 1)]] else
        switch(form.ext,
               "A" = randomhammer[[1]],
               "B" = randomhammer[[2]],
               "C" = randomhammer[[3]],
               "D" = randomhammer[[4]])

      # Plot A, B & C:
      for (i in 1:3) {
        if(directory != F & mode == "C")
          png(paste("item", m, abc[i], ".png", sep=""),
              width = width.single, height = width.single)
        plot.figure(items[[m]][[i + 8]])
        if(directory != F & mode == "C") dev.off()
      }

      # Plot question mark:
      if(directory != F & mode == "C")
        png(paste("item", m, "d(q).png", sep=""),
            width = width.single, height=width.single)
      par(mai = c(0, 0, 0, 0), pty="s")
      plot(xlim=c(-5, 5), ylim=c(-5, 5), 0, type="l",
           bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
      text(0, 0, "?", cex=plot.fa.size.q)
      if(directory != F & mode == "C") dev.off()

      # Plot options:
      for (i in 1:8){
        if(directory != F & mode == "C")
          png(paste("item", m, "op", i, ".png", sep=""), width = width.single, height=width.single)
        plot.figure(items[[m]][[i]])
        if(directory != F & mode == "C") dev.off()
      }

      # Plot 9th and 10th options:
      op9 <- if(langname[h] == "D") "Keine\nAntwort\nist richtig" else
             if(langname[h] == "E") "No\nanswer\nis correct" else
             if(langname[h] == "S") "Ninguna\nrespuesta\nes correcta"
      if(directory != F & mode == "C")
        png(paste("item", m, "op9", langname[h], ".png", sep=""), width = width.single,height=width.single)
      par(mai = c(0,0,0,0), pty = "s")
      plot(xlim = c(-5, 5), ylim = c(-5, 5),
           0, type = "l", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      text(0, 0, op9, cex=plot.fa.size.word)
      if(directory != F & mode == "C") dev.off()

      op10 <- if(langname[h] == "D") "Ich wei\u00DF\nnicht" else
              if(langname[h] == "E") "I don\u0027t\nknow" else
              if(langname[h] == "S") "No s\u00E9"
      if(directory != F & mode == "C")
        png(paste("item", m, "op10", langname[h], ".png", sep=""), width = width.single, height=width.single)
      par(mai = c(0, 0, 0, 0), pty = "s")
      plot(xlim = c(-5, 5), ylim = c(-5, 5),
           0, type = "l", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      text(0, 0, op10, cex = plot.fa.size.word)
      if(directory != F & mode == "C") dev.off()

      if(directory != F & mode != "C") dev.off()
    }}

  # Information about the items to be plotted:
  if (directory != F & info == T)
    write.table(information, "Info.csv", quote=F, sep=";", row.names=F)

  if (length(which) == 1 & which[1] > 0 &
      switch.from %in% 1:8 & switch.to %in% 1:8 & directory == F)
    return(items) else
      if (info == T) return(information)
}
