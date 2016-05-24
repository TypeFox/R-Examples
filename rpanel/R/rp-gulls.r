rp.gulls <- function(df.name = "", panel.plot = TRUE) {

click.capture <- function(panel, x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  distance <- sqrt((x - panel$lmks.x)^2 + (y - panel$lmks.y)^2)
  if (min(distance) <= 25) {
    lmk <- which(distance == min(distance))
    if (is.na(panel$lmk1)) {
      panel$lmk1 <- lmk
      rp.line(panel, gulls.image, 
              panel$lmks.x[panel$lmk1], panel$lmks.y[panel$lmk1], 
              panel$lmks.x[panel$lmk1] + 2, 
              panel$lmks.y[panel$lmk1] + 2, 
              width = 4, id = "current0", color = "red")
    }
    else if (is.na(panel$lmk2)) {
      panel$lmk2 <- lmk
      rp.line(panel, gulls.image, 
              panel$lmks.x[panel$lmk1], panel$lmks.y[panel$lmk1], 
              panel$lmks.x[panel$lmk2], panel$lmks.y[panel$lmk2], 
              width=4, color = "red", id = "current")
    }
    else {
      rp.clearlines(panel, gulls.image)
      if ("collected.lmks" %in% names(panel)) 
        for (i in 1:nrow(panel$collected.lmks))
          rp.line(panel, gulls.image, 
                  panel$lmks.x[panel$collected.lmks[i, 1]], 
                  panel$lmks.y[panel$collected.lmks[i, 1]], 
                  panel$lmks.x[panel$collected.lmks[i, 2]], 
                  panel$lmks.y[panel$collected.lmks[i, 2]], 
                  id = "current", width = 4, color = "green")
      panel$lmk1 <- lmk
      rp.line(panel, gulls.image, 
              panel$lmks.x[panel$lmk1], panel$lmks.y[panel$lmk1], 
              panel$lmks.x[panel$lmk1] + 2, 
              panel$lmks.y[panel$lmk1] + 2, 
              width = 4, id = "current0", color = "red")
      panel$lmk2 <- NA
    }
  }
  panel
}
   
collect.data <- function(panel) {
  lmk1 <- panel$lmk1
  lmk2 <- panel$lmk2
  if (is.na(lmk1) | is.na(lmk2)) 
    rp.messagebox("Two landmarks must be defined.")
  else if (lmk1 == lmk2) 
    rp.messagebox("The two landmarks must be different.")
  else {
    lmks <- sort(c(lmk1, lmk2))
    if      (all(lmks == c(1, 3))) vname <- "Wing.Length"
    else if (all(lmks == c(4, 6))) vname <- "Bill.Depth"
    else if (all(lmks == c(5, 8))) vname <- "Head.and.Bill.Length"
    else                           vname <- "none"
    if (vname == "none") {
      if (2 %in% lmks) 
        rp.messagebox(paste(
                "There will be considerable variability in the position", "\n",
                "of the tail feathers, so this is not a good measurement."))
      else if ((lmks[1] %in% 1:3) & (lmks[2] %in% 4:8))
        rp.messagebox(paste(
                "This will depend on the neck position, so it is not a ", "\n",
                "good measurement."))
      else {
        rp.messagebox(paste(
                "This measurement is difficult to make and is not", "\n",
                "easily reproduced."))
      }
      rp.clearlines(panel, gulls.image)
    }
    else {
      y <- panel$gulls.all[, vname]
      added <- FALSE
      if (!("gulls" %in% names(panel))) {
        panel$gulls        <- data.frame(y)
        names(panel$gulls) <- vname
        added              <- TRUE
      }
      else if (vname %in% names(panel$gulls)) {
        rp.messagebox(paste(vname, "is already in the dataset."))
        rp.deleteline(panel, gulls.image, "current")
        added <- FALSE
      }
      else {
        panel$gulls[, vname] <- y
        added <- TRUE
      }
      if (added) {
        rp.messagebox(paste(vname, "has been added to the dataset."))
        gulls <- panel$gulls
        # assign(panel$df.name, gulls, envir = .GlobalEnv)
        if (nchar(panel$df.name) > 0) save(gulls, file = panel$df.name)
        if (panel$panel.plot) {
          nvar <- ncol(panel$gulls)
          if (nvar == 1)
            rp.checkbox(panel, var1, labels = vname, grid = "middle", row = 0, column = 0)
          if (nvar == 2)
            rp.checkbox(panel, var2, labels = vname, grid = "middle", row = 1, column = 0)
          if (nvar == 3)
            rp.checkbox(panel, var3, labels = vname, grid = "middle", row = 2, column = 0)
          if (!("collected.lmks" %in% names(panel))) {
            panel$collected.lmks <- matrix(lmks, ncol = 2)
            rownames(panel$collected.lmks) <- vname
            
            #This function is nested inside another action function, 
            #again to keep its definition in the right place.
#############################################################################            
rp.gulls.plot <- function(panel) {
  with(panel, {
    var.vec <- numeric(0)
    if (var1) var.vec <- 1
    if (var2) var.vec <- c(var.vec, 2)
    if (var3) var.vec <- c(var.vec, 3)
    nvars <- length(var.vec)
    sex <- gulls.all[,1]
    txt <- rep("M", length(sex))
    txt[sex == 0] <- "F"
    clr <- rep("blue", length(sex))
    clr[sex == 0] <- "pink"
    if (nvars == 1) boxplot(gulls[,var.vec] ~ sex,
                        ylab = names(gulls)[var.vec],
                        names = c("Females", "Males"))
    if (nvars == 2) {
      plot(gulls[,var.vec], type = "n")
      text(gulls[,var.vec], txt, col = clr)
    }
    if (nvars == 3) {
       if (require(rgl)) {
          rp.plot3d(
          gulls[,var.vec[2]], gulls[,var.vec[1]], gulls[,var.vec[3]], 
          xlab  = names(gulls)[var.vec[2]], ylab = names(gulls)[var.vec[1]],
          zlab = names(gulls)[var.vec[3]], col = clr)
          } 
       else {
       	  rp.messagebox("You can only plot three variables if package RGL is installed.")
       	  warning("Package RGL is not installed.")
       	  }        
       }
   })
 panel
}
#############################################################################            
            rp.button(panel, title = "Plot data", action = rp.gulls.plot,
                      grid = "right", row = 0, column = 2)
          }
          else {
            panel$collected.lmks <- rbind(panel$collected.lmks, lmks)
            nvars <- nrow(panel$collected.lmks)
            rownames(panel$collected.lmks)[nvars] <- vname
          }
        }
      }
    }
  }
  rp.clearlines(panel, gulls.image)
  if ("collected.lmks" %in% names(panel))
    for (i in 1:nrow(panel$collected.lmks))
      rp.line(panel, gulls.image, 
              panel$lmks.x[panel$collected.lmks[i, 1]], 
              panel$lmks.y[panel$collected.lmks[i, 1]], 
              panel$lmks.x[panel$collected.lmks[i, 2]], 
              panel$lmks.y[panel$collected.lmks[i, 2]], 
              width=4, color = "green")   
# The following code is neater by deleting the current line only,
# rp.deleteline(panel, gulls.image, "current")
# rp.line(panel, gulls.image, 
#                        panel$lmks.x[lmk1], panel$lmks.y[lmk1], 
#                        panel$lmks.x[lmk2], panel$lmks.y[lmk2], 
#                        color = "green", width=4)
  panel$lmk1 <- NA
  panel$lmk2 <- NA
  panel
}
  
  # data("gulls", package = "rpanel")
  image.file <- file.path(system.file(package = "rpanel"), "images", "gulllmks.gif")
  gulls.panel <- rp.control("STEPS:  The Birds and the Bees",
                   gulls.all = gulls.all,
                   lmk.names = c("Wing tip", "Tail feathers", 
                                 "Wing joint", "Bottom of bill",
                                 "Tip of bill", "Top of bill", 
                                 "Top of head", "Back of head"),
                   lmks.x = c(401, 386, 268,  94, 69, 92, 165, 200),
                   lmks.y = c(277, 298, 251,  96, 86, 61,  42,  68),
                   lmk1 = NA, lmk2 = NA, df.name = df.name,
                   var1 = FALSE, var2 = FALSE, var3 = FALSE,
                   panel.plot = panel.plot, size = c(480, 550))
  rp.grid(gulls.panel, "image",  row = 0, column = 0, columnspan = 3)
  rp.grid(gulls.panel, "left",   row = 1, column = 0)
  rp.grid(gulls.panel, "middle", row = 1, column = 1)
  rp.grid(gulls.panel, "right",  row = 1, column = 2)
  rp.image(gulls.panel, image.file, gulls.image, action = click.capture,
            grid = "image", row = 0, column = 0)
  rp.button(gulls.panel, title = "Collect data", action = collect.data,
            grid = "left", row = 0, column = 0)
  invisible()
}

